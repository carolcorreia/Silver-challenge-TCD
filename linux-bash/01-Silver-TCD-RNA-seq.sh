#####################################################################
#   RNA-seq analysis of bronchoalveolar lavage (BAL) samples from   #
#        LTBI patients after bronchoscopic challenge with PPD       #
#                     Silver Challenge TCD                          #
#              --- Linux bioinformatics workflow ---                #
#####################################################################


# Author: Correia, C.N.
# Last updated on: 22/05/2019

##########################################
# Download raw FASTQ files from provider #
##########################################

# Create and enter working directory on Stampede:
mkdir /home/ccorreia/storage/Silver_Challenge_Seq
cd !$

# Download raw data:
scp -r carolina.correia@seqdownload.cwru.edu:/home/carolina.correia/111418_NextSeq/ .

# Check that all files were downloaded:
cd 111418_NextSeq/
ls -l | grep "CWRU.*_R.*\.fastq\.gz" | wc -l # 24 lines
ls -l | grep "CWRU.*_I.*\.fastq\.gz" | wc -l # 24 lines ### ask what are these.
ls -l | grep "Undetermined.*\.fastq\.gz" | wc -l # 16 lines ### ask what are these.

# In the MacBook's command line, download .html reports from Stampede:
scp -r ccorreia@stampede.ucd.ie:/home/ccorreia/storage/Silver_Challenge_Seq/111418_NextSeq/111418_*.html .

########################
# Perform MD5 checksum #
########################

# Enter working directory:
cd /home/ccorreia/storage/Silver_Challenge_Seq/111418_NextSeq

# Perform md5sum check:
md5sum -c md5sum.txt >> \
/home/ccorreia/storage/Silver_Challenge_Seq/111418_NextSeq/md5check_UCD.txt

# Check that all files passed the check:
wc -l md5sum.txt # 64 lines
grep -c 'OK' md5check_UCD.txt # 64 lines

# Modify folder permissions to read and execute only:
cd /home/ccorreia/storage/Silver_Challenge_Seq
chmod -R 500 111418_NextSeq

###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p /home/ccorreia/scratch/Silver_Challenge_Seq/quality_check/pre-filtering
cd !$

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /home/ccorreia/storage/Silver_Challenge_Seq/111418_NextSeq \
-name *.fastq.gz`; \
do echo "fastqc --noextract --nogroup -t 10 \
-o /home/ccorreia/scratch/Silver_Challenge_Seq/quality_check/pre-filtering $file" \
>> fastqc.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 32 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Analysis complete" >> successful_fastqc.txt;
done

wc -l succesful_fastqc.txt # 64 lines

for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt;
done

wc -l failed_fastqc.txt # Zero lines

# Deleted all HTML files:
rm -r *.html

# Collect FastQC stats:
mkdir /home/ccorreia/scratch/Silver_Challenge_Seq/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/ccorreia/scratch/Silver_Challenge_Seq/quality_check/pre-filtering/tmp; \
done

for file in \
`find /home/ccorreia/scratch/Silver_Challenge_Seq/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> summary_pre-filtering.txt; \
done

for file in \
`find c/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Check sequence quality:
grep 'Per base sequence quality' summary_pre-filtering.txt >> seq_quality.txt
wc -l seq_quality.txt # 64 lines
grep PASS seq_quality.txt | wc -l # 63 lines
grep WARN seq_quality.txt | wc -l # 1 line

# Check if sequence contain adapters:
grep 'Adapter Content' summary_pre-filtering.txt >> adapter_content.txt
wc -l adapter_content.txt # 64 lines
grep PASS adapter_content.txt | wc -l # 32 lines
grep WARN adapter_content.txt | wc -l # 32 lines

# Check sequence length:
grep 'Sequence length' basic_stats_pre-filtering.txt >> seq_length.txt
wc -l seq_length.txt # 64 lines
less seq_length.txt # Contains seqs with 76bp (R1 and R2) and 8bp (I1 and I2)
# I think that the files containg I1 or I2 are just the indexes.

# Transfer compressed folders to personal laptop via SCP
# and check HTML reports:
scp -r \
ccorreia@stampede.ucd.ie:/home/ccorreia/scratch/Silver_Challenge_Seq/quality_check/pre-filtering/tmp .

# Remove temporary folder and its files:
rm -r tmp

#################################################
# Download and index the Human reference genome #
#################################################

# Create and enter working directory:
mkdir /workspace/storage/genomes/homosapiens/h38_p12
cd !$

# Download the reference genome from Ensembl Release 95:
wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download the genome annotation from Ensembl Release 95:
wget ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz
gunzip Homo_sapiens.GRCh38.95.gtf.gz

# To index the reference genome the required software is STAR 2.7.0e,
# consult manual/tutorial for details:
# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Create and enter working directory:
mkdir /workspace/storage/genomes/homosapiens/h38_p12/STAR-2.7.0e_index_75
cd !$

# Generate genome index:
nohup STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /workspace/storage/genomes/homosapiens/h38_p12/STAR-2.7.0e_index_75 \
--genomeFastaFiles \
/workspace/storage/genomes/homosapiens/h38_p12/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /workspace/storage/genomes/homosapiens/h38_p12/Homo_sapiens.GRCh38.95.gtf \
--sjdbOverhang 75 --outFileNamePrefix \
/workspace/storage/genomes/homosapiens/h38_p12/STAR-2.7.0e_index_75 &

#########################################################################
# Alignment of FASTQ files against the Human reference genome with STAR #
#########################################################################

# Required software is STAR 2.7.0e, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Create and enter alignment working directory:
mkdir /home/ccorreia/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment
cd !$

# Create a bash script to perform alignment of paired FASTQ files (did not
# use I1, I2, or Undetermined files. That's why instead of 32 lines,
# there are 12 lines):
for file in `find /home/ccorreia/storage/Silver_Challenge_Seq/111418_NextSeq \
-name CWRU*_R1.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(\.fastq\.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1\.fastq\.gz//'`; \
echo "mkdir /home/ccorreia/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment/$sample; \
cd /home/ccorreia/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment/$sample; \
STAR --runThreadN 10 --runMode alignReads --genomeLoad LoadAndRemove \
--genomeDir /workspace/storage/genomes/homosapiens/h38_p12/STAR-2.7.0e_index_75 \
--readFilesIn $file $file2 --readFilesCommand gunzip -c \
--outFileNamePrefix ./${sample}_ --outSAMtype BAM Unsorted" \
>> alignment.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 6 alignment.sh alignment.sh.
for script in `ls alignment.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'finished successfully' alignment.sh.00.nohup # 6 lines
grep -c 'finished successfully' alignment.sh.01.nohup # 6 lines

# Merge all STAR log.final.out files into a single file:
for file in `find $HOME/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment \
-name *Log.final.out`; \
do perl /home/nnalpas/SVN/star_report_opener.pl -report $file; \
done

# Using my mac command line, transfer mapping stats via SCP:
scp \
ccorreia@stampede.ucd.ie:/home/ccorreia/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment/All_star_log_final_out.txt .

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir $HOME/scratch/Silver_Challenge_Seq/quality_check/post_alignment
cd !$

# Create a bash script to perform FastQC quality check on aligned BAM files:
for file in `find $HOME/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/Silver_Challenge_Seq/quality_check/post_alignment $file" >> \
fastqc_aligned.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 1 fastqc_aligned.sh fastqc_aligned.sh.
for script in `ls fastqc_aligned.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Delete all the HTML files:
rm -r *.html

# Check if all files were processed:
for file in `ls *.nohup`; \
do more $file | grep "Analysis complete" >> successful_fastqc.txt;
done

wc -l successful_fastqc.txt # 64 lines

for file in `ls *.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt;
done

wc -l failed_fastqc.txt # Zero lines

# Check all output from FastQC:
mkdir $HOME/scratch/Silver_Challenge_Seq/quality_check/post_alignment/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/Silver_Challenge_Seq/quality_check/post_alignment/tmp; \
done

for file in \
`find $HOME/scratch/Silver_Challenge_Seq/quality_check/post_alignment/tmp \
-name summary.txt`; do more $file >> reports_post-alignment.txt; done

for file in \
`find $HOME/scratch/Silver_Challenge_Seq/quality_check/post_alignment/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check sequence quality:
grep 'Per base sequence quality' reports_post-alignment.txt >> seq_quality.txt
wc -l seq_quality.txt # 12 lines
grep PASS seq_quality.txt | wc -l # 12 lines
grep WARN seq_quality.txt | wc -l # Zero lines

# Check if sequence contain adapters:
grep 'Adapter Content' reports_post-alignment.txt >> adapter_content.txt
wc -l adapter_content.txt # 12 lines
grep PASS adapter_content.txt | wc -l # 12 lines
grep WARN adapter_content.txt | wc -l # Zero lines

# Check sequence length:
grep 'Sequence length' basic_stats_post_alignment.txt >> seq_length.txt
wc -l seq_length.txt # 12 lines
less seq_length.txt # 76bp

# Transfer compressed folders to personal laptop via SCP
# and check HTML reports:
scp -r \
ccorreia@stampede.ucd.ie:/home/ccorreia/scratch/Silver_Challenge_Seq/quality_check/post_alignment/tmp .

# Remove temporary folder and its files:
rm -r tmp

#########################
# Calculate insert size #
#########################

# Create and go to working directory:
mkdir $HOME/scratch/Silver_Challenge_Seq/insert_size
cd !$

# Sort aligned BAM files:
for file in \
`find $HOME/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment \
-name *Aligned.out.bam`; \
do outfile=`basename $file | perl -p -e 's/Aligned.out.bam/Sorted.out.bam/'`; \
echo "java -jar /usr/local/src/picard-tools-1.137/picard.jar \
SortSam I=$file O=$outfile SORT_ORDER=coordinate" >> sort.sh; \
done

# Run script on Stampede:
chmod 755 sort.sh
nohup ./sort.sh > sort.sh.nohup &

# Collect insert sizes:
for file in `ls *_Sorted.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Sorted.out.bam//'`; \
echo "java -jar /usr/local/src/picard-tools-1.137/picard.jar \
CollectInsertSizeMetrics \
I=$file \
O=${sample}_insert_size_metrics.txt \
H=${sample}_insert_size_histogram.pdf M=0.5" >> collect_insert_size.sh; \
done

# Run script on Stampede:
chmod 755 collect_insert_size.sh
nohup ./collect_insert_size.sh > collect_insert_size.sh.nohup &

# Collect insert size metrics for all samples into one file:
for file in `ls $HOME/scratch/Silver_Challenge_Seq/insert_size/*_insert_size_metrics.txt`; \
do sample=`basename $file | perl -p -e 's/_insert_size_metrics.txt//'`; \
header=`grep 'MEDIAN_INSERT_SIZE' $file`; \
stats=`sed -n '/MEDIAN_INSERT_SIZE/{n;p;}' $file`; \
printf "Sample_id\t$header\n$sample\t$stats" \
>> All_insert_size.txt; \
done

# Transfer stats to laptop:
scp \
ccorreia@stampede.ucd.ie:/home/ccorreia/scratch/Silver_Challenge_Seq/insert_size/All_insert_size.txt .

###################################################
# Summarisation of gene counts with featureCounts #
###################################################

# Required package is featureCounts, which is part of Subread 1.6.4 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directory:
mkdir -p $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded
cd !$

# Create a bash script to run featureCounts on BAM files
# using the unstranded parameter:
for file in `find $HOME/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment \
-name *_Aligned.out.bam`; \
do sample=`basename $file | perl -p -e 's/_Aligned.out.bam//'`; \
echo "mkdir $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/$sample; \
cd $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/$sample; \
featureCounts -a \
/workspace/storage/genomes/homosapiens/h38_p12/Homo_sapiens.GRCh38.95.gtf \
-B -p -C -R BAM -s 0 -T 20 -t exon -g gene_id \
-o ${sample}_counts.txt $file" >> counts.sh; \
done

# Run script on Stampede:
chmod 755 counts.sh
nohup ./counts.sh > counts.sh.nohup &

# Check if all files were processed:
grep -c 'Summary of counting results can be found in file' counts.sh.nohup # 12 lines

# Copy all *counts.txt files to temporary folder:
mkdir $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/tmp

for file in `find $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/ \
-name *_counts.txt`; do cp $file \
-t $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/tmp; \
done

# Transfer count files to laptop:
scp -r ccorreia@stampede.ucd.ie:/home/ccorreia/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/tmp .

# Remove tmp folder:
rm -r tmp

# Copy all *summary files to temporary folder:
mkdir $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/tmp2

for file in `find $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/ \
-name *summary`; do cp $file \
-t $HOME/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/tmp2; \
done

# Transfer count files to laptop:
scp -r ccorreia@stampede.ucd.ie:/home/ccorreia/scratch/Silver_Challenge_Seq/Count_summarisation_unstranded/tmp2 .

# Remove tmp folder:
rm -r tmp2

# All featureCounts summary stats will be merged into one file using R: 01-Silver-TCD-featureCounts-stats.R

########################################
# R analysis of gene counts with edgeR #
########################################

# Subsequent DE analysis, performed in R: 02-Silver-TCD-DE-analysis-PartA.R


























