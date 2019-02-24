#####################################################################
#   RNA-seq analysis of bronchoalveolar lavage (BAL) samples from   #
#        LTBI patients after bronchoscopic challenge with PPD       #
#                     Silver Challenge TCD                          #
#              --- Linux bioinformatics workflow ---                #
#####################################################################


# Author: Correia, C.N.
# Last updated on: 24/02/2019

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
do more $file | grep "Analysis complete" >> succesful_fastqc.txt
done

wc -l succesful_fastqc.txt

# Deleted all HTML files:
rm -r *.html

# Collect FastQC stats:
mkdir /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp; \
done

for file in \
`find /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> summary_pre-filtering.txt; \
done

for file in \
`find /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Transfer compressed folders to personal laptop via SCP
# and check HTML reports:
scp -r \
ccorreia@rodeo.ucd.ie:/home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp .

# Remove temporary folder and its files:
rm -r tmp

