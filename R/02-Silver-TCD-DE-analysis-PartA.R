#####################################################################
#   RNA-seq analysis of bronchoalveolar lavage (BAL) samples from   #
#        LTBI patients after bronchoscopic challenge with PPD       #
#                     Silver Challenge TCD                          #
#        --- Differential Expression analysis - Part A ---          #
#####################################################################

# Author: Correia, C.N.
# Last updated on: 07/05/2019

############################################
# 01 Load and/or install required packages #
############################################

library(here)
library(biomaRt)
library(edgeR)
library(limma)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
library(Cairo)
library(extrafont)
library(biobroom)

# Uncomment functions below to install packages in case you don't have them
# Bioconductor packages
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("edgeR")
#BiocManager::install("biomaRt")
#BiocManager::install("biobroom")


# CRAN packages
#install.packages("here")
#install.packages("tidyverse")
#install.packages("Cairo")
#install.packages("extrafont")

##############################################
# 02 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()

# Define variables for subdirectories
countsDir <- here("unstranded_counts")
imgDir <- here("Figures/")
tablesDir <- here("Tables/")

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

#########################################
# 03 Import featureCounts sense counts  #
#########################################

# Create a vector of file names
files <- list.files(path        = countsDir,
                    pattern     = "^C",
                    all.files   = TRUE,
                    full.names  = FALSE,
                    recursive   = FALSE,
                    ignore.case = FALSE)

files
length(files)

# Create a dataframe with raw counts for all samples
rawCounts <- readDGE(path         = countsDir,
                     files        = files,
                     header       = TRUE,
                     comment.char = "#",
                     columns      = c(1, 7))
names(rawCounts)
head(rawCounts$samples)
head(rawCounts$counts)

# Edit sample names
colnames(rawCounts$counts) %<>%
  str_replace("_counts", "")

rownames(rawCounts$samples) %<>%
  str_replace("_counts", "")

# Check data frames
head(rawCounts$samples)
dim(rawCounts$samples)
head(rawCounts$counts)
dim(rawCounts$counts)

#########################################
# 04 Add sample information for DGElist #
#########################################

# Treatment group (avoid using underscores)
rawCounts$samples$group <- rownames(rawCounts$samples)
rawCounts$samples$group %<>%
  str_replace("C.+S(1|3|5|7|9|11)$", "MED") %>%
  str_replace("C.+S(2|4|6|8|10|12)$", "RV") %>%
  factor(levels = c("MED", "RV"))

# Check data frame
rawCounts$samples
levels(rawCounts$samples$group)

#####################
# 05 Create DGElist #
#####################

# Use newly assigned variables to create DGElist
silver_dgelist <- DGEList(counts       = as.data.frame(rawCounts$counts),
                          group        = rawCounts$samples$group,
                          lib.size     = NULL,
                          norm.factors = NULL,
                          remove.zeros = FALSE)

names(silver_dgelist)
dim(silver_dgelist)
head(silver_dgelist$counts)
head(silver_dgelist$samples)
levels(silver_dgelist$samples$group)

# Include addtional experimental information into DGElist
# Subject number
silver_dgelist$samples$subject <- rownames(silver_dgelist$samples)
silver_dgelist$samples$subject %<>%
  str_replace("C.+S(1|2|3|4)$", "P081") %>%
  str_replace("C.+S(5|6|7|8)$", "P054") %>%
  str_replace("C.+S(9|10|11|12)$", "P084") %>%
  factor(levels = c("P081", "P054", "P084"))

# Check data frame
silver_dgelist$samples
levels(silver_dgelist$samples$subject)

# Cell type (avoid using underscores)
silver_dgelist$samples$cell.type <- rownames(silver_dgelist$samples)
silver_dgelist$samples$cell.type %<>%
  str_replace("C.+S(1|2|5|6|9|10)$", "HAF") %>%
  str_replace("C.+S(3|4|7|8|11|12)$", "LAF") %>%
  factor(levels = c("HAF", "LAF"))

# Check data frame
silver_dgelist$samples
levels(silver_dgelist$samples$cell.type)

# Combine cell type and group
silver_dgelist$samples$cell.group <- factor(paste(silver_dgelist$samples$cell.type,
                                                  silver_dgelist$samples$group,
                                                  sep = "."),
                                            levels = c("HAF.MED", "LAF.MED",
                                                       "HAF.RV", "LAF.RV"))

# Check data frame
silver_dgelist$samples
levels(silver_dgelist$samples$cell.group)

# Export sample information
silver_dgelist$samples %>%
  rownames_to_column(var = "sample_name") %>%
  write_csv(file.path(paste0(tablesDir, "/silver-sample-info.csv")),
            col_names = TRUE)

################################################
# 06 Density plot: raw gene counts per library #
################################################

# Tidy DGElist and plot data
silver_dgelist %>%
  tidy() %>%
  ggplot() +
  geom_density(aes(x     = log10(count + 1),
                   group = sample), size = 0.1) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("Density of raw gene counts per sample") +
  ylab("Density of raw gene counts per sample") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_raw


density_raw

# Export image
ggsave("Silver-density_plot_raw_counts.png",
       plot      = density_raw,
       device    = "png",
       limitsize = FALSE,
       dpi       = 300,
       path      = imgDir)

###########################################
# 07 Remove zero and lowly expressed tags #
###########################################

# Filter non expressed tags (all genes that have zero counts in all samples)
silver_no_zeros <- silver_dgelist[rowSums(silver_dgelist$counts) > 0, ]
dim(silver_no_zeros$counts)
head(silver_no_zeros$counts)

# Filter lowly expressed tags, retaining only tags with
# more than 1 count per million in 3 or more libraries
# (3 libraries correspond to either one of the four groups)
silver_filt <- silver_no_zeros[rowSums(cpm(silver_no_zeros) > 1) >= 3, ]
dim(silver_filt$counts)
head(silver_filt$counts)

# Ouptut filtered counts
silver_filt$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "EntrezID") %>%
  write_csv(file.path(paste0(tablesDir, "/silver-sense_filt_counts.csv")),
            col_names = TRUE)

##############################
# 08 Recompute library sizes #
##############################

silver_filt$samples$lib.size <- colSums(silver_filt$counts)
head(silver_filt$samples)
head(silver_dgelist$samples)

###########################################################################
# 09 Calculate normalisation factors using Trimmed Mean of M-values (TMM) #
###########################################################################

# With edgeR, counts are not transformed in any way after
# calculating normalisation factors
silver_norm <- calcNormFactors(silver_filt, method = "TMM")
head(silver_norm$samples)

##########################
# 10 Save R session info #
##########################

devtools::session_info()

#######################
# 11 Save .RData file #
#######################

# Uncomment line below in case this package hasn't been installed yet:
#devtools::install_github("romainfrancois/nothing")

# Detach all loaded packages (except base R packages):
require(nothing, quietly = TRUE)

# Save all environment objects:
save.image(file = "Silver-TCD-DE-Analysis.RData")

######################################
# Proceed to Part 2 of this analysis #
######################################

# File: 03-Silver-TCD-DE-analysis-PartB.R









