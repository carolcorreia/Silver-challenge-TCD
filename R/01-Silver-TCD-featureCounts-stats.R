#####################################################################
#   RNA-seq analysis of bronchoalveolar lavage (BAL) samples from   #
#        LTBI patients after bronchoscopic challenge with PPD       #
#                     Silver Challenge TCD                          #
#                  --- featureCounts stats ---                      #
#####################################################################

# Author: Correia, C.N.
# Last updated on: 22/05/2019

############################################
# 01 Load and/or install required packages #
############################################

library(here)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)

# Uncomment functions below to install packages in case you don't have them
# CRAN packages
#install.packages("here")
#install.packages("tidyverse")

##############################################
# 02 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()

# Define variables for subdirectories
summaryDir <- here("featCounts_summary")
tablesDir <- here("Tables/")

# Set time zone
Sys.setenv(TZ = "Europe/London")

##################################################
# 03 Import and tidy featureCounts summary files #
##################################################

# List files
summ_files <- list.files(summaryDir, full.names = TRUE)

# Import, tidy and export data:
map_dfc(summ_files, ~ read_tsv(.x)) %>%
  dplyr::rename(S = Status) %>%
  dplyr::select(-contains("Status")) %>%
  dplyr::rename_at(vars(contains("/")),
                   funs(str_remove(., "/home/ccorreia/scratch/Silver_Challenge_Seq/STAR-2.7.0e_alignment/CWRU\\d+_S\\d+/"))) %>%
  dplyr::rename_at(vars(contains("_")),
                   funs(str_remove(., "_Aligned.out.bam"))) %>%
  column_to_rownames(var = "S") %>%
  as.matrix() %>%
  t() %>%
  write.csv(file.path(paste0(tablesDir, "/stats-featureCounts.csv")))

##########################
# 04 Save R session info #
##########################

devtools::session_info()

##########################
# Proceed to DE analysis #
##########################

# File: 02-Silver-TCD-DE-analysis-PartA.R



























