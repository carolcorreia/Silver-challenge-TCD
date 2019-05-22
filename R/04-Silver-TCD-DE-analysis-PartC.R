#####################################################################
#   RNA-seq analysis of bronchoalveolar lavage (BAL) samples from   #
#        LTBI patients after bronchoscopic challenge with PPD       #
#                     Silver Challenge TCD                          #
#        --- Differential Expression analysis - Part C ---          #
#####################################################################

# Author: Correia, C.N.
# Last updated on: 07/05/2019

############################################
# 21 Load and/or install required packages #
############################################

library(here)
library(edgeR)
library(limma)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
library(Cairo)
library(extrafont)
library(VennDiagram)
library(forcats)
library(biomaRt)

# Uncomment functions below to install packages in case you don't have them
# Bioconductor packages
#library(BiocManager)
#BiocManager::install("biomaRt")

# CRAN packages
#install.packages("VennDiagram")

####################################################
# 22 Working directory, data, fonts, and time zone #
####################################################

# Check working directory
here()

# Load previously saved data
load("Silver-TCD-DE-Analysis.RData")

# Check variables for subdirectories
imgDir
tablesDir

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

#######################
# 23 Define contrasts #
#######################

# Specifies which comparisons you would like to make
# between the RNA samples

### HAF.MvsI: Comparison between the infected RV and
# control medium-incubated groups in HAF cells
### LAF.IvsM: Comparison between the infected RV and
# control medium-incubated groups in LAF cells
### LAFvsHAF.I: Comparison between the infected LAF and HAF cells
### LAFvsHAF.m: Comparison between the medium-incubated LAF and HAF cells

contrastsMatrix <- makeContrasts(HAF.IvsM = HAF.RV-HAF.MED,
                                 LAF.IvsM = LAF.RV-LAF.MED,
                                 LAFvsHAF.I = LAF.RV-HAF.RV,
                                 LAFvsHAF.MED = LAF.MED-HAF.MED,
                                 levels = colnames(design))

# Check matrix
contrastsMatrix

# Output contrasts matrix
as.data.frame(contrastsMatrix) %>%
  rownames_to_column(var = "cell.condition") %>%
  write_csv(file.path(paste0(tablesDir, "Contrasts_Matrix.csv")),
            col_names = TRUE)

###########################
# 24 Fit the linear model #
###########################

# Fit linear model
silver_fit <- lmFit(silver_voom,
                    design,
                    block = silver_voom$targets$subject,
                    correlation = corfit$consensus.correlation)

# Compute these contrasts and moderated t-tests
silver_fitC <- contrasts.fit(silver_fit, contrastsMatrix)
silver_ebayes <- eBayes(silver_fitC)

# Plot mean-variance trend of final model
png(filename = "final-model-variance.png",
    width    = 1366,
    height   = 768,
    units    = "px")
plotSA(silver_ebayes)
dev.off()

########################################
# 25 Determine differential expression #
########################################

# Test for differential expression between the different comparisons
# using the groups defined in the contrasts matrix
contrastsMatrix

## 1 ##
### HAF.MvsI: Comparison between the infected RV and
# control medium-incubated groups in HAF cells
HAF_RV_DE <- topTable(silver_ebayes,
                      coef   = "HAF.IvsM",
                      number = "inf",
                      adjust.method = "BH",
                      confint = TRUE)

head(HAF_RV_DE)
dim(HAF_RV_DE)

## 2 ##
### LAF.IvsM: Comparison between the infected RV and
# control medium-incubated groups in LAF cells
LAF_RV_DE <- topTable(silver_ebayes,
                      coef   = "LAF.IvsM",
                      number = "inf",
                      adjust.method = "BH",
                      confint = TRUE)

head(LAF_RV_DE)
dim(LAF_RV_DE)

## 3 ##
### LAFvsHAF.I: Comparison between the infected LAF and HAF cells
DE_LAFvsHAF_RV <- topTable(silver_ebayes,
                           coef   = "LAFvsHAF.I",
                           number = "inf",
                           adjust.method = "BH",
                           confint = TRUE)

head(DE_LAFvsHAF_RV)
dim(DE_LAFvsHAF_RV)

## 4 ##
### LAFvsHAF.m: Comparison between the medium-incubated LAF and HAF cells
DE_LAFvsHAF_MED <- topTable(silver_ebayes,
                            coef   = "LAFvsHAF.MED",
                            number = "inf",
                            adjust.method = "BH",
                            confint = TRUE)

head(DE_LAFvsHAF_MED)
dim(DE_LAFvsHAF_MED)

#######################
# 26 Get gene symbols #
#######################

# Check list of the Ensembl mart available for release 95
listEnsembl(version = 95)
ensembl95 = useEnsembl(biomart = "ensembl", version = 95)
listDatasets(ensembl95)[55:60,]

# Connect to the Ensembl gene mart human dataset for release 95
ensembl95 = useEnsembl(biomart = "ensembl",
                       dataset = "hsapiens_gene_ensembl",
                       version = 95)

listAttributes(ensembl95)[57:60,1]

# Prepare list of IDs to query
ensembl_ids <- rownames(HAF_RV_DE)

# Get gene symbols and description
gene_annotation <- getBM(attributes = c("ensembl_gene_id",
                                        "hgnc_symbol",
                                        "description"),
                         filters    = "ensembl_gene_id",
                         values     = ensembl_ids,
                         mart       = ensembl95)

# Check data frame
head(gene_annotation)
dim(gene_annotation)
length(ensembl_ids)

# Get duplicated ID
gene_annotation$ensembl_gene_id[duplicated(gene_annotation$ensembl_gene_id)]
dplyr::filter(gene_annotation, ensembl_gene_id == "ENSG00000276085")

# Remove gene CCL3L3, since CCL3L1 shows up as the main gene in the
# Ensembl website when searching with the ID
gene_annotation %<>%
  dplyr::filter(!hgnc_symbol == "CCL3L3")

# Check data frame
dim(gene_annotation)

#####################################
# 27 Output all genes tested for DE #
#####################################

## 1 ##
# HAF/AM RV x HAF/AM MED
HAF_RV_DE %>%
  rownames_to_column(var = "Ensembl_ID") %>%
    dplyr::full_join(y  = gene_annotation,
                     by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(Ensembl_ID, hgnc_symbol, description, everything()) %>%
  dplyr::mutate(`is DE (adj.P.Val < 0.1)` = if_else(adj.P.Val < 0.1,
                                               "TRUE",
                                               "FALSE")) %>%
  dplyr::arrange(adj.P.Val) %>%
  write_csv(file.path(paste0(tablesDir, "HAF_RV_AllGenes.csv")),
            col_names = TRUE)

## 2 ##
# LAF/AM RV x LAF/AM MED
LAF_RV_DE %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::full_join(y  = gene_annotation,
                   by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(Ensembl_ID, hgnc_symbol, description, everything()) %>%
  dplyr::mutate(`is DE (adj.P.Val < 0.1)` = if_else(adj.P.Val < 0.1,
                                               "TRUE",
                                               "FALSE")) %>%
  dplyr::arrange(adj.P.Val) %>%
  write_csv(file.path(paste0(tablesDir, "LAF_RV_AllGenes.csv")),
            col_names = TRUE)


## 3 ##
# LAF/AM RV x HAF/AM RV
DE_LAFvsHAF_RV %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::full_join(y  = gene_annotation,
                   by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(Ensembl_ID, hgnc_symbol, description, everything()) %>%
  dplyr::mutate(`is DE (adj.P.Val < 0.1)` = if_else(adj.P.Val < 0.1,
                                               "TRUE",
                                               "FALSE")) %>%
  dplyr::arrange(adj.P.Val) %>%
  write_csv(file.path(paste0(tablesDir, "DE_LAFvsHAF_RV_AllGenes.csv")),
            col_names = TRUE)

## 4 ##
# LAF/AM MED x HAF/AM MED
DE_LAFvsHAF_MED %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::full_join(y  = gene_annotation,
                   by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(Ensembl_ID, hgnc_symbol, description, everything()) %>%
  dplyr::mutate(`is DE (adj.P.Val < 0.1)` = if_else(adj.P.Val < 0.1,
                                               "TRUE",
                                               "FALSE")) %>%
  dplyr::arrange(adj.P.Val) %>%
  write_csv(file.path(paste0(tablesDir, "DE_LAFvsHAF_MED_AllGenes.csv")),
            col_names = TRUE)

#############################################
# 28 Filter genes considered DE (FDR < 0.1) #
#############################################

## 1 ##
# HAF/AM RV x HAF/AM MED
HAF_RV_DE %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::full_join(y  = gene_annotation,
                   by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(Ensembl_ID, hgnc_symbol, description, everything()) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::arrange(adj.P.Val) %>%
  as_tibble()

## 2 ##
# LAF/AM MED x LAF/AM RV
LAF_RV_DE %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::full_join(y  = gene_annotation,
                   by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(Ensembl_ID, hgnc_symbol, description, everything()) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::arrange(adj.P.Val) %>%
  as_tibble()

## 3 ##
# LAF/AM RV x HAF/AM RV
DE_LAFvsHAF_RV %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::full_join(y  = gene_annotation,
                   by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(Ensembl_ID, hgnc_symbol, description, everything()) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::arrange(adj.P.Val) %>%
  as_tibble() -> DE_LAFvsHAF_RV_FDR_01

## 4 ##
# LAF/AM MED x HAF/AM MED
DE_LAFvsHAF_MED %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::full_join(y  = gene_annotation,
                   by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(Ensembl_ID, hgnc_symbol, description, everything()) %>%
  dplyr::filter(adj.P.Val < 0.1) %>%
  dplyr::arrange(adj.P.Val) %>%
  as_tibble() -> DE_LAFvsHAF_MED_FDR_01

### Output all lists of DE genes
DElists <- list(DE_LAFvsHAF_RV_FDR_01, DE_LAFvsHAF_MED_FDR_01)
DEfiles <- c(paste0(c("DE_LAFvsHAF_RV_FDR_01", "DE_LAFvsHAF_MED_FDR_01"),
                    "_genes.csv"))
DEpaths <- file.path(tablesDir, DEfiles)

pwalk(list(DElists, DEpaths),
      write_csv,
      col_names = TRUE)

##########################
# 29 Save R session info #
##########################

devtools::session_info()

#######################
# 30 Save .RData file #
#######################

# Detach all loaded packages (except base R packages):
require(nothing, quietly = TRUE)

# Save all environment objects:
save.image(file = "Silver-TCD-DE-Analysis.RData")


