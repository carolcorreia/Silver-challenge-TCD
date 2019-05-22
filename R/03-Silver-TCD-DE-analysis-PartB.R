#####################################################################
#   RNA-seq analysis of bronchoalveolar lavage (BAL) samples from   #
#        LTBI patients after bronchoscopic challenge with PPD       #
#                     Silver Challenge TCD                          #
#        --- Differential Expression analysis - Part B ---          #
#####################################################################

# Author: Correia, C.N.
# Last updated on: 07/05/2019

############################################
# 12 Load and/or install required packages #
############################################

library(here)
library(edgeR)
library(limma)
library(devtools)
library(tidyverse)
library(stringr)
library(magrittr)
library(biobroom)
library(ggridges)
library(Cairo)
library(extrafont)
library(ggrepel)
library(ggfortify)
library(statmod)

# Uncomment functions below to install packages in case you don't have them

#install.packages("ggridges")
#install.packages("ggrepel")
#install.packages("ggfortify")
#install.packages("statmod")

####################################################
# 13 Working directory, data, fonts, and time zone #
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

##########################################################
# 14 Tidy DGElist for exploratory data analysis plotting #
##########################################################

tidy_silverNorm <- tidy(silver_norm, addSamples = TRUE)

# Turn sample info into factors for plotting
tidy_silverNorm$sample %<>%
  factor()

# Check factors
levels(tidy_silverNorm$sample)
levels(tidy_silverNorm$group)
levels(tidy_silverNorm$subject)
levels(tidy_silverNorm$cell.type)

# Check data frame
tidy_silverNorm

########################################################
# 15 Plot: density of filtered gene counts per library #
########################################################

ggplot(tidy_silverNorm, aes(x = log10(count + 1),
                          y = sample)) +
  scale_y_discrete(limits = rev(levels(tidy_silverNorm$sample))) +
  geom_density_ridges(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("Condition", values = c("#b2b2b2", "#e06377")) +
  theme_bw(base_size = 12, base_family = "Calibri") +
  ggtitle("Density of filtered gene counts per sample") +
  facet_grid(. ~ group) +
  ylab("Sample code") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_norm

density_norm

# Export high quality image
ggsave("silver-density-filt.pdf",
       plot      = density_norm,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")

#########################################
# 16 Plot: principal component analysis #
#########################################

# Plot PCA using log2CPM filtered counts
pca <- prcomp(t(cpm(silver_norm, log = TRUE)),
              scale. = TRUE)

# Get PCA summary
summary_pca <- summary(pca)
summary_pca

# Tidy PCA summary and export it
tidy(summary_pca$importance) %>%
  write_csv(file.path(paste0(tablesDir, "/summary_pca.csv")),
            col_names = TRUE)

# Plot PC1 versus PC2
pc1vspc2 <- autoplot(pca,
                     data = silver_norm$samples,
                     colour = "group",
                     shape = "cell.type",
                     label = TRUE,
                     label.show.legend = FALSE,
                     label.repel = TRUE) +
            theme_bw(base_size = 12, base_family = "Calibri")

pc1vspc2

# Export high quality image
ggsave(filename  = "pc1vspc2.pdf",
       plot      = pc1vspc2,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 5,
       width     = 6,
       units     = "in")

# Plot PC2 versus PC3
pc2vspc3 <- autoplot(pca,
                     x = 2,
                     y = 3,
                     data = silver_norm$samples,
                     colour = "group",
                     shape = "cell.type",
                     label = TRUE,
                     label.show.legend = FALSE,
                     label.repel = TRUE) +
            theme_bw(base_size = 12, base_family = "Calibri")

pc2vspc3

# Export high quality image
ggsave(filename  = "pc2vspc3.pdf",
       plot      = pc2vspc3,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 5,
       width     = 6,
       units     = "in")

###############################
# 17 Create the design matrix #
###############################

# Create a design matrix
design <- model.matrix(~0 + cell.group,
                       data = silver_norm$samples)

# Check design matrix
design

# Rename design matrix columns for simplicity
colnames(design) %<>%
  str_replace("cell.group", "")

# Check design matrix
design
colnames(design)

#####################################################
# 18 Voom transformation and correlation estimation #
#####################################################

# Convert the read counts to log2-cpm, with associated weights,
# ready for linear modelling
silver_voom <- voom(silver_norm, design, plot = TRUE)

# To make comparisons both within and between subjects,
# it is necessary to treat Patient as a random effect
corfit <- duplicateCorrelation(silver_voom,
                               design,
                               block = silver_voom$targets$subject)
corfit$consensus.correlation

##########################
# 19 Save R session info #
##########################

devtools::session_info()

#######################
# 20 Save .RData file #
#######################

# Detach all loaded packages (except base R packages):
require(nothing, quietly = TRUE)

# Save all environment objects:
save.image(file = "Silver-TCD-DE-Analysis.RData")

######################################
# Proceed to Part 3 of this analysis #
######################################

# File: 04-Silver-TCD-DE-analysis-PartC.R










