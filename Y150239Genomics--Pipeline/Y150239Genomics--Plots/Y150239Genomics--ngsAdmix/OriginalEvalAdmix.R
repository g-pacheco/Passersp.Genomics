### The BEGINNING ~~~~~
##
# ~ Plots Y150239Genomics--evalAdmix | Written by George Pacheco with help from Jose Samaniego.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, ggnewscale)
source("visFuns.R")


# Loads Q file ~
q <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.qopt")


# Loads .corres ~
r <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.corres")


# Loads annotation ~
ids <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.labels", as.is = TRUE)


# Adds column ids names ~
colnames(ids) <- c("Sample_ID")


# Expands ids by adding Population ~
ids$Population <- ifelse(grepl("FR0", ids$Sample_ID), "Sales",
                  ifelse(grepl("KAZ", ids$Sample_ID), "Chokpak",
                  ifelse(grepl("Lesina", ids$Sample_ID), "Lesina",
                  ifelse(grepl("Crotone", ids$Sample_ID), "Crotone",
                  ifelse(grepl("Guglionesi", ids$Sample_ID), "Guglionesi",
                  ifelse(grepl("PI22NLD0001M", ids$Sample_ID), "Garderen Region",
                  ifelse(grepl("PD22NLD0146F", ids$Sample_ID), "Garderen Region",
                  ifelse(grepl("PD22NLD0147F", ids$Sample_ID), "Garderen Region",
                  ifelse(grepl("PDOM2022NLD0077M", ids$Sample_ID), "Garderen Region",
                  ifelse(grepl("PDOM2022NLD0", ids$Sample_ID), "Utrecht", "Error"))))))))))


# Orders according to population and plot the NGSadmix reults ~ 
ord <- orderInds(pop = as.vector(ids[, 2]), q = q)


# Plot correlation of residuals ~
png("Auto.K2_NoLimits.png", width = 1000, height = 1000)
plotCorRes(cor_mat = r, pop = as.vector(ids[, 2]), ord = ord,
           title = "Evaluation of 1000G admixture proportions with K = 2")
dev.off()


