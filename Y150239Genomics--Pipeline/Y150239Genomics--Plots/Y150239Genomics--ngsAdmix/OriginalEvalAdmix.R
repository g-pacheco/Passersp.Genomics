### The BEGINNING ~~~~~
##
# ~ Plots Y150239Genomics--evalAdmix | Written by George Pacheco with help from Jose Samaniego.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, ggnewscale)
source("visFuns_NEW.R")


# Loads Q file ~
q <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MalesOnly.MAFfiltered.Pruned.K7.qopt")


# Loads .corres ~
r <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MalesOnly.MAFfiltered.Pruned.K7.corres")


# Loads annotation ~
ids <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MalesOnly.MAFfiltered.Pruned.K7.labels", as.is = TRUE)


# Adds column ids names ~
colnames(ids) <- c("Sample_ID")


# Expands ids by adding Population ~
ids$Population <- ifelse(grepl("FR0", ids$Sample_ID), "Sales",
                  ifelse(grepl("KAZ", ids$Sample_ID), "Chokpak",
                  ifelse(grepl("Lesina", ids$Sample_ID), "Lesina",
                  ifelse(grepl("Crotone", ids$Sample_ID), "Crotone",
                  ifelse(grepl("Guglionesi", ids$Sample_ID), "Guglionesi",
                  ifelse(grepl("PI22NLD0001M", ids$Sample_ID), "Y150239",
                  ifelse(grepl("PDOM2022NLD0077M", ids$Sample_ID), "Meerkerk",
                  ifelse(grepl("PDOM2022NLD0", ids$Sample_ID), "Utrecht", "Error"))))))))


# Orders according to population and plot the NGSadmix reults ~ 
ord <- orderInds(pop = as.vector(ids[, 2]), q = q)


# Plot correlation of residuals ~
png("Allo.K7.png", width = 1000, height = 1000)
plotCorRes(cor_mat = r, pop = as.vector(ids[, 2]), ord = ord, min_z = -.3, max_z = .3, save_objects = TRUE, save_dir = ".",
           title = "Evaluation of 1000G admixture proportions with K=7")
dev.off()
