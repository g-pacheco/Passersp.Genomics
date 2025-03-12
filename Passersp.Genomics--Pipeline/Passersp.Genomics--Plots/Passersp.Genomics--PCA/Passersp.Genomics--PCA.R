### The BEGINNING ~~~~~
##
# Plots Passer sp. Genomics -- PCA | Written by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(optparse, tidyverse, ggh4x, plyr, RColorBrewer, extrafont, ggforce, ggstar, ggrepel, RcppCNPy, reshape2,
               gridExtra, grid, cowplot, patchwork, ggpubr, rphylopic)


# Loads data ~
data.auto <- as.matrix(read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.PCAone.SVD0.cov"),
                      header = FALSE, stringsAsFactors = FALSE)
data.allo <- as.matrix(read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.PCAone.SVD0.cov"),
                      header = FALSE, stringsAsFactors = FALSE)


# Loads annotations files ~
annot.auto <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.labels",
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)
annot.allo <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.labels",
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)


# Runs PCA ~
PCAauto <- eigen(data.auto)
PCAallo <- eigen(data.allo)


# Merges the first 3 PCs with annot ~
PCAauto_Annot <- as.data.frame(cbind(annot.auto, PCAauto$vectors[, c(1:3)])); colnames(PCAauto_Annot) <- c("Sample_ID", "PCA_1", "PCA_2", "PCA_3")
PCAallo_Annot <- as.data.frame(cbind(annot.allo, PCAallo$vectors[, c(1:3)])); colnames(PCAallo_Annot) <- c("Sample_ID", "PCA_1", "PCA_2", "PCA_3")


# Merges the first 3 PCs with annot ~
PCAauto_Annot$CHR <- "Autosomes"
PCAallo_Annot$CHR <- "Chromosome Z"


# Binds the 2 DFs based on common columns ~
fulldf <- rbind(PCAauto_Annot, PCAallo_Annot)


# Expands PCA_Annot by adding Population ~
fulldf$Population <- ifelse(grepl("FR0", fulldf$Sample_ID), "Sales",
                     ifelse(grepl("KAZ", fulldf$Sample_ID), "Chokpak",
                     ifelse(grepl("Lesina", fulldf$Sample_ID), "Lesina",
                     ifelse(grepl("Crotone", fulldf$Sample_ID), "Crotone",
                     ifelse(grepl("Guglionesi", fulldf$Sample_ID), "Guglionesi",
                     ifelse(grepl("PI22NLD0001M", fulldf$Sample_ID), NA,
                     ifelse(grepl("PD22NLD0146F", fulldf$Sample_ID), NA,
                     ifelse(grepl("PD22NLD0147F", fulldf$Sample_ID), NA,
                     ifelse(grepl("PDOM2022NLD0077M", fulldf$Sample_ID), NA,
                     ifelse(grepl("PDOM2022NLD0", fulldf$Sample_ID), "Utrecht", "Error"))))))))))


# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                        levels = c("Utrecht",
                                   "Sales",
                                   "Crotone",
                                   "Guglionesi",
                                   "Lesina",
                                   "Chokpak",
                                   NA))


# Expands PCA_Annot by adding Species ~
fulldf$Species <- ifelse(fulldf$Population %in% c("Utrecht", "Sales"), "House",
                  ifelse(fulldf$Population %in% c("Chokpak", "Lesina"), "Spanish",
                  ifelse(fulldf$Population %in% c("Crotone", "Guglionesi"), "Italian",
                  ifelse(fulldf$Population %in% NA, NA, "Error"))))


# Reorders Population ~
fulldf$Species <- factor(fulldf$Species, ordered = T,
                         levels = c("House",
                                    "Italian",
                                    "Spanish",
                                    NA))


# Defines the shapes to be used for each Group ~
shapes.legend <- as.vector(c(1, 9, 13, 21, 11, 23))


# Creates legend plot ~
MyLegend_Plot <-
  ggplot(data = fulldf, aes_string(x = "PCA_1", y = "PCA_2")) +
  geom_star(aes(starshape = Population, fill = Species), size = 2.8, starstroke = .15, alpha = .7) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000"), na.translate = FALSE) +
  scale_starshape_manual(values = shapes.legend, na.translate = FALSE) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.position = "top",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
        legend.box = "vertical",
        legend.box.margin = margin(t = 10, b = 15, r = 0, l = 0)) +
  guides(starshape = guide_legend(title = "Population", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                                  label.theme = element_text(size = 15, family = "Optima"),
                                  override.aes = list(starshape = shapes.legend, size = 5, starstroke = .15), nrow = 1, order = 2),
         fill = guide_legend(title = "Species", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                             label.theme = element_text(size = 15, family = "Optima"),
                             override.aes = list(starshape = 21, size = 5, starstroke = .15), nrow = 1, order = 1),
         colour = "none")


# Defines the shapes to be used for each Group ~
shapes.auto <- as.vector(c(1, 9, 13, 21, 11, 23, 14))


# Combines all populations from the Faroe Islands ~
fulldf$Species <- as.character(fulldf$Species)
fulldf$Population <- as.character(fulldf$Population)


# Expands fulldf by adding Labels ~
fulldf$Species <- ifelse(fulldf$Sample_ID %in% c("PI22NLD0001M", "PD22NLD0146F", "PD22NLD0147F", "PDOM2022NLD0077M"),
                                                 "Focal", fulldf$Species)


# Expands fulldf by adding Labels ~
fulldf$Population <- ifelse(fulldf$Sample_ID %in% c("PI22NLD0001M", "PD22NLD0146F", "PD22NLD0147F", "PDOM2022NLD0077M"),
                                                    "Focal", fulldf$Population)


# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                               levels = c("Utrecht",
                                          "Sales",
                                          "Crotone",
                                          "Guglionesi",
                                          "Lesina",
                                          "Chokpak",
                                          "Focal"))


# Reorders Population ~
fulldf$Species <- factor(fulldf$Species, ordered = T,
                            levels = c("House",
                                       "Italian",
                                       "Spanish",
                                       "Focal"))


# Expands PCA_Annot by adding Labels ~
fulldf$Labels <- ifelse(fulldf$Sample_ID %in% c("PI22NLD0001M"), "Focal Ind.",
                 ifelse(fulldf$Sample_ID %in% c("PD22NLD0146F"), "Garderen_01",
                 ifelse(fulldf$Sample_ID %in% c("PD22NLD0147F"), "Garderen_02",
                 ifelse(fulldf$Sample_ID %in% c("PDOM2022NLD0077M"), "Meerkerk_01", ""))))


# Gets Eigenvalues of each Eigenvectors (Allosome) ~
PCAauto_Eigenval_Sum <- sum(PCAauto$values)
(PCAauto$values[1]/PCAauto_Eigenval_Sum)*100
(PCAauto$values[2]/PCAauto_Eigenval_Sum)*100
(PCAauto$values[3]/PCAauto_Eigenval_Sum)*100


PCAauto_12 <-
  ggplot(data = subset(fulldf, CHR == "Autosomes"), aes_string(x = "PCA_1", y = "PCA_2")) +
  geom_star(aes(starshape = Population, fill = Species), alpha = .7, size = 2.15, starstroke = .15) +
  facet_grid2(CHR ~. , scales = "free_x", axes = "all", remove_labels = "x") +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#d9d9d9")) +
  scale_starshape_manual(values = shapes.auto) +
  geom_label_repel(data = subset(fulldf, CHR == "Autosomes" & Labels == c("Focal Ind.", "Garderen_02")), aes(label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = -.06, nudge_y = 0,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                   ends = "last", type = "open")) +
  geom_label_repel(data = subset(fulldf, CHR == "Autosomes" & Labels == c("Garderen_01", "Meerkerk_01")), aes(label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .04, nudge_y = -.1,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_mark_ellipse(aes(filter = Species == "House", label = "House\nSparrow"), con.colour = "#1E90FF", colour = "#1E90FF",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "straight", label.family = "Optima", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  geom_mark_ellipse(aes(filter = Species == "Spanish", label = "Spanish\nSparrow"), con.colour = "#ee0000", colour = "#ee0000",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "elbow", label.family = "Optima", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  geom_mark_ellipse(aes(filter = Species == "Italian", label = "Italian\nSparrow"), con.colour = "#FFD700", colour = "#FFD700",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "elbow", label.family = "Optima", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  scale_x_continuous("PC 1 (5.55%)",
                     breaks = c(-.1, 0, .1, .2),
                     labels = c("-0.1", "0", ".01", ".02"),
                     limits = c(-.19, .25),
                     expand = c(0, 0)) +
  scale_y_continuous("PC 2 (1.97%)",
                     #breaks = c(-.08, -.04, 0.00), 
                     #labels = c("-0.08", "-0.04", "0.00"),
                     limits = c(-.31, .35),
                     expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(.2, "cm"),
        legend.position = "none",
        axis.title.x = element_text(family = "Optima", size = 16, face = "bold", margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Optima", size = 16, face = "bold", margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.text = element_text(family = "Optima", color = "#000000", size = 11, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


# Gets Eigenvalues of each Eigenvectors (Autosomes) ~
PCAallo_Eigenval_Sum <- sum(PCAallo$values)
(PCAallo$values[1]/PCAallo_Eigenval_Sum)*100
(PCAallo$values[2]/PCAallo_Eigenval_Sum)*100
(PCAallo$values[3]/PCAallo_Eigenval_Sum)*100


# Defines the shapes to be used for each Group ~
shapes.allo <- as.vector(c(1, 9, 13, 21, 11, 23, 14))


PCAallo_12 <-
  ggplot(data =  subset(fulldf, CHR == "Chromosome Z"), aes_string(x = "PCA_1", y = "PCA_2")) +
  geom_star(aes(starshape = Population, fill = Species), alpha = .7, size = 2.15, starstroke = .15) +
  facet_grid2(CHR ~. , scales = "free_x", axes = "all", remove_labels = "x") +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#d9d9d9")) +
  scale_starshape_manual(values = shapes.allo) +
  geom_label_repel(data = subset(fulldf, CHR == "Chromosome Z"), aes(label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = -.03, nudge_y = .015,
                   point.padding = .6, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_mark_ellipse(aes(filter = Species == "House", label = "House\nSparrow"), con.colour = "#1E90FF", colour = "#1E90FF",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "straight", label.family = "Optima", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  geom_mark_ellipse(aes(filter = Species == "Spanish", label = "Spanish\nSparrow"), con.colour = "#ee0000", colour = "#ee0000",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "elbow", label.family = "Optima", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  geom_mark_ellipse(aes(filter = Species == "Italian", label = "Italian\nSparrow"), con.colour = "#FFD700", colour = "#FFD700",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "elbow", label.family = "Optima", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  scale_x_continuous("PC 1 (9.10%)",
                     limits = c(-.19, .25),
                     expand = c(0, 0)) +
  scale_y_continuous("PC 2 (4.80%)",
                     limits = c(-.31, .35),
                     expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(.2, "cm"),
        legend.position = "none",
        axis.title.x = element_text(family = "Optima", size = 16, face = "bold", margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Optima", size = 16, face = "bold", margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.text = element_text(family = "Optima", color = "#000000", size = 11, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


# Isolates legend ~
MyLegendBlog <- get_legend(MyLegend_Plot)


# Gets final plot ~
PCA_Plot <- ggarrange(PCAauto_12, PCAallo_12, nrow = 2, legend.grob = MyLegendBlog)


# Saves plot ~
ggsave(PCA_Plot, file = "Passersp.enomics--PCA.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 12, height = 12, dpi = 600)
ggsave(PCA_Plot, file = "Passersp.Genomics--PCA.jpeg",
      limitsize = FALSE, scale = 1, width = 12, height = 12, dpi = 600)


#
##
### The END ~~~~~