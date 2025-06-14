### The BEGINNING ~~~~~
##
# ~ Plots Passer sp. Genomics -- Mitophynogeny | Written by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(ggtree, tidyverse, ggrepel, extrafont, treeio, phytools, ape, ggtreeExtra, ggnewscale, ggstar, reshape2)


# Imports extra fonts ~
loadfonts(device = "win", quiet = TRUE)


# Loads mitophylogeny ~
Mitophylo <- read.tree(file = "Y150239Genomics_mtGenome.Phylogeny.ML.Renamed.WithBSs.treefile")


# Reads annotations ~
Mitophylo_Annot <- read.table("Y150239Genomics_mtGenome.Phylogeny.ML.Renamed.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Mitophylo_Annot) <- c("Sample_ID")


# Expands Mitophylo_Annot by adding Population ~
Mitophylo_Annot$Population <- ifelse(grepl("Sales", Mitophylo_Annot$Sample_ID), "Sales",
                              ifelse(grepl("Chokpak", Mitophylo_Annot$Sample_ID), "Chokpak",
                              ifelse(grepl("Lesina", Mitophylo_Annot$Sample_ID), "Lesina",
                              ifelse(grepl("Crotone", Mitophylo_Annot$Sample_ID), "Crotone",
                              ifelse(grepl("Guglionesi", Mitophylo_Annot$Sample_ID), "Guglionesi",
                              ifelse(grepl("Focal Ind.", Mitophylo_Annot$Sample_ID), "Focal Ind.",
                              ifelse(grepl("Garderen", Mitophylo_Annot$Sample_ID), "Garderen",
                              ifelse(grepl("Meerkerk", Mitophylo_Annot$Sample_ID), "Meerkerk",
                              ifelse(grepl("Utrecht", Mitophylo_Annot$Sample_ID), "Utrecht",
                              ifelse(grepl("Tree Sparrow", Mitophylo_Annot$Sample_ID), "Tree Sparrow","Error"))))))))))


# Expands PCA_Annot by adding Species ~
Mitophylo_Annot$Species <- ifelse(Mitophylo_Annot$Population %in% c("Utrecht", "Sales"), "House",
                           ifelse(Mitophylo_Annot$Population %in% c("Chokpak", "Lesina"), "Spanish",
                           ifelse(Mitophylo_Annot$Population %in% c("Guglionesi", "Crotone"), "Italian",
                           ifelse(Mitophylo_Annot$Population %in% c("Focal Ind.", "Garderen", "Meerkerk"), "Focal Group",
                           ifelse(Mitophylo_Annot$Population %in% c("Tree Sparrow"), "Tree Sparrow", "Error")))))


# Reorders Population ~
Mitophylo_Annot$Species <- factor(Mitophylo_Annot$Species, ordered = TRUE,
                                  levels = c("Tree Sparrow",
                                             "House",
                                             "Italian",
                                             "Spanish",
                                             "Focal Group"))


# Melts Mitophylo_Annot ~
Mitophylo_Annot <- melt(Mitophylo_Annot)

#ggtree(Mitophylo) +
#  geom_text2(aes(label = node), hjust = -0.3, size = 3) +
#  theme_tree2()


# Roots the phylogenies ~
Mitophylo_Rooted <- root(Mitophylo, node = 132)


# Selects clades to highlight ~
Mitophylo_Outgroup <- list(Mitophylo_Outgroup = c("TreeSparrow_01", "TreeSparrow_02", "TreeSparrow_03", "TreeSparrow_04", "TreeSparrow_05",
                                                  "TreeSparrow_06", "TreeSparrow_07", "TreeSparrow_08", "TreeSparrow_09", "TreeSparrow_10"))
Mitophylo_Rooted <- groupOTU(Mitophylo_Rooted, Mitophylo_Outgroup)


# Creates Phylo_1 plot ~
Mitophylo_Plot <-
 ggtree(Mitophylo_Rooted, aes(colour = group, show.legend = TRUE), layout = "fan", size = .125, branch.length = "none") +
  geom_tiplab(align = TRUE, linesize = .02, offset = .5, size = .7, show.legend = FALSE) +
  #scale_colour_manual(name = "Species", labels = c("Columba livia", "Columba rupestris", "Columba palumbus", "Streptopelia risoria"),
  #                    values = c("#000000", "#fb8072", "#984ea3", "#33a02c")) +
  theme(panel.spacing = margin(t = 0, b = 0, r = 0, l = 0),
        plot.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 5, b = -20, r = 0, l = 0)) +
  guides(colour = guide_legend(title.theme = element_text(size = 10, face = "bold"),
                                 label.theme = element_text(size = 7.5, face = "italic"),
                                 override.aes = list(size = .35, alpha = .9)))


# Saves Phylo_1 plot ~
ggsave(Mitophylo_Plot, file = "Passersp.Genomics--Mitophynogeny.pdf",
       device = cairo_pdf, width = 4, height = 4, scale = 1.5, dpi = 600)


# Creates Phylo_2 plot ~
Mitophylo_BasePhylo <- ggtree(Mitophylo_Rooted, layout = "circular", size = .125, colour = "#000000", branch.length = "none")


# Merges annotation to base phylogeny ~
Mitophylo_BasePhylo_Annot <- Mitophylo_BasePhylo %<+% Mitophylo_Annot


# Creates Phylo_2 plot ~
Mitophylo_2_Plot <- 
  Mitophylo_BasePhylo_Annot +
  geom_point2(aes(label = label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 85), shape = 21, size = 4, fill = "#155211", colour = "#155211", alpha = .85, stroke = .075, show.legend = FALSE) +
  geom_tiplab(aes(label = gsub("TreeSparrow", "Tree Sparrow",
                               gsub("FocalInd\\.", "Focal Ind.", 
                                    gsub("_", " ", label)))), family = "Optima", colour = "#000000", align = FALSE, linesize = .02, offset = .5, size = 6, show.legend = FALSE) +
  geom_tippoint(aes(fill = Species, subset = !is.na(Species)), size = 4.5, stroke = .075, colour = "#000000", alpha = 1, shape = 21, na.rm = TRUE) +
  scale_fill_manual(values = c("#d9d9d9", "#1E90FF", "#FFD700", "#ee0000", "#c994c7")) +
  theme(panel.spacing = margin(t = 0, b = 0, r = 0, l = 0),
        panel.border = element_blank(),
        plot.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.position = c(.07, .935),
        legend.spacing.y = unit(.4, "cm"),
        legend.key.height = unit(.45, "cm"),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 5, b = -20, r = 0, l = 30)) +
  guides(fill = guide_legend(title = "Species", title.theme = element_text(family = "Optima", size = 25, face = "bold"),
                             label.theme = element_text(family = "Optima", size = 22),
                             override.aes = list(starshape = 21, size = 5, alpha = 1, starstroke = .0015), order = 2))


# Saves Mitophylo_2_Plot plot ~
ggsave(Mitophylo_2_Plot, file = "Passersp.Genomics--Mitophynogeny.pdf",
       device = cairo_pdf, width = 17, height = 17, scale = 1, dpi = 600)
ggsave(Mitophylo_2_Plot, file = "Passersp.Genomics--Mitophynogeny.jpeg",
       device = "jpeg", width = 17, height = 17, scale = 1, dpi = 600)


#
##
### The END ~~~~~