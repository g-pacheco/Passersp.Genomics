### The BEGINNING ~~~~~
##
# ~ Plots Passer sp. Genomics -- PopStructure Panel | Written by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, ggpattern, ggh4x, scales, optparse, plyr, RColorBrewer, extrafont, gtable, grid, ggtext, glue,
               patchwork, cowplot, grImport2, pdftools, magick, gtable)


# Loads PCA_Plot to create article panel ~
PCA_Plot <- readRDS("PCA_Plot.rds")
ngsAdmix <- readRDS("Passersp.Genomics--ngsAdmix_Abridged.rds")


# Adds taps ~
PCA_Plot <- PCA_Plot + labs(tag = "A)") + theme(plot.tag.position = c(.015, .98), plot.tag = element_text(family = "Optima", face = "bold", size = 26))
ngsAdmix <- ngsAdmix + labs(tag = "B)") + theme(plot.tag.position = c(.015, 1.075), plot.tag = element_text(family = "Optima", face = "bold", size = 26))


# Creates PopStructure panel ~
vertical_spacer <- plot_spacer()
top_row <- PCA_Plot
top_row_alt <- plot_spacer() + PCA_Plot + plot_spacer() + plot_layout(widths = c(1, 4, 1))
PopStructurePanel_Plot <- top_row / vertical_spacer / ngsAdmix + plot_layout(heights = c(3, .08, 1))


# Saves PopStructure panel ~
ggsave(PopStructurePanel_Plot, file = "Passersp.Genomics--PopStructurePanel.pdf", device = cairo_pdf,
       width = 26, height = 20, scale = 1, dpi = 600)
ggsave(PopStructurePanel_Plot, file = "Passersp.Genomics--PopStructurePanel.jpeg",
       width = 26, height = 20, scale = 1, dpi = 600)


#
##
### The END ~~~~~