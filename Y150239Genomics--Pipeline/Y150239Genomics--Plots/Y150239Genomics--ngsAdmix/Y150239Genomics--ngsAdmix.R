### The BEGINNING ~~~~~
##
# ~ Plots Y150239Genomics--ngsAdmix <> Written by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, ggpattern, ggh4x, scales, optparse, plyr, RColorBrewer, extrafont, gtable, grid, ggtext, glue)


# Creates colour palette ~
nb.cols <- 15
MyColours <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)


# Loads the data ~
samples.auto <- read.table("Y150239Genomics--ngsAdmix.Autosomes.popfile", stringsAsFactors = FALSE, sep = "\t")
samples.allo <- read.table("Y150239Genomics--ngsAdmix.Allosome.popfile", stringsAsFactors = FALSE, sep = "\t")


# Reads the annotation file ~
ids.auto <- read.table("./LABELs/AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.labels",
                       stringsAsFactors = FALSE, sep = "\t", header = FALSE, col.names = c("Sample_ID"))
ids.allo <- read.table("./LABELs/AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MalesOnly.MAFfiltered.Pruned.K2.labels",
                       stringsAsFactors = FALSE, sep = "\t", header = FALSE, col.names = c("Sample_ID"))


# Expands ids.auto by adding Population ~
ids.auto$Population <- ifelse(grepl("FR0", ids.auto$Sample_ID), "Sales",
                       ifelse(grepl("KAZ", ids.auto$Sample_ID), "Chokpak",
                       ifelse(grepl("Lesina", ids.auto$Sample_ID), "Lesina",
                       ifelse(grepl("Crotone", ids.auto$Sample_ID), "Crotone",
                       ifelse(grepl("Guglionesi", ids.auto$Sample_ID), "Guglionesi",
                       ifelse(grepl("PI22NLD0001M", ids.auto$Sample_ID), "Focal Ind.",
                       ifelse(grepl("PD22NLD0146F", ids.auto$Sample_ID), "Garderen",
                       ifelse(grepl("PD22NLD0147F", ids.auto$Sample_ID), "Garderen",
                       ifelse(grepl("PDOM2022NLD0077M", ids.auto$Sample_ID), "Meerkerk",
                       ifelse(grepl("PDOM2022NLD0", ids.auto$Sample_ID), "Utrecht", "Error"))))))))))


# Expands ids.allo by adding Population ~
ids.allo$Population <- ifelse(grepl("FR0", ids.allo$Sample_ID), "Sales",
                       ifelse(grepl("KAZ", ids.allo$Sample_ID), "Chokpak",
                       ifelse(grepl("Lesina", ids.allo$Sample_ID), "Lesina",
                       ifelse(grepl("Crotone", ids.allo$Sample_ID), "Crotone",
                       ifelse(grepl("Guglionesi", ids.allo$Sample_ID), "Guglionesi",
                       ifelse(grepl("PI22NLD0001M", ids.allo$Sample_ID), "Focal Ind.",
                       ifelse(grepl("PDOM2022NLD0077M", ids.allo$Sample_ID), "Meerkerk",
                       ifelse(grepl("PDOM2022NLD0", ids.allo$Sample_ID), "Utrecht", "Error"))))))))


# Reorders Population ~
ids.auto$Population <- factor(ids.auto$Population, ordered = T,
                              levels = c("Utrecht",
                                         "Sales",
                                         "Crotone",
                                         "Guglionesi",
                                         "Lesina",
                                         "Chokpak",
                                         "Focal Ind.",
                                         "Garderen",
                                         "Meerkerk"))


# Reorders Population ~
ids.allo$Population <- factor(ids.allo$Population, ordered = T,
                              levels = c("Utrecht",
                                         "Sales",
                                         "Crotone",
                                         "Guglionesi",
                                         "Lesina",
                                         "Chokpak",
                                         "Focal Ind.",
                                         "Garderen",
                                         "Meerkerk"))


# Expands PCA_Annot by adding Species ~
ids.auto$Species <- ifelse(ids.auto$Population %in% c("Utrecht", "Sales", "Garderen", "Meerkerk"), "House",
                    ifelse(ids.auto$Population %in% c("Chokpak", "Lesina"), "Spanish",
                    ifelse(ids.auto$Population %in% c("Guglionesi", "Crotone"), "Italian",
                    ifelse(ids.auto$Population %in% c("Y150239"), "Focal Ind.", "Error"))))


# Expands PCA_Annot by adding Species ~
ids.allo$Species <- ifelse(ids.allo$Population %in% c("Utrecht", "Sales", "Meerkerk"), "House",
                    ifelse(ids.allo$Population %in% c("Chokpak", "Lesina"), "Spanish",
                    ifelse(ids.allo$Population %in% c("Guglionesi", "Crotone"), "Italian",
                    ifelse(ids.allo$Population %in% c("Y150239"), "Focal Ind.", "Error"))))


# Reorders Population ~
ids.auto$Species <- factor(ids.auto$Species, ordered = T,
                           levels = c("House",
                                      "Italian",
                                      "Spanish",
                                      "Focal Ind."))


# Reorders Population ~
ids.allo$Species <- factor(ids.allo$Species, ordered = T,
                           levels = c("House",
                                      "Italian",
                                      "Spanish",
                                      "Focal Ind."))


# Recognises chromosomes ~
ids.auto$chrtype <- "Autosomes"
ids.allo$chrtype <- "Chromosome Z"


# Creates data frame ~
fulldf.auto <- data.frame()
fulldf.allo <- data.frame()


x.auto <- list(c(3, 2, 1, 4, 6, 5, 7),
               #c(4, 5, 1, 2, 6, 3),
               #c(4, 2, 5, 3, 1),
               #c(4, 3, 1, 2),
               c(1, 3, 2),
               c(1, 2))


x.allo <- list(c(3, 5, 1, 6, 7, 4, 2),
               #c(1, 4, 2, 6, 3, 5),
               #c(4, 3, 2, 1, 5),
               #c(3, 2, 1, 4),
               c(3, 1, 2),
               c(1, 2))


# Defines samples' IDs ~
sampleid = "Sample_ID"


# Loops over all Ks while adding labels and reordering clusters ~
for (j in 1:length(samples.auto[, 1])){
  data <- read.table(samples.auto[j, 1])[, x.auto[[j]]]
  for (i in 1:dim(data)[2]) { 
    temp <- data.frame(Ancestry = data[, i])
    temp$K <- as.factor(rep(i, times = length(temp$Ancestry)))
    temp[sampleid] <- as.factor(ids.auto[sampleid][, 1])
    temp$K_Value <- as.factor(rep(paste("K = ", dim(data)[2], sep = ""), times = length(temp$Ancestry)))
    temp <- merge(ids.auto, temp)
    fulldf.auto <- rbind(fulldf.auto, temp)}}


for (j in 1:length(samples.allo[, 1])){
  data <- read.table(samples.allo[j, 1])[, x.allo[[j]]]
  for (i in 1:dim(data)[2]) { 
    temp <- data.frame(Ancestry = data[, i])
    temp$K <- as.factor(rep(i, times = length(temp$Ancestry)))
    temp[sampleid] <- as.factor(ids.allo[sampleid][, 1])
    temp$K_Value <- as.factor(rep(paste("K = ", dim(data)[2], sep = ""), times = length(temp$Ancestry)))
    temp <- merge(ids.allo, temp)
    fulldf.allo <- rbind(fulldf.allo, temp)}}


# Finds missing rows in fulldf.allo ~ 
fulldf.allo_missing_rows <- anti_join(fulldf.auto, fulldf.allo, by = c("Sample_ID", "Population", "Species", "K", "K_Value"))


# Adds missing rows to fulldf.allo ~
fulldf.allo_missing_rows <- fulldf.allo_missing_rows %>%
                            group_by(K_Value) %>%
                            mutate(Ancestry = ifelse(K == 1, 1, 0),
                            chrtype = "Chromosome Z",
                            Status = "Missing") %>%
                            ungroup()


# Adds Status in fulldf.allo_present_rows ~ 
fulldf.allo_present_rows <- fulldf.allo %>%
                            mutate(Status = "Present")


# Combines fulldf.allo_present_rows with fulldf.allo_missing_rows ~ 
fulldf.allo <- bind_rows(fulldf.allo_missing_rows, fulldf.allo_present_rows)


# Adds Status in full.auto ~
fulldf.auto <- fulldf.auto %>%
               mutate(Status = "Present")


# Combines fulldf.auto with fulldf.allo ~
fulldf <- rbind(fulldf.auto, fulldf.allo)


# Expands fulldf by adding Group ~
fulldf$Group <- ifelse(fulldf$Population %in% c("Utrecht", "Sales"), "House Sparrow",
                ifelse(fulldf$Population %in% c("Guglionesi", "Crotone"), "Italian Sparrow",
                ifelse(fulldf$Population %in% c("Chokpak", "Lesina"), "Spanish Sparrow",
                ifelse(fulldf$Population %in% c("Focal Ind.", "Meerkerk", "Garderen"), "Focal Area", "Error"))))


# Reorders Population ~
fulldf$Group <- factor(fulldf$Group, ordered = T,
                           levels = c("House Sparrow",
                                      "Italian Sparrow",
                                      "Spanish Sparrow",
                                      "Focal Area"))


# Reorders chrtype ~
fulldf$chrtype <- factor(fulldf$chrtype, ordered = T,
                        levels = c("Autosomes",
                                   "Chromosome Z"))


# Defines the target to be plotted ~
target = "Population"


# Define colors for Present status based on K
color_palette <- c("#ee0000", "#1E90FF", "#FFD700", "#88419d", "#c994c7", "#FF6347", "#00BFFF")


# Sets bar fill colour ~
fulldfUp <- fulldf %>%
            mutate(K = as.factor(K),
            fill_color = if_else(Status == "Present", color_palette[as.numeric(K)], "#ffffff"))


# Creates the plot ~
ngsAdmix <-
  ggplot(fulldfUp, aes(x = Sample_ID, y = Ancestry, fill = fill_color, pattern = Status), colour = "#000000") +
  geom_col_pattern(width = .85, alpha = .7, pattern_size = .1, pattern_density = .01, pattern_spacing = .075, pattern_units = "in", pattern_colour = "#000000", pattern_fill = "#000000") +
  facet_nested(chrtype + K_Value ~ Group + Population, scales = "free_x", space = "free",
               strip = strip_nested(text_x = elem_list_text(size = c(16, 13), family = c("Optima", "Optima"), face = c("bold", "bold"), angle = c(0, 90), margins = c(1, 2, 3, 4)),
                                    background_x = elem_list_rect(fill = c("#d6d6d6", "#FAFAFA"), colour = c("#000000", "#000000"), linewidth = c(.3, .3)),
                                    by_layer_x = TRUE,
                                    text_y = elem_list_text(size = c(16, 13), family = c("Optima", "Optima"), face = c("bold", "bold")),
                                    background_y = elem_list_rect(fill = c("#d6d6d6", "#FAFAFA"), colour = c("#000000", "#000000"), linewidth = c(.3, .3)),
                                    by_layer_y = TRUE)) +
  scale_fill_identity() +
  scale_pattern_manual(values = c("stripe", "none")) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(.1, "cm"),
        plot.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x.bottom = element_blank(),
        #axis.text.x.bottom = element_text(colour = "#000000", face = "bold", angle = 90, vjust = .5, hjust = .5),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
        
#strip.text.x = element_text(family = "Optima", colour = "#000000", face = "bold", size = 14, angle = 90, margin = margin(0.5, 0.1, 0.5, 0.1, "cm")),
#trip.text.y = element_text(family = "Optima", colour = "#000000", face = "bold", size = 14, angle = 90, margin = margin(0, 0.1, 0, 0.1, "cm")))


# Adds grob ~
#ngsAdmix_G <- ggplotGrob(ngsAdmix)
#ngsAdmix_G <- gtable_add_rows(ngsAdmix_G, unit(1.25, "cm"), pos = 5)


# Adds top strips ~
#ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#1E90FF", alpha = .7, size = .75, lwd = .25)),
#               textGrob("House Sparrow", gp = gpar(cex = 1.5, fontface = 'bold', fontfamily = "Optima", col = "black"))),
#               t = 6, l = 4, b = 6, r = 12, name = c("a", "b"))
#ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#FFD700", alpha = .7, size = .5, lwd = .25)),
#               textGrob("Italian Sparrow", gp = gpar(cex = 1.5, fontface = 'bold', fontfamily = "Optima", col = "black"))),
#               t = 6, l = 14, b = 6, r = 20, name = c("a", "b"))
#ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#ee0000", alpha = .7, size = .75, lwd = .25)),
#               textGrob("Spanish Sparrow", gp = gpar(cex = 1.5, fontface = 'bold', fontfamily = "Optima", col = "black"))),
#               t = 6, l = 22, b = 6, r = 28, name = c("a", "b"))
#ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#c994c7", alpha = .7, size = .75, lwd = .25)),
#              textGrob("Focal Area", gp = gpar(cex = 1.5, fontface = 'bold', fontfamily = "Optima", col = "black"))),
#              t = 6, l = 30, b = 6, r = 39, name = c("a", "b"))


# Controls separation ~
#ngsAdmix_G <- gtable_add_rows(ngsAdmix_G, unit(2 / 10, "line"), 6)


# Creates the final plot ~
#grid.newpage()
#grid.draw(ngsAdmix_G)


# Saves the final plot ~
ggsave(ngsAdmix, file = "Y150239Genomics--ngsAdmix.pdf",
       device = cairo_pdf, width = 20, height = 8, scale = 1, dpi = 600)
ggsave(ngsAdmix, file = "Y150239Genomics--ngsAdmix.jpeg",
       width = 20, height = 8, scale = 1, dpi = 600)


#
##
### The END ~~~~~