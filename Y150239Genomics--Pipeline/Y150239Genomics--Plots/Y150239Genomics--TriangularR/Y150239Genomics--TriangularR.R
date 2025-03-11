### The BEGINNING ~~~~~
##
# Plots Passer sp. Genomics -- TriangularR | Written by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads packages ~
devtools::install_github("omys-omics/triangulaR")
pacman::p_load(tidyverse, ggstar, ggforce, vcfR, triangulaR, ggh4x, ggrepel, grid, gtable)


# Loads VCF data ~
VCF_auto <- read.vcfR("../../../../LargeFiles/Y150239--TriangularR/AllSamples_bcftools.raw.vcf.Filtered.Focal.Autosomes.ALL.vcf", verbose = TRUE)
VCF_allo <- read.vcfR("../../../../LargeFiles/Y150239--TriangularR/AllSamples_bcftools.raw.vcf.Filtered.Focal.Allosome.ALL.vcf", verbose = TRUE)


# Loads annotation file ~
annot_auto <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Focal.Autosomes.ALL.annot",  sep = " ", header = FALSE, stringsAsFactors = FALSE, col.names = c("id", "pop"))
annot_allo <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Focal.Allosome.ALL.annot",  sep = " ", header = FALSE, stringsAsFactors = FALSE, col.names = c("id", "pop"))


# Gets AIMs ~
VCF_auto.diff9 <- alleleFreqDiff(vcfR = VCF_auto, pm = annot_auto, p1 = "House", p2 = "Spanish", difference = 0.9)
VCF_allo.diff9 <- alleleFreqDiff(vcfR = VCF_allo, pm = annot_allo, p1 = "House", p2 = "Spanish", difference = 0.9)


# Gets AIMs´ genotypes ~
m_auto <- extract.gt(VCF_auto.diff9)
m_allo <- extract.gt(VCF_allo.diff9)


# Recodes to allele counts ~ 
m_auto[m_auto=="0|0"] <- 0
m_auto[m_auto=="0|1"] <- 1
m_auto[m_auto=="1|0"] <- 1
m_auto[m_auto=="1|1"] <- 2
m_auto[m_auto=="0/0"] <- 0
m_auto[m_auto=="0/1"] <- 1
m_auto[m_auto=="1/0"] <- 1
m_auto[m_auto=="1/1"] <- 2
m_allo[m_allo=="0|0"] <- 0
m_allo[m_allo=="0|1"] <- 1
m_allo[m_allo=="1|0"] <- 1
m_allo[m_allo=="1|1"] <- 2
m_allo[m_allo=="0/0"] <- 0
m_allo[m_allo=="0/1"] <- 1
m_allo[m_allo=="1/0"] <- 1
m_allo[m_allo=="1/1"] <- 2


# Filters $ subsets the genotypes for the two populations
p1.gts_auto <- m_auto[, annot_auto[annot_auto$pop == "House",]$id]
p2.gts_auto <- m_auto[, annot_auto[annot_auto$pop == "Spanish",]$id]
p1.gts_allo <- m_allo[, annot_allo[annot_allo$pop == "House",]$id]
p2.gts_allo <- m_allo[, annot_allo[annot_allo$pop == "Spanish",]$id]


# Converts to numeric ~
p1.gts_auto[] <- sapply(p1.gts_auto, as.numeric)
p2.gts_auto[] <- sapply(p2.gts_auto, as.numeric)
p1.gts_allo[] <- sapply(p1.gts_allo, as.numeric)
p2.gts_allo[] <- sapply(p2.gts_allo, as.numeric)


# Calculates allele frequencies for P1 & P2 ~
af_p1_auto <- (rowSums(p1.gts_auto == 1, na.rm = TRUE) + (2 * rowSums(p1.gts_auto == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts_auto)))
af_p2_auto <- (rowSums(p2.gts_auto == 1, na.rm = TRUE) + (2 * rowSums(p2.gts_auto == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts_auto)))
af_p1_allo <- (rowSums(p1.gts_allo == 1, na.rm = TRUE) + (2 * rowSums(p1.gts_allo == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts_allo)))
af_p2_allo <- (rowSums(p2.gts_allo == 1, na.rm = TRUE) + (2 * rowSums(p2.gts_allo == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts_allo)))


# Determines P1 & P2 alleles based on allele frequencies ~
p1.allele_auto <- ifelse(af_p1_auto > af_p2_auto, 2, 0)
p2.allele_auto <- ifelse(af_p2_auto > af_p1_auto, 2, 0)
p1.allele_allo <- ifelse(af_p1_allo > af_p2_allo, 2, 0)
p2.allele_allo <- ifelse(af_p2_allo > af_p1_allo, 2, 0)


# Creates a matrix to store hybrid index scores ~
n_auto <- matrix(nrow = nrow(m_auto), ncol = ncol(m_auto))
n_allo <- matrix(nrow = nrow(m_allo), ncol = ncol(m_allo))


# Compares genotypes & assigns scores ~
n_auto[m_auto == p1.allele_auto] <- 0
n_auto[m_auto == 1] <- 1
n_auto[m_auto == p2.allele_auto] <- 2
n_auto[is.na(m_auto)] <- NA
n_auto[m_auto == -9] <- NA
n_allo[m_allo == p1.allele_allo] <- 0
n_allo[m_allo == 1] <- 1
n_allo[m_allo == p2.allele_allo] <- 2
n_allo[is.na(m_allo)] <- NA
n_allo[m_allo == -9] <- NA


# Names columns & rows ~
colnames(n_auto) <- colnames(m_auto)
rownames(n_auto) <- rownames(m_auto)
colnames(n_allo) <- colnames(m_allo)
rownames(n_allo) <- rownames(m_allo)


# Fills matrix ~
n_auto <- as.data.frame(n_auto)
n_allo <- as.data.frame(n_allo)


# Selects focal individual & controls ~
fulldf_auto <- n_auto %>% select(PI22NLD0001M_SAMPLE, PD22NLD0146F_SAMPLE, PD22NLD0147F_SAMPLE, PDOM2022NLD0077M_SAMPLE)
fulldf_allo <- n_allo %>% select(PI22NLD0001M_SAMPLE, PDOM2022NLD0077M_SAMPLE)


# Convert row names to a column ~
fulldf_auto <- fulldf_auto %>%
               mutate(CHR = sub("_.*", "", rownames(fulldf_auto))) %>%
               mutate(POS = sub(".*_", "", rownames(fulldf_auto))) %>%
               tibble::rownames_to_column(var = "SNP") %>%
               select(SNP, CHR, POS, everything())

fulldf_allo <- fulldf_allo %>%
               mutate(CHR = sub("_.*", "", rownames(fulldf_allo))) %>%
               mutate(POS = sub(".*_", "", rownames(fulldf_allo))) %>%
               tibble::rownames_to_column(var = "SNP") %>%
               select(SNP, CHR, POS, everything())


# Creates Index per CHR ~
fulldf_auto$Index <- with(fulldf_auto, ave(seq_along(CHR), CHR, FUN = seq_along))

fulldf_allo$Index <- with(fulldf_allo, ave(seq_along(CHR), CHR, FUN = seq_along))


# Converts to wide data frame ~
fulldf_auto <- gather(fulldf_auto, Individual, Ancestry,
                 "PI22NLD0001M_SAMPLE", "PD22NLD0146F_SAMPLE", "PD22NLD0147F_SAMPLE", "PDOM2022NLD0077M_SAMPLE")


fulldf_allo <- gather(fulldf_allo, Individual, Ancestry,
                      "PI22NLD0001M_SAMPLE", "PDOM2022NLD0077M_SAMPLE")


# Combines the DFs ~
fulldf <- rbind(fulldf_auto, fulldf_allo)


# Expands PCA_Annot by adding Population ~
fulldf$Ancestry <- ifelse(grepl("0", fulldf$Ancestry), "House",
                   ifelse(grepl("1", fulldf$Ancestry), "Heterozygous",
                   ifelse(grepl("2", fulldf$Ancestry), "Spanish", "Error")))


# Expands PCA_Annot by adding Population ~
fulldf$Individual <- ifelse(grepl("PI22NLD0001M_SAMPLE", fulldf$Individual), "Y150239",
                     ifelse(grepl("PD22NLD0146F_SAMPLE", fulldf$Individual), "Garderen_01",
                     ifelse(grepl("PD22NLD0147F_SAMPLE", fulldf$Individual), "Garderen_02",
                     ifelse(grepl("PDOM2022NLD0077M_SAMPLE", fulldf$Individual), "Meerkerk_01", "Error"))))


# Reorders BioStatus ~
fulldf$Individual <- factor(fulldf$Individual, ordered = TRUE,
                            levels = c("Meerkerk_01",
                                       "Garderen_02",
                                       "Garderen_01",
                                       "Y150239"))


# Reorders Ancestry ~
fulldf$Ancestry <- factor(fulldf$Ancestry, ordered = TRUE,
                           levels = c("House",
                                      "Heterozygous",
                                      "Spanish"))


# Reorders CHR ~
fulldf$CHR <- factor(fulldf$CHR, ordered = TRUE,
                     levels = c("chr1", "chr1A", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                "chr11", "chr12", "chr13", "chr14", "chr15", "chr17", "chr18", "chr19", "chr20", "chr21", 
                                "chr22", "chr23", "chr24", "chr26", "chr27", "chr28", "chrZ", "scaffold00239"))


# Fixes CHRs´ names ~
y_strip_labels <- setNames(c("CHR 01", "CHR 01A", "CHR 02", "CHR 03", "CHR 04", "CHR 05", "CHR 06", "CHR 07", 
                             "CHR 08", "CHR 09", "CHR 10", "CHR 11", "CHR 12", "CHR 13", "CHR 14", "CHR 15", 
                             "CHR 17", "CHR 18", "CHR 19", "CHR 20", "CHR 21", "CHR 22", "CHR 23", "CHR 24", 
                             "CHR 26", "CHR 27", "CHR 28", "CHR Z", "SD00169", "SD00221", "SD00223", "SD00224", 
                             "SD00238", "SD00239", "SD00242"), 
                           c("chr1", "chr1A", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                             "chr11", "chr12", "chr13", "chr14", "chr15", "chr17", "chr18", "chr19", "chr20", "chr21", 
                             "chr22", "chr23", "chr24", "chr26", "chr27", "chr28", "chrZ", "scaffold00169", 
                             "scaffold00221", "scaffold00223", "scaffold00224", "scaffold00238", "scaffold00239", 
                             "scaffold00242"))


# Marks y-axis labels for no display ~ 
fulldf <- fulldf %>% mutate(Individual = if_else(!(CHR %in% c("chr1", "chrZ")), paste0(as.character(Individual), "no_display"), 
                            as.character(Individual)))


# Marks y-axis labels for no display ~ 
fulldf$POSMb <- as.numeric(fulldf$POS) / 1000000


# Little function to suppress y-axis labels ~
delete_no_display <- function(v) {
  if_else(str_detect(v, 'no_display'), '', v)}


# Creates Index plot ~
AncestryPlot_Index <-
ggplot(fulldf, aes(x = Index, y = Individual, fill = as.factor(Ancestry))) +
  geom_point(shape = 21, size = 2, colour = "#000000", stroke = 0) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000"), na.translate = FALSE) +
  ggtitle("Ancestry-informative Markers") +
  scale_x_discrete(expand = c(.005, .005)) +
  scale_y_discrete(labels = delete_no_display) +
  facet_grid2(CHR ~ ., scales = "free", axes = "all", remove_labels = "x", labeller = labeller(CHR = y_strip_labels)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(.05, "cm"),
        legend.position = c(.8, .875),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 0, b = 0, r = 0, l = 0),
        plot.title = element_text(family = "Optima", size = 20, face = "bold", color = "#000000", hjust = .5, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "Optima", color = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 8, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Ancestry", title.theme = element_text(family = "Optima", size = 12, face = "bold"),
                             label.theme = element_text(family = "Optima", size = 10), override.aes = list(shape = 21, size = 4, stroke = .15)))


# Saves Index plot ~
ggsave(AncestryPlot_Index, file = "Y150239Genomics--AncestryHeatmap_AIMs.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 10, height = 12, dpi = 600)
ggsave(AncestryPlot_Index, file = "Y150239Genomics--AncestryHeatmap_AIMs.jpeg",
       limitsize = FALSE, scale = 1, width = 10, height = 12, dpi = 600)


# Calculates differentiation indexes ~
HI_HET_auto.diff9 <- hybridIndex(vcfR = VCF_auto.diff9, pm = annot_auto, p1 = "House", p2 = "Spanish")
HI_HET_allo.diff9 <- hybridIndex(vcfR = VCF_allo.diff9, pm = annot_allo, p1 = "House", p2 = "Spanish")


# Expands HI_HET ~
HI_HET_auto.diff9$CHRType <- "Autosomes"
HI_HET_auto.diff9$Diff <- "0.90"
HI_HET_auto.diff9$SNPs <- nrow(VCF_auto.diff9@fix)
HI_HET_allo.diff9$CHRType <- "Chromosome Z"
HI_HET_allo.diff9$Diff <- "0.90"
HI_HET_allo.diff9$SNPs <- nrow(VCF_allo.diff9@fix)


# Combines fulldf_auto and fulldf_allo ~
HI_HET <- rbind(HI_HET_auto.diff9, HI_HET_allo.diff9)


# Expands PCA_Annot by adding Population ~
HI_HET$Population <- ifelse(grepl("FR0", HI_HET$id), "Sales",
                     ifelse(grepl("KAZ", HI_HET$id), "Chokpak",
                     ifelse(grepl("Lesina", HI_HET$id), "Lesina",
                     ifelse(grepl("Crotone", HI_HET$id), "Crotone",
                     ifelse(grepl("Guglionesi", HI_HET$id), "Guglionesi",
                     ifelse(grepl("PI22NLD0001M", HI_HET$id), NA,
                     ifelse(grepl("PD22NLD0146F", HI_HET$id), NA,
                     ifelse(grepl("PD22NLD0147F", HI_HET$id), NA,
                     ifelse(grepl("PDOM2022NLD0077M", HI_HET$id), NA,
                     ifelse(grepl("PDOM2022NLD0", HI_HET$id), "Utrecht", "Error"))))))))))


# Reorders Population ~
HI_HET$Population <- factor(HI_HET$Population, ordered = T,
                            levels = c("Utrecht",
                                       "Sales",
                                       "Crotone",
                                       "Guglionesi",
                                       "Lesina",
                                       "Chokpak",
                                       NA))


# Expands fulldf by adding Species ~
HI_HET$Species <- ifelse(HI_HET$Population %in% c("Utrecht", "Sales"), "House",
                  ifelse(HI_HET$Population %in% c("Chokpak", "Lesina"), "Spanish",
                  ifelse(HI_HET$Population %in% c("Crotone", "Guglionesi"), "Italian",
                  ifelse(HI_HET$Population %in% c(NA), NA, "Error"))))


# Reorders Population ~
HI_HET$Species <- factor(HI_HET$Species, ordered = T,
                         levels = c("House",
                                    "Italian",
                                    "Spanish",
                                    NA))


# Expands PCA_Annot by adding Labels ~
HI_HET$Labels <- ifelse(HI_HET$pop %in% c("Y150239"), "Focal Ind.",
                 ifelse(HI_HET$pop %in% c("Meerkerk"), "Meerkerk_01",
                 ifelse(HI_HET$id %in% c("PD22NLD0146F_SAMPLE"), "Garderen_01",
                 ifelse(HI_HET$id %in% c("PD22NLD0147F_SAMPLE"), "Garderen_02", ""))))


# Creates triangle ~
triangle <- data.frame(x = c(0, 1, 0.5, 0), y = c(0, 0, 1, 0))


# Creates semicircle ~
x_vals <- seq(0, 1, length.out = 100)
semicircle <- data.frame(x = x_vals, y = 2 * x_vals * (1 - x_vals))


# Create the plot
Panel <-
ggplot() +
  geom_path(data = triangle, aes(x = x, y = y), color = "#000000", linetype = 2, linewidth = .35) +
  geom_path(data = semicircle, aes(x = x, y = y), color = "#000000", linetype = 2, linewidth = .35) +
  geom_star(data = HI_HET, aes(x = hybrid.index, y = heterozygosity, fill = Species), starshape = 15, colour = "#000000", size = 2, starstroke = .15, alpha = .75) +
  #geom_star(data = fulldf, aes(x = hybrid.index, y = heterozygosity, fill = perc.missing), starshape = 15, colour = "#000000", size = 2, starstroke = .15, alpha = .75) +
  facet_grid2(CHRType ~ ., scales = "free_y", axes = "all", remove_labels = "x") +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000"), na.translate = FALSE) +
  geom_label(data = HI_HET, aes(x = .15, y = .85, label = paste0("# of AIMs: ", scales::comma(SNPs))), alpha = 1,
             size = 4.5, fontface = "bold", fill = "#d6d6d6", label.padding = unit(.5, "lines"), family = "Optima", show.legend = FALSE) +
  geom_label_repel(data = subset(HI_HET, CHRType == "Autosomes" & Labels == "Focal Ind."), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = -.05, nudge_y = .2,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(HI_HET, CHRType == "Autosomes" & Labels == "Garderen_01"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = -.01, nudge_y = .16,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(HI_HET, CHRType == "Autosomes" & Labels == "Garderen_02"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .2, nudge_y = 0,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(HI_HET, CHRType == "Autosomes" & Labels == "Meerkerk_01"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .2, nudge_y = .08,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(HI_HET, CHRType == "Chromosome Z" & Labels == "Focal Ind."), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .1, nudge_y = 0,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(HI_HET, CHRType == "Chromosome Z" & Labels == "Meerkerk_01"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = 0, nudge_y = .175,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  scale_x_continuous("Hybird Index",
                     breaks = c(.25, .5, .75),
                     labels = c("0.25", "0.50", "0.75"),
                     limits = c(0, 1),
                     expand = c(.01, .01)) +
  scale_y_continuous("Interclass Heterozygozity",
                     breaks = c(.25, .5, .75),
                     labels = c("0.25", "0.50", "0.75"),
                     limits = c(0, 1),
                     expand = c(.01, .01)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#E5E7E9", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.y = unit(.2, "cm"),
        legend.position = "top",
        legend.box = "vertical",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 10, r = 0, l = 0),
        axis.title.x = element_text(family = "Optima", size = 16, face = "bold", margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Optima", size = 16, face = "bold", margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.text = element_text(family = "Optima", color = "#000000", size = 11, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Species", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                                      label.theme = element_text(size = 14, family = "Optima"),
                                      override.aes = list(starshape = 15, size = 5, starstroke = .15), nrow = 1, order = 1),
         starshape = "none",
         colour = "none")


# Saves plot ~
ggsave(Panel, file = "Y150239Genomics--Triangular.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 10, height = 12, dpi = 600)
ggsave(Panel, file = "Y150239Genomics--Triangular.jpeg",
       limitsize = FALSE, scale = 1, width = 10, height = 12, dpi = 600)


#
##
### The END ~~~~~