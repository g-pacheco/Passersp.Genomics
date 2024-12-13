### The BEGINNING ~~~~~
##
# Y150239Genomics--TriangularR by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads packages ~
pacman::p_load(tidyverse, ggstar, ggforce, vcfR, triangulaR, lemon, ggrepel, grid, gtable)
devtools::install_github("omys-omics/triangulaR")


# Loads VCF data ~
VCF_auto <- read.vcfR("AllSamples_bcftools.raw.vcf.Filtered.Y150239.Autosomes.ALL.vcf", verbose = TRUE)
VCF_allo <- read.vcfR("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.TriangularR.vcf", verbose = TRUE)


# Loads annotation file ~
pm <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Y150239.Autosomes.ALL.annot",  sep = " ", header = FALSE, stringsAsFactors = FALSE)
colnames(pm) <- c("id", "pop")
annot_allo <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.TriangularR.annot",  sep = " ", header = FALSE, stringsAsFactors = FALSE)
colnames(annot_allo) <- c("id", "pop")

# Process data ~
VCF_auto.diff9 <- alleleFreqDiff(vcfR = VCF_auto, pm = pm, p1 = "House", p2 = "Spanish", difference = 0.9)


m <- extract.gt(VCF_auto.diff9)


# recode to allele counts
m[m=="0|0"] <- 0
m[m=="0|1"] <- 1
m[m=="1|0"] <- 1
m[m=="1|1"] <- 2
m[m=="0/0"] <- 0
m[m=="0/1"] <- 1
m[m=="1/0"] <- 1
m[m=="1/1"] <- 2


# Filter and subset the genotypes for the two populations
p1.gts <- m[, pm[pm$pop == "House",]$id]
p2.gts <- m[, pm[pm$pop == "Spanish",]$id]


# convert to numeric
p1.gts[] <- sapply(p1.gts, as.numeric)
p2.gts[] <- sapply(p2.gts, as.numeric)


# Calculate allele frequencies for p1 and p2
af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))


# Determine p1 and p2 allele based on allele frequencies
p1.allele <- ifelse(af_p1 > af_p2, 2, 0)
p2.allele <- ifelse(af_p2 > af_p1, 2, 0)


# Create a matrix to store hybrid index scores
n <- matrix(nrow = nrow(m), ncol = ncol(m))


# Compare genotypes and assign scores
n[m == p1.allele] <- 0
n[m == 1] <- 1
n[m == p2.allele] <- 2
n[is.na(m)] <- NA
n[m == -9] <- NA


colnames(n) <- colnames(m)
rownames(n) <- rownames(m)


n <- as.data.frame(n)


Layka <- n %>% select(PI22NLD0001M_SAMPLE, PD22NLD0146F_SAMPLE, PD22NLD0147F_SAMPLE, PDOM2022NLD0077M_SAMPLE)


Layka <- Layka %>% 
  filter(!is.na(PDOM2022NLD0077M_SAMPLE) & !is.na(PD22NLD0146F_SAMPLE) & !is.na(PD22NLD0147F_SAMPLE))


# Convert row names to a column
Layka <- Layka %>%
         mutate(CHR = sub("_.*", "", rownames(Layka))) %>%
         mutate(POS = sub(".*_", "", rownames(Layka))) %>%
         tibble::rownames_to_column(var = "SNP") %>%
         select(SNP, CHR, POS, everything())


# Create the Index column, restarting for each new chromosome
Layka <- Layka %>%
         group_by(CHR) %>%
         mutate(Index = row_number()) %>%
         ungroup()


fulldf <- gather(Layka, Individual, Ancestry,
                 "PI22NLD0001M_SAMPLE", "PD22NLD0146F_SAMPLE", "PD22NLD0147F_SAMPLE", "PDOM2022NLD0077M_SAMPLE")


# Expands PCA_Annot by adding Population ~
fulldf$Ancestry <- ifelse(grepl("0", fulldf$Ancestry), "Spanish",
                   ifelse(grepl("1", fulldf$Ancestry), "Heterozygous",
                   ifelse(grepl("2", fulldf$Ancestry), "House", "Error")))


# Expands PCA_Annot by adding Population ~
fulldf$Individual <- ifelse(grepl("PI22NLD0001M_SAMPLE", fulldf$Individual), "Y150239",
                     ifelse(grepl("PD22NLD0146F_SAMPLE", fulldf$Individual), "Garderen_01",
                     ifelse(grepl("PD22NLD0147F_SAMPLE", fulldf$Individual), "Garderen_02",
                     ifelse(grepl("PDOM2022NLD0077M_SAMPLE", fulldf$Individual), "Meerkerk_01", "Error"))))


# Reorders BioStatus ~
fulldf$Individual <- factor(fulldf$Individual, ordered = T,
                          levels = c("Meerkerk_01",
                                     "Garderen_02",
                                     "Garderen_01",
                                     "Y150239"))


# Reorders BioStatus ~
fulldf$Ancestry <- factor(fulldf$Ancestry, ordered = T,
                           levels = c("House",
                                      "Heterozygous",
                                      "Spanish"))


# Step 1: Create a custom label column and remove "no_display"
fulldf$Individual_Label <- ifelse(fulldf$CHR == "chr1", as.character(fulldf$Individual), NA)


fulldf <- readRDS("Layka.rds")


# Step 2: Filter out NA values from Individual_Label (those we don't want displayed)
valid_labels <- fulldf %>%
  filter(!is.na(Individual_Label)) %>%
  select(Individual_Label, Individual) %>%
  distinct()

# Step 3: Reorder CHR as a factor
fulldf$CHR <- factor(fulldf$CHR, ordered = TRUE,
                     levels = c("chr1", "chr1A", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                "chr11", "chr12", "chr13", "chr14", "chr15", "chr17", "chr18", "chr19", "chr20", "chr21", 
                                "chr22", "chr23", "chr24", "chr26", "chr27", "chr28", "scaffold00239"))


# Step 4: Create a label mapping for facet grid
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
fulldf <- fulldf %>%
          mutate(Individual = if_else(CHR != "chr1", paste0(as.character(Individual), "no_display"), as.character(Individual)))


# Little function to suppress y-axis labels ~
delete_no_display <- function(v) {
  if_else(str_detect(v, 'no_display'), '', v)}


saveRDS(fulldf, "Layka.rds")

# Creates Index plot ~
AncestryPlot_Index <-
ggplot(fulldf, aes(x = Index, y = Individual, fill = as.factor(Ancestry))) +
  geom_point(shape = 21, size = 2, colour = "#000000", stroke = 0) +
  scale_fill_manual(values = c("#ee0000", "#FFD700", "#1E90FF"), na.translate = FALSE) +
  ggtitle("Ancestry-informative Markers") +
  scale_x_discrete(expand = c(.005, .005)) +
  scale_y_discrete(labels = delete_no_display) +
  facet_grid(CHR ~ ., scales = "free", labeller = labeller(CHR = y_strip_labels)) +
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
        axis.text.x = element_text(family = "Optima", color = "#000000", size = 5, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(family = "Optima", color = "#000000", size = 6.5, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .2),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 7, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .2),
        axis.line = element_line(colour = "#000000", linewidth = .2)) +
  guides(fill = guide_legend(title = "Ancestry", title.theme = element_text(family = "Optima", size = 12, face = "bold"),
                             label.theme = element_text(family = "Optima", size = 10), override.aes = list(shape = 21, size = 4, stroke = .15)))


# Saves Index plot ~
ggsave(AncestryPlot_Index, file = "Y150239Genomics--AncestryHeatmap_AIMs.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 10, height = 12, dpi = 600)
ggsave(AncestryPlot_Index, file = "Y150239Genomics--AncestryHeatmap_AIMs.png",
       limitsize = FALSE, scale = 1, width = 10, height = 12, dpi = 600)


# Creates POS plot ~
AncestryPlot_POS <-
  ggplot(fulldf, aes(x = POS, y = Individual, fill = as.factor(Ancestry))) +
  geom_point(shape = 21, size = 1.25, colour = "#000000", stroke = 0) +
  scale_fill_manual(values = c("#ee0000", "#FFD700", "#1E90FF"), na.translate = FALSE) +
  ggtitle("Ancestry-informative Markers") +
  scale_x_continuous(breaks = c(25000000, 50000000, 75000000, 100000000, 125000000),
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     expand = c(.005, .005)) +
  scale_y_discrete(labels = delete_no_display) +
  facet_grid(CHR ~ ., scales = "free", labeller = labeller(CHR = y_strip_labels)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(.05, "cm"),
        legend.position = c(.3, .25),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 0, b = 0, r = 0, l = 0),
        plot.title = element_text(family = "Optima", size = 20, face = "bold", color = "#000000", hjust = .5, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title = element_blank(),
        axis.text = element_text(family = "Optima", color = "#000000", size = 6.5, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .2),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 7, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .2),
        axis.line = element_line(colour = "#000000", linewidth = .2)) +
  guides(fill = guide_legend(title = "Ancestry", title.theme = element_text(family = "Optima", size = 12, face = "bold"),
                             label.theme = element_text(family = "Optima", size = 10), override.aes = list(shape = 21, size = 4, stroke = .15)))


# Saves Index plot ~
ggsave(AncestryPlot_POS, file = "Y150239Genomics--AncestryHeatmap_GenomicAIMs.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 20, height = 12, dpi = 600)
ggsave(AncestryPlot, file = "Y150239Genomics--AncestryHeatmapGenomicAIMs.png",
       limitsize = FALSE, scale = 1, width = 40, height = 12, dpi = 600)



# Calculates differentiation indexes ~
HI_HET_auto.diff9 <- hybridIndex(vcfR = VCF_auto.diff9, pm = pm, p1 = "House", p2 = "Spanish")


# Converts VCF into matrix ~
VCF_auto.diff9_matrix <- extract.gt(VCF_auto, element = "GT", as.numeric = FALSE)


# Creates a .GENO file based on the matrix ~
VCF_auto.diff9_geno <- as.data.frame(apply(VCF_auto.diff9_matrix, c(1, 2), function(gt) {
                             if (is.na(gt)) {
                               return(-1)}
                             else if (gt == "0/0") {
                               return(0)}
                             else if (gt == "0/1" || gt == "1/0") {
                               return(1)}
                             else if (gt == "1/1") {
                               return(2)}
                             else {return(-1)}}))


# Expands VCF_auto.diff1_geno by creating a the Population row ~
Population_Row <- ifelse(grepl("FR0", colnames(VCF_auto.diff9_geno)), "Sales",
                  ifelse(grepl("KAZ", colnames(VCF_auto.diff9_geno)), "Chokpak",
                  ifelse(grepl("Lesina", colnames(VCF_auto.diff9_geno)), "Lesina",
                  ifelse(grepl("Crotone", colnames(VCF_auto.diff9_geno)), "Crotone",
                  ifelse(grepl("Guglionesi", colnames(VCF_auto.diff9_geno)), "Guglionesi",
                  ifelse(grepl("PI22NLD0001M", colnames(VCF_auto.diff9_geno)), "Y150239",
                  ifelse(grepl("PD22NLD0146F", colnames(VCF_auto.diff9_geno)), "Garderen",
                  ifelse(grepl("PD22NLD0147F", colnames(VCF_auto.diff9_geno)), "Garderen",
                  ifelse(grepl("PDOM2022NLD0077M", colnames(VCF_auto.diff9_geno)), "Meerkerk",
                  ifelse(grepl("PDOM2022NLD0", colnames(VCF_auto.diff9_geno)), "Utrecht", "Error"))))))))))
VCF_auto.diff9_geno <- rbind(VCF_auto.diff9_geno, Population_Row)
rownames(VCF_auto.diff9_geno)[nrow(VCF_auto.diff9_geno)] <- "Population"


# Expands VCF_auto.diff1_geno by creating a the Species row ~
Species_Row <- ifelse(VCF_auto.diff9_geno["Population", ] %in% c("Utrecht", "Sales", "Garderen"), "House",
               ifelse(VCF_auto.diff9_geno["Population", ] %in% c("Chokpak", "Lesina"), "Spanish",
               ifelse(VCF_auto.diff9_geno["Population", ] %in% c("Crotone", "Guglionesi"), "Italian",
               ifelse(VCF_auto.diff9_geno["Population", ] %in% c("Meerkerk"), "Meerkerk",
               ifelse(VCF_auto.diff9_geno["Population", ] %in% c("Y150239"), "Y150239", "Error")))))
VCF_auto.diff9_geno <- rbind(VCF_auto.diff9_geno, Species_Row)
rownames(VCF_auto.diff9_geno)[nrow(VCF_auto.diff9_geno)] <- "Species"


# Define the desired order
desired_order <- c("House", "Italian", "Spanish", "Meerkerk", "Y150239")


# Extract the Population row and ensure it's a character vector
species_row <- as.character(VCF_auto.diff9_geno["Species", ])


# Convert the row to a factor with the specified order
species_row <- factor(species_row, levels = desired_order, ordered = TRUE)


# Order the columns based on the custom order
ordered_columns <- order(species_row, na.last = TRUE)


# Rearrange the columns of the data frame
VCF_auto.diff9_geno <- VCF_auto.diff9_geno[, ordered_columns]


# Expands VCF_auto.diff1_geno by creating the Check columns ~
VCF_auto.diff9_geno$CheckY150239 <- rep(NA, nrow(VCF_auto.diff9_geno))
VCF_auto.diff9_geno$CheckMeerkerk <- rep(NA, nrow(VCF_auto.diff9_geno))
VCF_auto.diff9_geno$HetY150239 <- rep(NA, nrow(VCF_auto.diff9_geno))
VCF_auto.diff9_geno$HetMeerkerk <- rep(NA, nrow(VCF_auto.diff9_geno))


# Fills the CheckY150239 column ~
chr_rows <- grep("chr", rownames(VCF_auto.diff9_geno))
for (ROW in chr_rows) {
  background <- which(grepl("House", VCF_auto.diff9_geno["Species", ]))
  Y150239 <- which(grepl("Y150239", VCF_auto.diff9_geno["Species", ]))
  if (length(background) > 0 & length(Y150239) > 0) {
    VCF_auto.diff9_geno[ROW, "CheckY150239"] <-
      ifelse(all(VCF_auto.diff9_geno[ROW, background] %in% c(0, 1, -1)) &
             all(VCF_auto.diff9_geno[ROW, Y150239] %in% c(2)), "PassAlt",
      ifelse(all(VCF_auto.diff9_geno[ROW, background] %in% c(2, 1, -1)) &
             all(VCF_auto.diff9_geno[ROW, Y150239] %in% c(0)), "PassRef", "Other"))}}


# Fills the HetMeerkerk column ~
chr_rows <- grep("chr", rownames(VCF_auto.diff9_geno))
for(ROW in chr_rows) {
  background <- which(grepl("House", VCF_auto.diff9_geno["Species", ]))
  Y150239 <- which(grepl("Y150239", VCF_auto.diff9_geno["Species", ]))
  if (length(background) > 0 & length(Y150239) > 0) {
    VCF_auto.diff9_geno[ROW, "HetY150239"] <-
      ifelse(all(VCF_auto.diff9_geno[ROW, Y150239] %in% c(1)), "PassHet", "Other")}}


table(VCF_auto.diff9_geno$HetY150239)


# Fills the CheckMeerkerk column ~
chr_rows <- grep("chr", rownames(VCF_auto.diff9_geno))
for(ROW in chr_rows) {
  background <- which(grepl("House", VCF_auto.diff9_geno["Species", ]))
  Meerkerk <- which(grepl("Meerkerk", VCF_auto.diff9_geno["Species", ]))
  if (length(background) > 0 & length(Meerkerk) > 0) {
    VCF_auto.diff9_geno[ROW, "CheckMeerkerk"] <-
      ifelse(all(VCF_auto.diff9_geno[ROW, background] %in% c(0, -1)) &
             all(VCF_auto.diff9_geno[ROW, Meerkerk] %in% c(1)), "PassAlt",
      ifelse(all(VCF_auto.diff9_geno[ROW, background] %in% c(0, -1)) &
             all(VCF_auto.diff9_geno[ROW, Meerkerk] %in% c(1)), "PassRef", "Other"))}}


# Fills the HetMeerkerk column ~
chr_rows <- grep("chr", rownames(VCF_auto.diff9_geno))
for(ROW in chr_rows) {
  background <- which(grepl("House", VCF_auto.diff9_geno["Species", ]))
  Meerkerk <- which(grepl("Meerkerk", VCF_auto.diff9_geno["Species", ]))
  if (length(background) > 0 & length(Meerkerk) > 0) {
    VCF_auto.diff9_geno[ROW, "HetMeerkerk"] <-
      ifelse(all(VCF_auto.diff9_geno[ROW, Meerkerk] %in% c(1)), "PassHet", "Other")}}


table(VCF_auto.diff9_geno$HetMeerkerk)

VCF_auto.diff95 <- alleleFreqDiff(vcfR = VCF_auto, pm = annot_auto, p1 = "House", p2 = "Spanish", difference = 0.95)
HI_HET_auto.diff95 <- hybridIndex(vcfR = VCF_auto.diff95, pm = annot_auto, p1 = "House", p2 = "Spanish")

VCF_auto.diff9 <- alleleFreqDiff(vcfR = VCF_auto, pm = annot_auto, p1 = "House", p2 = "Spanish", difference = 0.9)
HI_HET_auto.diff9 <- hybridIndex(vcfR = VCF_auto.diff9, pm = annot_auto, p1 = "House", p2 = "Spanish")

VCF_auto.diff8 <- alleleFreqDiff(vcfR = VCF_auto, pm = annot_auto, p1 = "House", p2 = "Spanish", difference = 0.8)
HI_HET_auto.diff8 <- hybridIndex(vcfR = VCF_auto.diff8, pm = annot_auto, p1 = "House", p2 = "Spanish")

VCF_auto.diff7 <- alleleFreqDiff(vcfR = VCF_auto, pm = annot_auto, p1 = "House", p2 = "Spanish", difference = 0.7)
HI_HET_auto.diff7 <- hybridIndex(vcfR = VCF_auto.diff7, pm = annot_auto, p1 = "House", p2 = "Spanish")

VCF_allo.diff1 <- alleleFreqDiff(vcfR = VCF_allo, pm = annot_allo, p1 = "House", p2 = "Spanish", difference = 1)
HI_HET_allo.diff1 <- hybridIndex(vcfR = VCF_allo.diff1, pm = annot_allo, p1 = "House", p2 = "Spanish")

VCF_allo.diff95 <- alleleFreqDiff(vcfR = VCF_allo, pm = annot_allo, p1 = "House", p2 = "Spanish", difference = 0.95)
HI_HET_allo.diff95 <- hybridIndex(vcfR = VCF_allo.diff95, pm = annot_allo, p1 = "House", p2 = "Spanish")

VCF_allo.diff9 <- alleleFreqDiff(vcfR = VCF_allo, pm = annot_allo, p1 = "House", p2 = "Spanish", difference = 0.9)
HI_HET_allo.diff9 <- hybridIndex(vcfR = VCF_allo.diff9, pm = annot_allo, p1 = "House", p2 = "Spanish")

VCF_allo.diff8 <- alleleFreqDiff(vcfR = VCF_allo, pm = annot_allo, p1 = "House", p2 = "Spanish", difference = 0.8)
HI_HET_allo.diff8 <- hybridIndex(vcfR = VCF_allo.diff8, pm = annot_allo, p1 = "House", p2 = "Spanish")

VCF_allo.diff7 <- alleleFreqDiff(vcfR = VCF_allo, pm = annot_allo, p1 = "House", p2 = "Spanish", difference = 0.7)
HI_HET_allo.diff7 <- hybridIndex(vcfR = VCF_allo.diff7, pm = annot_allo, p1 = "House", p2 = "Spanish")


# Loads data ~
HI_HET_auto.diff1$CHRType <- "Autosomes"
HI_HET_auto.diff1$Diff <- "1.00"
HI_HET_auto.diff1$SNPs <- nrow(VCF_auto.diff1@fix)

HI_HET_auto.diff95$CHRType <- "Autosomes"
HI_HET_auto.diff95$Diff <- "0.95"
HI_HET_auto.diff95$SNPs <- nrow(VCF_auto.diff95@fix)

HI_HET_auto.diff9$CHRType <- "Autosomes"
HI_HET_auto.diff9$Diff <- "0.90"
HI_HET_auto.diff9$SNPs <- nrow(VCF_auto.diff9@fix)

HI_HET_auto.diff8$CHRType <- "Autosomes"
HI_HET_auto.diff8$Diff <- "0.80"
HI_HET_auto.diff8$SNPs <- nrow(VCF_auto.diff8@fix)

HI_HET_auto.diff7$CHRType <- "Autosomes"
HI_HET_auto.diff7$Diff <- "0.70"
HI_HET_auto.diff7$SNPs <- nrow(VCF_auto.diff7@fix)

HI_HET_allo.diff1$CHRType <- "Chromosome Z"
HI_HET_allo.diff1$Diff <- "1.00"
HI_HET_allo.diff1$SNPs <- nrow(VCF_allo.diff1@fix)

HI_HET_allo.diff95$CHRType <- "Chromosome Z"
HI_HET_allo.diff95$Diff <- "0.95"
HI_HET_allo.diff95$SNPs <- nrow(VCF_allo.diff95@fix)

HI_HET_allo.diff9$CHRType <- "Chromosome Z"
HI_HET_allo.diff9$Diff <- "0.90"
HI_HET_allo.diff9$SNPs <- nrow(VCF_allo.diff9@fix)

HI_HET_allo.diff8$CHRType <- "Chromosome Z"
HI_HET_allo.diff8$Diff <- "0.80"
HI_HET_allo.diff8$SNPs <- nrow(VCF_allo.diff8@fix)

HI_HET_allo.diff7$CHRType <- "Chromosome Z"
HI_HET_allo.diff7$Diff <- "0.70"
HI_HET_allo.diff7$SNPs <- nrow(VCF_allo.diff7@fix)


# Combines fulldf_auto and fulldf_allo ~
fulldf <- rbind(HI_HET_auto.diff9)
#HI_HET_auto.diff95, HI_HET_auto.diff9, HI_HET_auto.diff8, HI_HET_auto.diff7, HI_HET_allo.diff1, HI_HET_allo.diff95, HI_HET_allo.diff9, HI_HET_allo.diff8, HI_HET_allo.diff7


#cols <- c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b783")
#triangle.plot(HI_HET, colors = cols)


# Expands PCA_Annot by adding Population ~
fulldf$Population <- ifelse(grepl("FR0", fulldf$id), "Sales",
                     ifelse(grepl("KAZ", fulldf$id), "Chokpak",
                     ifelse(grepl("Lesina", fulldf$id), "Lesina",
                     ifelse(grepl("Crotone", fulldf$id), "Crotone",
                     ifelse(grepl("Guglionesi", fulldf$id), "Guglionesi",
                     ifelse(grepl("PI22NLD0001M", fulldf$id), NA,
                     ifelse(grepl("PD22NLD0146F", fulldf$id), "Garderen",
                     ifelse(grepl("PD22NLD0147F", fulldf$id), "Garderen",
                     ifelse(grepl("PDOM2022NLD0077M", fulldf$id), "Meerkerk",
                     ifelse(grepl("PDOM2022NLD0", fulldf$id), "Utrecht", "Error"))))))))))


# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                            levels = c("Utrecht",
                                       "Garderen",
                                       "Meerkerk",
                                       "Sales",
                                       "Crotone",
                                       "Guglionesi",
                                       "Lesina",
                                       "Chokpak",
                                       NA))


# Expands fulldf by adding Species ~
fulldf$Species <- ifelse(fulldf$Population %in% c("Utrecht", "Sales"), "House",
                  ifelse(fulldf$Population %in% c("Chokpak", "Lesina"), "Spanish",
                  ifelse(fulldf$Population %in% c("Crotone", "Guglionesi"), "Italian",
                  ifelse(fulldf$Population %in% c(NA, "Garderen", "Meerkerk"), NA, "Error"))))


# Reorders Population ~
fulldf$Species <- factor(fulldf$Species, ordered = T,
                         levels = c("House",
                                    "Italian",
                                    "Spanish",
                                    NA))


# Expands PCA_Annot by adding Labels ~
fulldf$Labels <- ifelse(fulldf$pop %in% c("Y150239"), "Y150239",
                 ifelse(fulldf$pop %in% c("Meerkerk"), "Meerkerk_01",
                 ifelse(fulldf$id %in% c("PD22NLD0146F_SAMPLE"), "Garderen_01",
                 ifelse(fulldf$id %in% c("PD22NLD0147F_SAMPLE"), "Garderen_02", ""))))


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
  geom_star(data = fulldf, aes(x = hybrid.index, y = heterozygosity, fill = Species), starshape = 15, colour = "#000000", size = 2, starstroke = .15, alpha = .75) +
  #geom_star(data = fulldf, aes(x = hybrid.index, y = heterozygosity, fill = perc.missing), starshape = 15, colour = "#000000", size = 2, starstroke = .15, alpha = .75) +
  facet_rep_grid(CHRType ~ Diff, scales = "free_y") +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#d9d9d9", "#d9d9d9"), na.translate = FALSE) +
  geom_label(data = fulldf, aes(x = .15, y = .85, label = paste0("# of SNPs: ", scales::comma(SNPs))), alpha = 1,
             size = 4.5, fontface = "bold", fill = "#d6d6d6", label.padding = unit(.5, "lines"), family = "Optima", show.legend = FALSE) +
  geom_label_repel(data = subset(fulldf, CHRType == "Autosomes" & Labels == "Y150239"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .15, nudge_y = 0,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(fulldf, CHRType == "Autosomes" & Labels == "Meerkerk_01"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .09, nudge_y = -.02,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(fulldf, CHRType == "Autosomes" & Labels == "Garderen_01"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .09, nudge_y = .02,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(fulldf, CHRType == "Autosomes" & Labels == "Garderen_02"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .09, nudge_y = 0,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(fulldf, CHRType == "Chromosome Z" & Labels == "Y150239"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .2, nudge_y = .2,
                   point.padding = .6, force_pull = 10, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                                 ends = "last", type = "open")) +
  geom_label_repel(data = subset(fulldf, CHRType == "Chromosome Z" & Labels == "Meerkerk"), aes(x = hybrid.index, y = heterozygosity, label = Labels),
                   family = "Optima", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = 0, nudge_y = .3,
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
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#E5E7E9", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.box = "vertical",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 10, r = 0, l = 0),
        axis.title.x = element_text(family = "Optima", size = 15, face = "bold", margin = margin(t = 18, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Optima", size = 15, face = "bold", margin = margin(t = 0, r = 18, b = 0, l = 0)),
        axis.text.x = element_text(family = "Optima", color = "#000000", size = 11, face = "bold"),
        axis.text.y = element_text(family = "Optima", color = "#000000", size = 11, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 14, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Species", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                                      label.theme = element_text(size = 14, family = "Optima"),
                                      override.aes = list(starshape = 15, size = 5, starstroke = .15), nrow = 1, order = 1),
         starshape = "none",
         colour = "none")


# Saves plot ~
ggsave(Panel, file = "Y150239Genomics--Triangular.Discord.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 10, height = 10, dpi = 600)
ggsave(Panel, file = "Y150239Genomics--Triangular.Discord.png",
       limitsize = FALSE, scale = 1, width = 10, height = 10, dpi = 600)


#
##
### The END ~~~~~