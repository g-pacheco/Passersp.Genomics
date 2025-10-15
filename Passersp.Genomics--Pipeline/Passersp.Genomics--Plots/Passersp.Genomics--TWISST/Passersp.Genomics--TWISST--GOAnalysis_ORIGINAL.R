### The BEGINNING ~~~~~
##
# Passersp.Genomics--TWISST by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads packages ~
pacman::p_load(tidyverse, scales, reshape2, ggh4x, lemon, patchwork, GenomicRanges, txdbmaker, rtracklayer, GOstats, GSEABase, outliers, clusterProfiler)


# Imports weights ~
Wlgz <- list()
Wlistgz <- dir(pattern = ".gz")
for (k in 1:length(Wlistgz)){
  Wlgz[[k]] <- read.table(gzfile(Wlistgz[k]))[-(1:1), ]
  Wlgz[[k]]$TotalWeight <- as.numeric(Wlgz[[k]]$V1) + as.numeric(Wlgz[[k]]$V2) + as.numeric(Wlgz[[k]]$V3)
  col_prefix <- ifelse(grepl("WithMeerkerk", Wlistgz[k]), "Meerkerk",
                       ifelse(grepl("WithY150239", Wlistgz[k]), "Y150239", "Error"))
  colnames(Wlgz[[k]]) <- c(paste(col_prefix, "~ Spanish"),
                           paste(col_prefix, "~ House"),
                           paste("House ~ Spanish"), "TotalWeight")
  Wlgz[[k]]$CHR <- gsub(".*(WithMeerkerk|WithY150239)\\.Phased\\.MinMaf\\.", "", Wlistgz[k])
  Wlgz[[k]]$CHR <- gsub("\\.SW150\\.Weights\\.csv\\.gz", "", Wlgz[[k]]$CHR)
  Wlgz[[k]]$Phylo <- ifelse(grepl("WithMeerkerk", Wlistgz[k]), "Meerkerk",
                            ifelse(grepl("WithY150239", Wlistgz[k]), "Y150239", "Error"))}


# Subsets list of data frames ~
Y150239_df <- Wlgz[36:70]
Meerkerk_df <- Wlgz[1:35]


# Expands list of data frames 
Y150239_df_WeightsDF <- bind_rows(Y150239_df)
Meerkerk_df_WeightsDF <- bind_rows(Meerkerk_df)


# Selects relevant columns ~
Meerkerk_df_WeightsDF <- dplyr::select(Meerkerk_df_WeightsDF, "Meerkerk ~ House", "Meerkerk ~ Spanish")


# Binds the data frames ~
fulldf <- bind_cols(Meerkerk_df_WeightsDF, Y150239_df_WeightsDF)


# Calculates the percentage of weights & delta ~
fulldf$PercY150239 <- as.numeric(fulldf$"Y150239 ~ Spanish") / (as.numeric(fulldf$"Y150239 ~ House") + as.numeric(fulldf$"Y150239 ~ Spanish")) * 100
fulldf$PercMeerkerk <- as.numeric(fulldf$"Meerkerk ~ Spanish") / (as.numeric(fulldf$"Meerkerk ~ House") + as.numeric(fulldf$"Meerkerk ~ Spanish")) * 100
fulldf$Delta <- (fulldf$PercY150239 - fulldf$PercMeerkerk) / 100


# Reorders columns ~
WeightsDF <- dplyr::select(fulldf, "Phylo", "CHR", "TotalWeight", "Y150239 ~ House", "Y150239 ~ Spanish", "Meerkerk ~ House", "Meerkerk ~ Spanish", "PercY150239", "PercMeerkerk", "Delta")


# Imports windows´ data ~
Wilgz <- list()
Wilistgz <- dir(pattern = ".tsv")
for (k in 1:length(Wilistgz)){
  Wilgz[[k]] <- read.table(Wilistgz[k])[-1, ]
  colnames(Wilgz[[k]]) <- c("Scaffold", "Start", "End", "Mid", "Sites", "lnL")}


# Melts windows ~
WindowsDF <- reshape2::melt(Wilgz)
WindowsDF <- WindowsDF[, -ncol(WindowsDF)]


# Merges WeightsDF & WindowsDF ~
fulldf <- cbind(WeightsDF, WindowsDF)


# Converts DF from wide into long ~
fulldfUp <- gather(fulldf, Estimation, Value, "Y150239 ~ House", "Y150239 ~ Spanish", "Meerkerk ~ House", "Meerkerk ~ Spanish", "Delta")


# Reorders CHR ~
fulldfUp$CHR <- factor(fulldfUp$CHR, ordered = TRUE,
                       levels = c("chr1", "chr1A", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                                  "chr15", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr26", "chr27", "chr28", "chrZ",
                                  "scaffold00169", "scaffold00221", "scaffold00223", "scaffold00224", "scaffold00238", "scaffold00239", "scaffold00242"))


# Corrects the y-strip facet labels ~
y_strip_labels <- c("chr1" = "CHR 01", "chr1A" = "CHR 01A", "chr2" = "CHR 02", "chr3" = "CHR 03", "chr4" = "CHR 04", "chr5" = "CHR 05", "chr6" = "CHR 06", "chr7" = "CHR 07",
                    "chr8" = "CHR 08", "chr9" = "CHR 09", "chr10" = "CHR 10", "chr11" = "CHR 11", "chr12" = "CHR 12", "chr13" = "CHR 13", "chr14" = "CHR 14", "chr15" = "CHR 15",
                    "chr17" = "CHR 17", "chr18" = "CHR 18", "chr19" = "CHR 19", "chr20" = "CHR 20", "chr21" = "CHR 21", "chr22" = "CHR 22", "chr23" = "CHR 23",
                    "chr24" = "CHR 24", "chr26" = "CHR 26", "chr27" = "CHR 27", "chr28" = "CHR 28", "chrZ" = "CHR Z",
                    "scaffold00169" = "Scaffold00169", "scaffold00221" = "Scaffold00221", "scaffold00223" = "Scaffold00223", "scaffold00224" = "Scaffold00224",
                    "scaffold00238" = "Scaffold00238", "scaffold00239" = "Scaffold00239", "scaffold00242" = "Scaffold00242")


generate_dynamic_breaks_and_labels <- function(min_val, max_val) {
  if (max_val <= 0) {return(list(breaks = c(1), labels = c("1Mb")))}
  
  # Determines appropriate step size based on data range ~
  data_range <- max_val - min_val
  if (data_range <= 1300000) {
    step_size <- 100000
  } else if (data_range <= 10000000) {
    step_size <- 1000000
  } else if (data_range <= 100000000) {
    step_size <- 5000000
  } else {
    step_size <- 20000000}
  
  # Generates breaks starting from the rounded min_val ~
  rounded_min_val <- floor(min_val / step_size) * step_size
  breaks <- seq(from = rounded_min_val, to = max_val, by = step_size)
  
  # Ensures breaks are within the actual data range ~
  breaks <- breaks[breaks >= min_val & breaks <= max_val]
  
  # Generates labels corresponding to the breaks ~
  if (step_size >= 1e6) {
    labels <- paste0(breaks / 1e6, "Mb")
  } else {
    labels <- paste0(round(breaks / 1e5) / 10, "Mb")}
  
  # Ensures breaks and labels have the same length ~
  if (length(breaks) != length(labels)) {
    stop("Breaks and labels have different lengths.")}
  return(list(breaks = breaks, labels = labels))}


# Applies the function to each unique CHR ~
patterns <- unique(fulldfUp$CHR)


# Creates an empty lists to store the filtered results ~ 
filtered_positions_upper <- list()
filtered_positions_lower <- list()
filtered_positions_outliers <- list()


# Sets function ~ 
for (x in patterns) {
  subset_df <- subset(fulldfUp, CHR == x)
  max_mid <- max(as.numeric(subset_df$Mid), na.rm = TRUE)
  min_mid <- min(as.numeric(subset_df$Mid), na.rm = TRUE)
  breaks_and_labels <- generate_dynamic_breaks_and_labels(min_mid, max_mid)
  Mean_Delta <- mean(as.numeric(subset_df$Value)[subset_df$Estimation == "Delta" & as.numeric(subset_df$Value) >= 0], na.rm = TRUE)
  delta_subset <- subset(subset_df, Estimation == "Delta" & !is.na(Value))
  if (nrow(delta_subset) > 0) {
    delta_values <- as.numeric(delta_subset$Value)
    Q1 <- quantile(delta_values, .25, na.rm = TRUE)
    Q3 <- quantile(delta_values, .75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_fence <- Q1 - 1.5 * IQR
    upper_fence <- Q3 + 1.5 * IQR
    upper_subset <- subset(delta_subset, delta_values > upper_fence)
    lower_subset <- subset(delta_subset, delta_values < lower_fence)
    outlier_subset <- subset(delta_subset, delta_values < lower_fence | delta_values > upper_fence)
    for (df in list(upper_subset, lower_subset, outlier_subset)) {
      df$Start <- as.numeric(df$Mid) - 1
      df$End <- as.numeric(df$Mid)
      df$Start <- as.numeric(df$Start)
      df$End <- as.numeric(df$End)
      df <- df[, c("CHR", "Start", "End", "Value")]}
    filtered_positions_upper[[x]] <- upper_subset
    filtered_positions_lower[[x]] <- lower_subset
    filtered_positions_outliers[[x]] <- outlier_subset}
  
  
  # Combines filtered positions into separate data frames ~
  filtered_positions_upper_df <- do.call(rbind, filtered_positions_upper)
  filtered_positions_lower_df <- do.call(rbind, filtered_positions_lower)
  filtered_positions_outliers_df <- do.call(rbind, filtered_positions_outliers)
  
  
  # Saves the BED file with the filtered positions ~
  write.table(filtered_positions_outliers_df, "Y150239Genomics--TWISST_Outliers.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(filtered_positions_upper_df, "Y150239Genomics--TWISST_Upper.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(filtered_positions_lower_df, "Y150239Genomics--TWISST_Lower.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  # Creates Y150239 plot ~
  Y150239_Plot <- 
    ggplot() +
    geom_area(data = subset_df %>% filter(Estimation %in% c("Y150239 ~ House", "Y150239 ~ Spanish")),
              aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Estimation, group = Estimation),
              position = "fill", colour = "#000000", alpha = .3, linetype = 1, linewidth = .2) +
    scale_x_continuous("Genomic Position",
                       breaks = breaks_and_labels$breaks, 
                       labels = breaks_and_labels$labels,
                       limits = c(min_mid, max_mid + 1000),
                       expand = c(0, 0)) +
    scale_y_continuous("Weights",
                       breaks = c(.25, .50, .75), 
                       labels = c("25%", "50%", "75%"),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
    theme(panel.background = element_rect(fill = "#ffffff"),
          panel.border = element_blank(),
          panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(family = "Optima", size = 36, face = "bold", color = "#000000", margin = margin(t = 0, r = 40, b = 0, l = 15)),
          axis.text.x = element_blank(),
          axis.text.y = element_text(family = "Optima", size = 24, colour = "#000000", face = "bold"),
          axis.ticks = element_line(color = "#000000", linewidth = .3),
          axis.line.x = element_line(colour = "#000000", linewidth = .3),
          axis.line.y = element_line(colour = "#000000", linewidth = .3),
          legend.position = "top",
          legend.title = element_text(margin = margin(r = 20)),
          legend.text = element_text(margin = margin(r = 15, l = 15)),
          legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
          legend.box.margin = margin(t = 15, b = 40, r = 0, l = 0),
          legend.key = element_rect(fill = NA),
          legend.background = element_blank()) +
    guides(colour = guide_legend(title = "Topologies", title.theme = element_text(family = "Optima", size = 34, face = "bold"),
                                 label.theme = element_text(family = "Optima", size = 30), override.aes = list(linewidth = 1.75, linetype = 1)),
           fill = guide_legend(title = "Topologies", title.theme = element_text(family = "Optima", size = 34, face = "bold"),
                               label.theme = element_text(family = "Optima", size = 30), override.aes = list(linewidth = .3, linetype = 1)))
  
  
  # Create Meerkerk plot ~
  Meerkerk_Plot <- ggplot() +
    geom_area(data = subset_df %>% filter(Estimation %in% c("Meerkerk ~ House", "Meerkerk ~ Spanish")),
              aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Estimation, group = Estimation),
              position = "fill", colour = "#000000", alpha = .3, linetype = 1, linewidth = .2) +
    scale_x_continuous("Genomic Position",
                       breaks = breaks_and_labels$breaks,
                       labels = breaks_and_labels$labels,
                       limits = c(min_mid, max_mid + 1000),
                       expand = c(0, 0)) +
    scale_y_continuous("Weights",
                       breaks = c(.25, .50, .75), 
                       labels = c("25%", "50%", "75%"),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
    theme(panel.background = element_rect(fill = "#ffffff"),
          panel.border = element_blank(),
          panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(family = "Optima", size = 36, face = "bold", color = "#000000", margin = margin(t = 0, r = 40, b = 0, l = 15)),
          axis.text.x = element_blank(),
          axis.text.y = element_text(family = "Optima", size = 24, colour = "#000000", face = "bold"),
          axis.ticks = element_line(color = "#000000", linewidth = .3),
          axis.line.x = element_line(colour = "#000000", linewidth = .3),
          axis.line.y = element_line(colour = "#000000", linewidth = .3),
          legend.position = "top",
          legend.title = element_text(margin = margin(r = 20)),
          legend.text = element_text(margin = margin(r = 15, l = 15)),
          legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
          legend.box.margin = margin(t = 15, b = 40, r = 0, l = 0),
          legend.key = element_rect(fill = NA),
          legend.background = element_blank()) +
    guides(colour = guide_legend(title = "Topologies", title.theme = element_text(family = "Optima", size = 34, face = "bold"),
                                 label.theme = element_text(family = "Optima", size = 30), override.aes = list(linewidth = 1.75, linetype = 1)),
           fill = guide_legend(title = "Topologies", title.theme = element_text(family = "Optima", size = 34, face = "bold"),
                               label.theme = element_text(family = "Optima", size = 30), override.aes = list(linewidth = .3, linetype = 1)))
  
  
  # Create Delta plot ~
  Delta_Plot <- 
    ggplot() +
    geom_line(data = subset_df %>% filter(Estimation %in% c("Delta")), aes(x = as.numeric(Mid), y = as.numeric(Value)), colour = "#000000",
              position = "identity", linetype = 1, linewidth = .2) +
    geom_hline(yintercept = Mean_Delta, linetype = "twodash", color = "#33a02c", linewidth = 1) +
    scale_x_continuous("Genomic Position",
                       breaks = breaks_and_labels$breaks, 
                       labels = breaks_and_labels$labels,
                       limits = c(min_mid, max_mid + 1000),
                       expand = c(0, 0)) +
    scale_y_continuous("Delta",
                       limits = c(-1.1, 1.1),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
    theme(panel.background = element_rect(fill = "#ffffff"),
          panel.border = element_blank(),
          panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(family = "Optima", size = 36, face = "bold", color = "#000000", margin = margin(t = 40, r = 0, b = 15, l = 0)),
          axis.title.y = element_text(family = "Optima", size = 36, face = "bold", color = "#000000", margin = margin(t = 0, r = 40, b = 0, l = 15)),
          axis.text = element_text(family = "Optima", size = 24, colour = "#000000", face = "bold"),
          axis.ticks = element_line(color = "#000000", linewidth = .3),
          axis.line = element_line(colour = "#000000", linewidth = .3))
  
  
  # Arranges the plots into a single panel ~
  Panel_Plot <- Y150239_Plot + Meerkerk_Plot + Delta_Plot + plot_layout(ncol = 1)}
  
  
  # Save the panel ~
#ggsave(paste("Y150239Genomics--TWISST_SW150_FillArea_", x, ".pdf", sep = ""), plot = Panel_Plot,
#       device = cairo_pdf, limitsize = FALSE, width = 40, height = 25, scale = 1, dpi = 600)}
#ggsave(paste("Y150239Genomics--TWISST_FillArea_", x, ".jpeg", sep = ""), plot = Panel_Plot,
#       limitsize = FALSE, width = 40, height = 25, scale = 1, dpi = 600)}


# Creates Delta boxplots per CHR ~
Delta_Boxplots <- 
  ggplot(data = fulldfUp %>% filter(Estimation == "Delta"), 
         aes(x = CHR, y = as.numeric(Value)), fill = "#000000") +
  geom_boxplot(outlier.shape = 21, alpha = .7) +
  scale_y_continuous("Delta", 
                     labels = scales::percent_format(scale = 1), 
                     expand = c(0.02, 0)) +
  scale_x_discrete("Chromosome", labels = y_strip_labels) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#dddddd", linetype = "dashed", linewidth = .0005),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(family = "Optima", face = "bold", size = 18),
        axis.text.x = element_text(family = "Optima", size = 12, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(family = "Optima", size = 12, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


# Save Delta boxplots as a PDF ~
#ggsave(plot = Delta_Boxplots, "Y150239Genomics--TWISST_Delta_Boxplots_Per_CHR_B.pdf",
#      device = cairo_pdf, limitsize = FALSE, width = 18, height = 8, dpi = 600)
#ggsave(plot = Delta_Boxplots, "Y150239Genomics--TWISST_Delta_Boxplots_Per_CHR.jpeg",
#       limitsize = FALSE, width = 18, height = 8, dpi = 600)


###################################################################################################################################################################################
###################################################################################################################################################################################


# Imports the House Sparrow annotation ~
HouseGFF <- import("house_sparrow.gff")


# Loads GOTerm table ~
GOTerms <- read.delim("./GOTerms/house_sparrow_genome_assembly-18-11-14_masked.Protein_gffreads.fasta.Edited.tsv", header = FALSE, sep = "\t", col.names = c("Gene_ID", "GO_Term"))
GOTermsOrtho <- read.delim("./GOTerms/Orthogroups.tsv", header = TRUE, sep = "\t", col.names = c("Orthogroup", "Gene_ID_Zebra", "Gene_ID_House"))


# Edits GOTerms lightly ~
GOTerms <- GOTerms |> 
  mutate(GO_Term = GO_Term |> 
           str_replace_all("\\(InterPro\\)", "") |> 
           str_replace_all("\\|", ", "),
         Gene_ID = str_replace_all(Gene_ID, "-.*", ""))


# Creates GOTermsWithData  ~
GOTermsWithData <- GOTerms %>%
  filter(GO_Term != "-") %>%
  separate_rows(GO_Term, sep = ", ") %>%
  mutate(Evidence = "IEA") %>%
  dplyr::select(GO_Term, Evidence, Gene_ID)


# Defines categories
categories <- c("upper", "lower", "outliers")


# Initializes lists ~ 
Regions_GR_list <- list()
GeneOverlaps_list <- list()
GenesInRange_df_list <- list()
GenesInRangeOnlyGenes_df_list <- list()
FocalGenes_list <- list()
GOTermsOrtho_Edited_list <- list()
GO_Params_list <- list()
GO_Enrich_list <- list()
GO_Enrich_Top50_list <- list()


for (cat in categories) {df_name <- paste0("filtered_positions_", cat, "_df")
Regions_GR_list[[cat]] <- GRanges(seqnames = get(df_name)$CHR,
                                  ranges = IRanges(start = as.numeric(get(df_name)$Start),
                                                   end = as.numeric(get(df_name)$End)))

# Finds overlaps ~
GeneOverlaps_list[[cat]] <- findOverlaps(HouseGFF, Regions_GR_list[[cat]])


# Extracts genes ~
GenesInRange <- HouseGFF[queryHits(GeneOverlaps_list[[cat]])]
GenesInRange_df_list[[cat]] <- data.frame(
  Chromosome = as.character(seqnames(GenesInRange)),
  Start = start(GenesInRange),
  End = end(GenesInRange),
  Gene_ID = mcols(GenesInRange)$ID,
  Gene_Name = mcols(GenesInRange)$Name,
  Type = mcols(GenesInRange)$type)
GenesInRange_df_list[[cat]]$Gene_ID <- sub("-.*", "", GenesInRange_df_list[[cat]]$Gene_ID)


# Filters for genes only ~
GenesInRangeOnlyGenes_df_list[[cat]] <- GenesInRange_df_list[[cat]] %>%
  filter(Type == "gene") %>%
  unique()


# Extracts focal genes ~
FocalGenes_list[[cat]] <- as.data.frame(GenesInRangeOnlyGenes_df_list[[cat]]$Gene_ID)
colnames(FocalGenes_list[[cat]]) <- "Gene_ID"


# Saves focal genes ~
saveRDS(FocalGenes_list, file = "../Passersp.Genomics--GOAnalysis/TWISST_FocalGenes_list.rds")}


###################################################################################################################################################################################
###################################################################################################################################################################################


# Creates plot ~
Weights_Area_Real_Plot <-
  ggplot() +
  geom_area(data = fulldfUp, aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Weight, group = Weight),
            position = "identity", colour = "#000000", alpha = .3, linetype = 1, linewidth = .2) +
  facet_grid(CHR ~., scales = "free_x", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     #limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     #breaks = c(50000, 100000, 150000, 200000), 
                     #labels = c("50K", "100K", "150K", "200K"),
                     #limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Area_Real_Plot, file = "Y150239Genomics--TWISST_AreaReal.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 40, height = 30, scale = 1, dpi = 600)
#ggsave(Weights_Area_Real_Plot, file = "Meerkerkgenomics--TWISST_AreaReal.jpeg",
#       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 600)

# Creates plot ~
Weights_Area_Real_Smooth_Plot <-
  fulldfUp %>%
  ggplot(aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Weight, group = Weight)) +
  stat_smooth(geom = 'area', method = 'loess', position = "identity", linetype = 1, linewidth = .2, colour = "#000000", span = .15, alpha = .3) +
  facet_grid(CHR ~., scales = "free_y", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     #breaks = c(50000, 100000, 150000, 200000), 
                     #labels = c("50K", "100K", "150K", "200K"),
                     #limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Area_Real_Smooth_Plot, file = "Meerkerkgenomics--TWISST_WithMeerkerk_AreaRealSmooth.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 40, height = 30, scale = 1, dpi = 600)
#ggsave(Weights_Area_Real_Smooth_Plot, file = "Meerkerkgenomics--TWISST_AreaRealSmooth.jpeg",
#       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 600)


# Creates plot ~
Weights_Line_Plot <-
  ggplot() +
  geom_line(data = fulldfUp, aes(x = as.numeric(Mid), y = as.numeric(Value), colour = Weight, group = Weight),
            position = "identity", linetype = 1, linewidth = .2) +
  facet_rep_grid(CHR ~., scales = "free_y", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     #breaks = c(50000, 100000, 150000, 200000), 
                     #labels = c("50K", "100K", "150K", "200K"),
                     #limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Line_Plot, file = "Meerkerkgenomics--TWISST_WithMeerkerk_LineReal.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 40, height = 30, scale = 1, dpi = 600)
#ggsave(Weights_Line_Plot, file = "Meerkerkgenomics--TWISST_LineReal.jpeg",
#       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 600)


#
##
### The END ~~~~~