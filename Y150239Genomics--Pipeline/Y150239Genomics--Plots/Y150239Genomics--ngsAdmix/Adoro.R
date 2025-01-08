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


# Imports data while incorporating annotation ~
corres <- list()
annot <- list()
corres_files <- dir(pattern = ".corres")
annot_files <- dir(pattern = ".labels")
Sample_Order <- c("Utrecht_01", "Utrecht_02", "Utrecht_03", "Utrecht_04", "Utrecht_05", "Utrecht_06", "Utrecht_07", "Utrecht_08", "Utrecht_09", "Utrecht_10",
                  "Utrecht_11", "Utrecht_12", "Utrecht_13", "Utrecht_14", "Utrecht_15",
                  "Garderen_01", "Garderen_02", "Meerkerk_01",
                  "Sales_01", "Sales_02", "Sales_03", "Sales_04", "Sales_05", "Sales_06", "Sales_07", "Sales_08", "Sales_09", "Sales_10",
                  "Crotone_01", "Crotone_02", "Crotone_03", "Crotone_04", "Crotone_05", "Crotone_06", "Crotone_07", "Crotone_08", "Crotone_09", "Crotone_10",
                  "Guglionesi_01", "Guglionesi_02", "Guglionesi_03", "Guglionesi_04", "Guglionesi_05", "Guglionesi_06", "Guglionesi_07", "Guglionesi_08", "Guglionesi_09", "Guglionesi_10",
                  "Lesina_01", "Lesina_02", "Lesina_03", "Lesina_04", "Lesina_05", "Lesina_06", "Lesina_07", "Lesina_08", "Lesina_09", "Lesina_10",
                  "Chokpak_01", "Chokpak_02", "Chokpak_03", "Chokpak_04", "Chokpak_05", "Chokpak_06", "Chokpak_07", "Chokpak_08", "Chokpak_09", "Chokpak_10",
                  "Y150239_01")

# Remove suffixes from Sample_Order to get Pop_Order
Pop_Order <- gsub("_\\d+$", "", Sample_Order)

for (k in seq_along(annot_files)) {
  annot[[k]] <- read.table(annot_files[k], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(annot[[k]]) <- c("Annot")
  annot[[k]]$Population <- ifelse(grepl("FR0", annot[[k]]$Annot), "Sales",
                                  ifelse(grepl("KAZ", annot[[k]]$Annot), "Chokpak",
                                         ifelse(grepl("Lesina", annot[[k]]$Annot), "Lesina",
                                                ifelse(grepl("Crotone", annot[[k]]$Annot), "Crotone",
                                                       ifelse(grepl("Guglionesi", annot[[k]]$Annot), "Guglionesi",
                                                              ifelse(grepl("PI22NLD0001M", annot[[k]]$Annot), "Y150239",
                                                                     ifelse(grepl("PD22NLD0146F", annot[[k]]$Annot), "Garderen",
                                                                            ifelse(grepl("PD22NLD0147F", annot[[k]]$Annot), "Garderen",
                                                                                   ifelse(grepl("PDOM2022NLD0077M", annot[[k]]$Annot), "Meerkerk",
                                                                                          ifelse(grepl("PDOM2022NLD0", annot[[k]]$Annot), "Utrecht", "Error"))))))))))
  
  
  # Sort annot based on Population_1
  #annot[[k]] <- annot[[k]][order(annot[[k]]$Population_1), ]
  
  # Create Ind_1 without duplicates or NA values
  annot[[k]]$Ind <- with(annot[[k]], ave(Population, Population, FUN = function(x) sprintf("%s_%02d", x, seq_along(x))))
  annot[[k]] <- annot[[k]][!duplicated(annot[[k]]$Ind) & !is.na(annot[[k]]$Ind), ]
  
  # Load the corresponding file
  corres_df <- as.data.frame(read.table(corres_files[k]))
  
  # Ensure rownames and colnames match Sample_Order
  rownames(corres_df) <- annot[[k]]$Ind
  colnames(corres_df) <- annot[[k]]$Ind
  
  # Order corres_df rows and columns according to Sample_Order
  corres_df <- corres_df[match(Sample_Order, rownames(corres_df)), match(Sample_Order, colnames(corres_df))]
  
  # Remove rows and columns that are NA
  corres_df <- corres_df[complete.cases(rownames(corres_df)), complete.cases(colnames(corres_df))]
  
  corres_df$CHRType <- str_extract(corres_files[k], "(Allosome|Autosomes)")
  corres_df$K <- str_extract(corres_files[k], "(K2|K3|K4|K5|K6|K7)")
  corres_df$K <- ifelse(grepl("K2", corres_df$K), "K = 2",
                        ifelse(grepl("K3", corres_df$K), "K = 3",
                               ifelse(grepl("K4", corres_df$K), "K = 4",
                                      ifelse(grepl("K5", corres_df$K), "K = 5",
                                             ifelse(grepl("K6", corres_df$K), "K = 6",
                                                    ifelse(grepl("K7", corres_df$K), "K = 7", "Error"))))))
  corres[[k]] <- corres_df}


# Defines compute_individual_correlations function ~
compute_individual_correlations <- function(corres_ind_list, pop_vector, ord_list) {
  result <- lapply(1:length(corres_ind_list), function(i) {
    corres_ind <- corres_ind_list[[i]]
    annotation_df <- data.frame(Sample_ID_1 = rownames(corres_ind), Population_1 = pop_vector)
    for (pop_name in unique(pop_vector)) {
      if (pop_name %in% names(ord_list[[i]])) {
        pop_ord <- ord_list[[i]][[pop_name]]
        chr_type <- unique(pop_ord$CHRType)
        k_val <- unique(pop_ord$K)
        annotation_df$CHRType[annotation_df$Population_1 == pop_name] <- chr_type
        annotation_df$K[annotation_df$Population_1 == pop_name] <- k_val}}
    indiv_matrix <- cbind(annotation_df, corres_ind)
    return(indiv_matrix)})
  return(result)}


# Prepare the ord_list
ord_list <- lapply(corres, orderInds)


# Applies compute_on_matrices function ~
ind_list <- compute_individual_correlations(corres_ind_list = corres,
                                            pop_vector = Pop_Order,
                                            ord_list = ord_list)


# Defines compute_mean_correlations function ~
compute_mean_correlations <- function(cor_mat_list, pop) {
  unique_pops <- unique(pop)
  num_pops <- length(unique_pops)
  process_single_matrix <- function(cor_mat) {
    # Extract and remove annotation columns
    annotations <- cor_mat[, c("Sample_ID_1", "Population_1", "CHRType", "K")]
    cor_mat <- cor_mat[, !(colnames(cor_mat) %in% c("Sample_ID_1", "Population_1", "CHRType", "K"))]
    mean_cor_df <- data.frame(matrix(ncol = num_pops, nrow = num_pops))
    rownames(mean_cor_df) <- unique_pops
    colnames(mean_cor_df) <- unique_pops
    # Calculate mean correlations for each pair of populations
    for (i1 in 1:num_pops) {
      for (i2 in 1:num_pops) {
        p1 <- unique_pops[i1]
        p2 <- unique_pops[i2]
        indices_p1 <- which(pop == p1)
        indices_p2 <- which(pop == p2)
        cor_values <- cor_mat[indices_p1, indices_p2]
        mean_cor_df[p1, p2] <- mean(cor_values[!is.na(cor_values)])}}
    # Replace NaN values with 0
    mean_cor_df[is.na(mean_cor_df)] <- 0
    # Update cor_mat using mean_cor_df
    for (i1 in 1:(nrow(cor_mat) - 1)) {
      for (i2 in (i1 + 1):nrow(cor_mat)) {
        cor_mat[i1, i2] <- mean_cor_df[pop[i1], pop[i2]]}}
    # Reattach annotation columns
    cor_mat <- cbind(annotations, cor_mat)
    return(cor_mat)}
  # Apply the process_single_matrix function to each matrix in the list and update them in place
  for (i in seq_along(cor_mat_list)) {
    cor_mat_list[[i]] <- process_single_matrix(cor_mat_list[[i]])}
  return(cor_mat_list)}


# Applies compute_mean_correlation function ~
final_list <- compute_mean_correlations(ind_list, Pop_Order)


# Convert all matrices to numeric type
finall_list <- lapply(final_list, as.matrix)
final_list <- lapply(final_list, as.data.frame)


# Combine all matrices into a single dataset
final_combined <- bind_rows(final_list, .id = "Matrix_ID")


# Select relevant columns and convert data from wide to long format
final_long_format <- final_combined %>%
  pivot_longer(cols = -c(Matrix_ID, Sample_ID_1, Population_1, CHRType, K),
               names_to = "Var2",
               values_to = "value")


# Renames some columns and select the relevant ones ~ 
fulldf_final <- final_long_format %>%
  rename_with(~"Sample_ID_2", Var2) %>%
  rename_with(~"Value", value) %>%
  select(Sample_ID_1, Sample_ID_2, Population_1, CHRType, K, Value) %>%
  ungroup()


# Gets Population_2 ~
fulldf_final <- fulldf_final %>%
  mutate(Population_2 = ifelse(grepl("Sales", Sample_ID_2), "Sales",
                               ifelse(grepl("Chokpak", Sample_ID_2), "Chokpak",
                                      ifelse(grepl("Lesina", Sample_ID_2), "Lesina",
                                             ifelse(grepl("Crotone", Sample_ID_2), "Crotone",
                                                    ifelse(grepl("Guglionesi", Sample_ID_2), "Guglionesi",
                                                           ifelse(grepl("Y150239", Sample_ID_2), "Y150239",
                                                                  ifelse(grepl("Garderen", Sample_ID_2), "Garderen",
                                                                         ifelse(grepl("Meerkerk", Sample_ID_2), "Meerkerk",
                                                                                ifelse(grepl("Utrecht", Sample_ID_2), "Utrecht", "Error")))))))))) %>%
  select(1:3, Population_2, everything())


# Reorders columns ~
fulldf <- fulldf_final %>%
  select(Sample_ID_1, Sample_ID_2, Population_1, Population_2, CHRType, K, Value)


# Corrects Population names ~
levels(fulldf$Sample_ID_1 <- sub("Y150239_01", "Y150239", fulldf$Sample_ID_1))
levels(fulldf$Sample_ID_2 <- sub("Y150239_01", "Y150239", fulldf$Sample_ID_2))


# Reorders CHRType ~
fulldf$CHRType <- factor(fulldf$CHRType, ordered = TRUE,
                         levels = c("Autosomes",
                                    "Allosome"))


# Reorders K ~
fulldf$K <- factor(fulldf$K, ordered = TRUE,
                   levels = c("K = 7",
                              "K = 6", 
                              "K = 5",
                              "K = 4",
                              "K = 3",
                              "K = 2"))


fulldf$Sample_ID_1 <- factor(fulldf$Sample_ID_1, levels = fulldf$Sample_ID_1)
#fulldf$Sample_ID_2 <- factor(fulldf$Sample_ID_2, levels = fulldf$Sample_ID_2)

# Defines color palette and breaks ~
color_palette <- c("#023858", "#ffffff", "#a50f15")
nHalf <- 4
Min <- -.3
Max <- .3
Thresh <- 0

rc1 <- colorRampPalette(colors = color_palette[1:2], space = "Lab")(nHalf)
rc2 <- colorRampPalette(colors = color_palette[2:3], space = "Lab")(nHalf)
rampcols <- c(rc1, rc2)
rampcols[c(nHalf, nHalf+1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue = 256) 

rb1 <- seq(Min, Thresh, length.out = nHalf + 1)
rb2 <- seq(Thresh, Max, length.out = nHalf + 1)[-1]
rampbreaks <- c(rb1, rb2)


# Creates heatmap ~
evalAdmix_Plot <- 
  ggplot(data = fulldf) +
  geom_tile(aes(x = Sample_ID_1, y = Sample_ID_2, fill = as.numeric(Value)), linewidth = .15, colour = "#000000") +
  scale_fill_gradientn(colors = rampcols, na.value = "#d6d6d6", breaks = rampbreaks, limits = c(-.3, .3)) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(K ~ CHRType, scales = "free", space = "free") +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
        legend.box = "vertical",
        legend.box.margin = margin(t = 20, b = 30, r = 0, l = 0),
        axis.title = element_blank(),
        axis.text.x = element_text(color = "#000000", family = "Optima", size = 9, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "#000000", family = "Optima", size = 9, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .15),
        strip.text = element_text(colour = "#000000", size = 26, face = "bold", family = "Optima"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .15),
        axis.line = element_line(colour = "#000000", linewidth = .15)) +
  guides(fill = guide_legend(title = "", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15), reverse = TRUE))


# Saves plot (Boxplot) ~
ggsave(evalAdmix_Plot, file = "Y150239Genomics--evalAdmix.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 30, height = 50, dpi = 600)