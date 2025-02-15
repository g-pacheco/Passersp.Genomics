### The BEGINNING ~~~~~
##
# ~ Plots Y150239Genomics--evalAdmix | Written by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, ggh4x, ggstar, ggrepel, ggnewscale, lemon, data.table, tidytext, patchwork)


# Defines the orderInds function ~
orderInds <- function(q=NULL, pop=NULL, popord=NULL){
  ordpop <- function(x, pop, q){
    idx <- which(pop==x)
    main_k <- which.max(apply(as.matrix(q[idx,]),2,mean))
    ord <- order(q[idx,main_k])
    idx[ord]} 
  
  if(!is.null(pop)){
    
    if(is.null(popord)) popord <- unique(pop)
    
    if(!is.null(q)){ 
      
      ord <- unlist(sapply(popord, ordpop, pop=pop, q=q))
      
    } else if (is.null(q)) {
      
      ord <- unlist(sapply(popord, function(x) which(pop==x)))
      
    }
  } else if (is.null(pop)&!is.null(q)) {
    
    # Gets index of k with max value per individual
    main_k <- apply(q,1, which.max)
    # Get max q per indivdiual ~
    main_q <- q[cbind(1:nrow(q),main_k)]
    ord <- order(main_k, main_q)
    
  } else {stop("Need at least an argument to order.")}
  return(ord)}


# Defines compute_mean_correlations function ~
compute_mean_correlations <- function(cor_mat_list, ord_list, pop) {
  pop <- pop[ord]
  unique_pops <- unique(pop)
  num_pops <- length(unique_pops)
  process_single_matrix <- function(cor_mat, ord, pop) {
    pop <- pop[ord]
    annotations <- cor_mat[, c("Sample_ID_1", "Population_1", "CHRType", "K")]
    cor_mat <- cor_mat[, !(colnames(cor_mat) %in% c("Sample_ID_1", "Population_1", "CHRType", "K"))]
    mean_cor_df <- data.frame(matrix(ncol = num_pops, nrow = num_pops))
    rownames(mean_cor_df) <- unique_pops
    colnames(mean_cor_df) <- unique_pops
    for (i1 in 1:num_pops) {
      for (i2 in 1:num_pops) {
        p1 <- unique_pops[i1]
        p2 <- unique_pops[i2]
        indices_p1 <- which(pop == p1)
        indices_p2 <- which(pop == p2)
        cor_values <- cor_mat[indices_p1, indices_p2]
        mean_cor_df[p1, p2] <- mean(cor_values[!is.na(cor_values)])}}
    mean_cor_df[is.na(mean_cor_df)] <- 0
    for (i1 in 1:(nrow(cor_mat) - 1)) {
      for (i2 in (i1 + 1):nrow(cor_mat)) {
        cor_mat[i2, i1] <- mean_cor_df[pop[i1], pop[i2]]
        cor_mat[i1, i2] <- cor_mat[i1, i2]}}
    cor_mat <- cbind(annotations, cor_mat)
    return(cor_mat)}
  
  final_list <- list()
  
  for (i in seq_along(cor_mat_list)) {
    cor_mat_list[[i]] <- process_single_matrix(cor_mat_list[[i]], ord_list[[i]], pop)
    current_K <- unique(cor_mat_list[[i]]$K)
    final_list[[as.character(current_K)]] <- cor_mat_list[[i]]}
  return(final_list)}


# Initialize empty lists for storing data separately for Allosome and Autosomes
corres_allosome <- list()
corres_autosomes <- list()
final_list_allosome <- list()
final_list_autosomes <- list()
annot_allosome <- list()
annot_autosomes <- list()
ord_list_allosome <- list()
ord_list_autosomes <- list()

# Define the two folder paths
folder_paths <- c("./Autosomes", "./Allosome")

# Loop over each folder (for Allosome and Autosomes)
for (folder in folder_paths) {
  
  # Get the files in the folder
  corres_files <- dir(folder, pattern = ".corres")
  annot_files <- dir(folder, pattern = ".labels")
  qopt_files <- dir(folder, pattern = ".qopt")
  
  # Check if files exist
  if (length(corres_files) == 0 || length(annot_files) == 0 || length(qopt_files) == 0) {
    stop(paste("Missing files in folder:", folder))
  }
  
  # Process each file
  for (k in seq_along(annot_files)) {
    # Read annotation file
    annot <- read.table(file.path(folder, annot_files[k]), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(annot) <- c("Annot")
    
    # Assign populations based on the annotations
    annot$Population <- ifelse(grepl("FR0", annot$Annot), "Sales",
                        ifelse(grepl("KAZ", annot$Annot), "Chokpak",
                        ifelse(grepl("Lesina", annot$Annot), "Lesina",
                        ifelse(grepl("Crotone", annot$Annot), "Crotone",
                        ifelse(grepl("Guglionesi", annot$Annot), "Guglionesi",
                        ifelse(grepl("PI22NLD0001M", annot$Annot), "Focal Area",
                        ifelse(grepl("PD22NLD0146F", annot$Annot), "Focal Area",
                        ifelse(grepl("PD22NLD0147F", annot$Annot), "Focal Area",
                        ifelse(grepl("PDOM2022NLD0077M", annot$Annot), "Focal Area",
                        ifelse(grepl("PDOM2022NLD0", annot$Annot), "Utrecht", "Error"))))))))))
    annot$Ind <- with(annot, ave(Population, Population, FUN = function(x) sprintf("%s_%02d", x, seq_along(x))))
    
    # Read qopt file if it exists
    qopt_df <- NULL
    if (length(qopt_files) >= k && file.exists(file.path(folder, qopt_files[k]))) {
      qopt_df <- as.matrix(read.table(file.path(folder, qopt_files[k]), header = FALSE))
    } else {
      warning(paste("Missing or empty qopt file:", qopt_files[k], "in folder:", folder))
    }
    
    # Read corres file
    corres_df <- as.data.frame(read.table(file.path(folder, corres_files[k])))
    if (nrow(corres_df) == 0 || ncol(corres_df) == 0) {
      stop(paste("Empty or invalid corres file:", corres_files[k], "in folder:", folder))
    }
    
    labels <- annot$Annot
    pop <- annot$Population
    ord <- orderInds(q = qopt_df, pop = pop)
    
    # Check if ord is valid
    if (length(ord) != nrow(corres_df)) {
      stop(paste("Invalid ordering vector (ord) for file:", corres_files[k], "in folder:", folder))
    }
    
    # Update the correct ord_list and annot_list depending on the CHRType
    if (grepl("Allosome", folder)) {
      ord_list_allosome[[k]] <- ord
      annot_allosome[[k]] <- annot  # Store annot separately for Allosome
      pop_allo <- annot$Population
    } else {
      ord_list_autosomes[[k]] <- ord
      annot_autosomes[[k]] <- annot  # Store annot separately for Autosomes
      pop_auto <- annot$Population
    }
    
    # Reorder corres_df based on ord
    corres_df <- corres_df[ord, ord]
    ordered_labels <- labels[ord]
    rownames(corres_df) <- ordered_labels
    colnames(corres_df) <- ordered_labels
    corres_df$Sample_ID_1 <- rownames(corres_df)
    corres_df$Population_1 <- annot$Population
    corres_df$CHRType <- str_extract(corres_files[k], "(Allosome|Autosomes)")
    corres_df$K <- str_extract(corres_files[k], "(K2|K3|K4|K5|K6|K7)")
    
    # Format K values
    corres_df$K <- ifelse(grepl("K2", corres_df$K), "K = 2",
                   ifelse(grepl("K3", corres_df$K), "K = 3",
                   ifelse(grepl("K4", corres_df$K), "K = 4",
                   ifelse(grepl("K5", corres_df$K), "K = 5",
                   ifelse(grepl("K6", corres_df$K), "K = 6",
                   ifelse(grepl("K7", corres_df$K), "K = 7", "Error"))))))
    
    # Add the data to the corresponding CHRType list
    if (grepl("Allosome", corres_df$CHRType[1])) {
      corres_allosome[[k]] <- corres_df
      # Apply compute_mean_correlations for Allosome
      final_list_allosome[[k]] <- compute_mean_correlations(cor_mat_list = list(corres_df), ord_list = list(ord), pop = pop_allo)
    } else {
      corres_autosomes[[k]] <- corres_df
      # Apply compute_mean_correlations for Autosomes
      final_list_autosomes[[k]] <- compute_mean_correlations(cor_mat_list = list(corres_df), ord_list = list(ord), pop = pop_auto)
    }
  }
}

                                            
# Combines all matrices for different Ks into data frame ~
final_combined_autosomes <- bind_rows(unlist(final_list_autosomes, recursive = FALSE), .id = "K_Value")
final_combined_allosome <- bind_rows(unlist(final_list_allosome, recursive = FALSE), .id = "K_Value")


# Converts data frame into long ~ 
final_long_format_autosomes <- final_combined_autosomes %>%
  pivot_longer(cols = -c(K_Value, Sample_ID_1, Population_1, CHRType, K),
               names_to = "Sample_ID_2",
               values_to = "Value") %>%
  select(Sample_ID_1, Sample_ID_2, Population_1, CHRType, K, Value)

final_long_format_allosome <- final_combined_allosome %>%
  pivot_longer(cols = -c(K_Value, Sample_ID_1, Population_1, CHRType, K),
               names_to = "Sample_ID_2",
               values_to = "Value") %>%
  select(Sample_ID_1, Sample_ID_2, Population_1, CHRType, K, Value)


# Gets Population_2 ~
fulldf_autosomes <- final_long_format_autosomes %>%
  mutate(Population_2 = ifelse(grepl("FR0", Sample_ID_2), "Sales",
                        ifelse(grepl("KAZ", Sample_ID_2), "Chokpak",
                        ifelse(grepl("Lesina", Sample_ID_2), "Lesina",
                        ifelse(grepl("Crotone", Sample_ID_2), "Crotone",
                        ifelse(grepl("Guglionesi", Sample_ID_2), "Guglionesi",
                        ifelse(grepl("PI22NLD0001M", Sample_ID_2), "Focal Area",
                        ifelse(grepl("PD22NLD0146F", Sample_ID_2), "Focal Area",
                        ifelse(grepl("PD22NLD0147F", Sample_ID_2), "Focal Area",
                        ifelse(grepl("PDOM2022NLD0077M", Sample_ID_2), "Focal Area",
                        ifelse(grepl("PDOM2022NLD0", Sample_ID_2), "Utrecht", "Error"))))))))))) %>%
  select(1:3, Population_2, everything())

fulldf_allosome <- final_long_format_allosome %>%
  mutate(Population_2 = ifelse(grepl("FR0", Sample_ID_2), "Sales",
                        ifelse(grepl("KAZ", Sample_ID_2), "Chokpak",
                        ifelse(grepl("Lesina", Sample_ID_2), "Lesina",
                        ifelse(grepl("Crotone", Sample_ID_2), "Crotone",
                        ifelse(grepl("Guglionesi", Sample_ID_2), "Guglionesi",
                        ifelse(grepl("PI22NLD0001M", Sample_ID_2), "Focal Area",
                        ifelse(grepl("PD22NLD0146F", Sample_ID_2), "Focal Area",
                        ifelse(grepl("PD22NLD0147F", Sample_ID_2), "Focal Area",
                        ifelse(grepl("PDOM2022NLD0077M", Sample_ID_2), "Focal Area",
                        ifelse(grepl("PDOM2022NLD0", Sample_ID_2), "Utrecht", "Error"))))))))))) %>%
  select(1:3, Population_2, everything())


# Defines the generate_ordered_permutations function ~ 
generate_ordered_permutations <- function(individuals, k) {
  perm <- do.call(rbind, lapply(individuals, function(id1) {
    data.frame(Sample_ID_1 = id1, Sample_ID_2 = individuals, K = k)}))
  return(perm)}

all_permutations_autosomes <- do.call(rbind, lapply(seq_along(corres_autosomes), function(i) {
  individuals <- corres_autosomes[[i]]$Sample_ID_1
  k <- corres_autosomes[[i]]$K[1]
  generate_ordered_permutations(individuals, k)}))

all_permutations_allosome <- do.call(rbind, lapply(seq_along(corres_allosome), function(i) {
  individuals <- corres_allosome[[i]]$Sample_ID_1
  k <- corres_allosome[[i]]$K[1]
  generate_ordered_permutations(individuals, k)}))


# Sets the Order column per K ~
all_permutations_autosomes <- all_permutations_autosomes %>%
  group_by(K) %>%
  mutate(Order = match(Sample_ID_1, unique(Sample_ID_1))) %>%
  ungroup()

all_permutations_allosome <- all_permutations_allosome %>%
  group_by(K) %>%
  mutate(Order = match(Sample_ID_1, unique(Sample_ID_1))) %>%
  ungroup()


# Defines the reorder_fulldf function ~
reorder_fulldf <- function(df, permutations) {
  permutations$order <- seq_len(nrow(permutations))
  merged <- merge(permutations, df, by = c("Sample_ID_1", "Sample_ID_2", "K"), all.x = TRUE)
  reordered <- merged[order(merged$order), ]
  reordered$order <- NULL
  return(reordered)}


# Splits fulldf and all_permutations by K ~
split_fulldf_autosomes <- split(fulldf_autosomes, fulldf_autosomes$K)
split_permutations_autosomes <- split(all_permutations_autosomes, all_permutations_autosomes$K)

split_fulldf_allosome <- split(fulldf_allosome, fulldf_allosome$K)
split_permutations_allosome <- split(all_permutations_allosome, all_permutations_allosome$K)


# Applies the reordering function to each subset of fulldf ~
fulldfUp_autosomes <- do.call(rbind, lapply(names(split_fulldf_autosomes), function(k) {
  reordered_df <- reorder_fulldf(split_fulldf_autosomes[[k]], split_permutations_autosomes[[k]])
  reordered_df$Order <- split_permutations_autosomes[[k]]$Order
  return(reordered_df)}))

fulldfUp_allosome <- do.call(rbind, lapply(names(split_fulldf_allosome), function(k) {
  reordered_df <- reorder_fulldf(split_fulldf_allosome[[k]], split_permutations_allosome[[k]])
  reordered_df$Order <- split_permutations_allosome[[k]]$Order
  return(reordered_df)}))


# Sets factor levels for Sample_ID_1 & Sample_ID_2 per K ~ 
fulldfUp_autosomes <- fulldfUp_autosomes %>%
  group_by(K) %>%
  mutate(Sample_ID_1 = factor(Sample_ID_1, levels = unique(Sample_ID_1)),
         Sample_ID_2 = factor(Sample_ID_2, levels = unique(Sample_ID_1))) %>%
  ungroup()


fulldfUp_allosome <- fulldfUp_allosome %>%
  group_by(K) %>%
  mutate(Sample_ID_1 = factor(Sample_ID_1, levels = unique(Sample_ID_1)),
         Sample_ID_2 = factor(Sample_ID_2, levels = unique(Sample_ID_1))) %>%
  ungroup()

# Safely convert to numeric factor
fulldfUp_autosomes <- fulldfUp_autosomes %>% mutate(Sample_ID_Factor = as.numeric(Sample_ID_1))
fulldfUp_allosome <- fulldfUp_allosome %>% mutate(Sample_ID_Factor = as.numeric(Sample_ID_1))


fulldfUp <- rbind(fulldfUp_autosomes, fulldfUp_allosome)


# Calculate population positions
population_positions <- fulldfUp %>%
  filter(!is.na(Sample_ID_Factor)) %>%
  group_by(Population_1) %>%
  summarise(center = (min(Sample_ID_Factor) + max(Sample_ID_Factor)) / 2)


# Corrects CHRType ~
levels(fulldfUp$CHRType <- sub("Allosome", "Chromosome Z", fulldfUp$CHRType))


# Reorders CHRType ~
fulldfUp$CHRType <- factor(fulldfUp$CHRType, ordered = TRUE,
                           levels = c("Autosomes",
                                      "Chromosome Z"))


# Function to set the Triangle column ~
assign_triangle <- function(df) {
  df <- df %>%
    mutate(Sample_ID_1 = as.character(Sample_ID_1),
           Sample_ID_2 = as.character(Sample_ID_2),
           Pair = paste0(pmin(Sample_ID_1, Sample_ID_2), "_", pmax(Sample_ID_1, Sample_ID_2)),
           Triangle = case_when(Sample_ID_1 == Sample_ID_2 ~ "Diagonal", !duplicated(Pair) ~ "Individual", duplicated(Pair) ~ "Population")) %>%
    select(-Pair)
  return(df)}


# Applies the assign_triangle function per K ~
fulldfUp <- fulldfUp %>%
  group_by(CHRType, K) %>%
  do({df <- assign_triangle(.)
  df$K <- unique(df$K)
  df}) %>%
  ungroup()


# Creates fulldf_points ~
fulldf_points <- fulldfUp %>%
                 filter((CHRType == "Autosomes" & K %in% c("K = 2", "K = 3", "K = 4")) | (CHRType == "Chromosome Z" & K %in% c("K = 5", "K = 6", "K = 7"))) %>%
                 filter(Triangle %in% c("Population", "Individual")) %>%
                 group_by(CHRType, K, Triangle, Population_1, Population_2) %>%
                 mutate(Value = if_else(Triangle == "Population", round(mean(Value, na.rm = TRUE), 9), Value)) %>%
                 ungroup() %>%
                 #mutate(temp = if_else(str_detect(Population_2, "Focal Region"), Population_1, Population_2),
                 #Population_1 = if_else(str_detect(Population_2, "Focal Region"), Population_2, Population_1),
                 #Population_2 = temp) %>%
                 #select(-temp) %>%
                 mutate(PopPair = paste(Population_1, "Vs", Population_2)) %>%
                 mutate(Labels = paste("Vs", Population_2)) %>%
                 group_by(CHRType, K, Triangle, Population_1, Population_2) %>%
                 filter(if_else(Triangle == "Population", !duplicated(PopPair), TRUE)) %>%
                 ungroup()


# Defines color palette and breaks ~
color_palette <- c("#4575b4", "#EAEDE9", "#d73027")  
nHalf <- 10
Min <- -.15
Max <- .15
Thresh <- 0
rc1 <- colorRampPalette(colors = color_palette[1:2], space = "Lab")(nHalf)
rc2 <- colorRampPalette(colors = color_palette[2:3], space = "Lab")(nHalf)
rampcols <- c(rc1, rc2)
rampcols[c(nHalf, nHalf+1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue = 256) 
rb1 <- seq(Min, Thresh, length.out = nHalf + 1)
rb2 <- seq(Thresh, Max, length.out = nHalf + 1)[-1]
rampbreaks <- c(rb1, rb2)


# Reorders Population_1 in fulldf_points ~
fulldf_points$Population_1 <- factor(fulldf_points$Population_1, ordered = T,
                                     levels = c("Focal Area",
                                                "Chokpak",
                                                "Lesina",
                                                "Guglionesi",
                                                "Crotone",
                                                "Sales",
                                                "Utrecht"))


# Reorders chrtype ~
fulldf_points$K <- factor(fulldf_points$K, ordered = T,
                         levels = c("K = 7",
                                    "K = 6",
                                    "K = 5",
                                    "K = 4",
                                    "K = 3",
                                    "K = 2"))

# Creates the panel ~
Y150239Genomics_evalAdmix_Points_Plot <-
  ggplot(fulldf_points, aes(x = Value, y = Population_1)) +
  geom_vline(xintercept = 0, linewidth = .35, linetype = 4, color = "#000000") +
  geom_violin(data = subset(fulldf_points, Triangle == "Individual"), width = .6, linewidth = .2) + 
  geom_star(data = subset(fulldf_points, Triangle == "Population"), aes(fill = as.numeric(Value)),
            size = 3.25, starshape = 15, alpha = .85, starstroke = .1, color = "#000000") +
  geom_star(data = subset(fulldf_points, CHRType == "Chromosome Z" & Population_1 == "Focal Area" & Triangle == "Individual"),
            size = 3.25, starshape = 15, alpha = .85, starstroke = .25, fill = NA, color = "#000000") +
  geom_label_repel(data = subset(fulldf_points, Triangle == "Population" & Value >= .015), aes(label = Labels),
                   family = "Optima", size = 4.25, fontface = "bold", nudge_x = .06, nudge_y = .5,
                   point.padding = 1, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  geom_label_repel(data = subset(fulldf_points, Triangle == "Population" & Value >= .09), aes(label = Labels),
                   family = "Optima", size = 4.25, fontface = "bold", nudge_y = -1.5,
                   point.padding = 2, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  geom_label_repel(data = subset(fulldf_points, Triangle == "Population" & Value <= -.01), aes(label = Labels),
                   family = "Optima", size = 4.25, fontface = "bold", nudge_x = -.05, nudge_y = .5,
                   point.padding = 1, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  scale_x_continuous(limits = c(-.15, .15),
                     breaks = c(-.1, 0, .1),
                     labels = c(-0.10, 0, 0.10)) +
  scale_fill_gradientn(colors = rampcols,
                       limits = c(-.15, .15),
                       breaks = c(-.15, 0, .15),
                       labels = c(-0.15, 0, 0.15)) +
  facet_nested(CHRType + K ~ Triangle, scales = "free_x", remove_labels = "y", 
               strip = strip_nested(text_x = elem_list_text(size = 16, family = "Optima", face = "bold"),
                                    background_x = elem_list_rect(fill = "#d6d6d6", colour = "#000000", linewidth = .3),
                                    by_layer_x = TRUE,
                                    text_y = elem_list_text(size = c(16, 14), family = c("Optima", "Optima"), face = c("bold", "bold")),
                                    background_y = elem_list_rect(fill = c("#d6d6d6", "#FAFAFA"), colour = c("#000000", "#000000"), linewidth = c(.3, .3)),
                                    by_layer_y = TRUE)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#E5E7E9", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.spacing = unit(.1, "cm"),
        legend.position = "top",
        legend.title.align = .5,
        legend.title = element_text(family = "Optima", size = 16, face = "bold"),
        legend.text = element_text(family = "Optima", size = 10, face = "bold"),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box = "horizontal",
        legend.box.margin = margin(t = 10, b = 5, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title = element_blank(),
        axis.text = element_text(family = "Optima", colour = "#000000", size = 13, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.placement = "outside") +
  guides(fill = guide_colourbar(title = "Residual Correlation", 
                                title.position = "top", title.theme = element_text(family = "Optima", size = 22, face = "bold"),
                                label.theme = element_text(family = "Optima", size = 14, face = "bold"), 
                                label.position = "bottom", barwidth = 12, barheight = 1, order = 1, frame.linetype = 1, 
                                frame.colour = NA, ticks.colour = NA, direction = "horizontal", 
                                even.steps = TRUE, draw.ulim = TRUE, draw.llim = TRUE))


# Saves the panel ~
ggsave(Y150239Genomics_evalAdmix_Points_Plot, file = "Y150239Genomics--evalAdmix_Points.pdf",
       device = cairo_pdf, width = 12, height = 14, scale = 1, dpi = 600)
ggsave(Y150239Genomics_evalAdmix_Points_Plot, file = "Y150239Genomics--evalAdmix_Points.jpeg",
       width = 12, height = 14, scale = 1, dpi = 600)


# Defines plotting function ~ 
plot_for_K <- function(k_value, facet_type = "default", show_x_labels = FALSE, is_last_plot = FALSE) {
  subset_data <- fulldfUp %>%
    filter(K == k_value) %>%
    mutate(
      Sample_ID_1_ordered = factor(Sample_ID_1, levels = unique(Sample_ID_1[order(Order)])),
      Sample_ID_2_ordered = factor(Sample_ID_2, levels = unique(Sample_ID_2[order(Order)])))
  
  facet_formula <- if (facet_type == "CHRType") {as.formula("K ~ CHRType")}
  else {as.formula("K ~ .")}
  
  ggplot(subset_data, aes(x = Sample_ID_1_ordered, y = Sample_ID_2_ordered, fill = as.numeric(Value))) +
    geom_tile(data = subset(subset_data, Triangle == "Population"), linewidth = 0) +
    geom_tile(data = subset(subset_data, Triangle != "Population"), linewidth = .15, colour = "#000000") +
    scale_fill_gradientn(colors = rampcols, na.value = "#d6d6d6", breaks = c(-.3, 0, .3), rampbreaks, limits = c(-.3, .3)) +
    scale_x_discrete(labels = if (show_x_labels) {
      function(x) {
        labels <- rep("", length(x))
        for (i in seq_along(population_positions$center)) {
          labels[population_positions$center[i]] <- population_positions$Population_1[i]
        }
        labels
      }
    } else NULL,
    expand = c(0, 0),
    drop = FALSE) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    facet_grid(facet_formula, scales = "free", space = "free") +
    theme(panel.background = element_rect(fill = "#ffffff"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1, "lines"),
          legend.position = ifelse(is_last_plot, "right", "none"),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
          legend.box = "vertical",
          legend.box.margin = margin(t = 20, b = 30, r = 0, l = 0),
          axis.title = element_blank(),
          axis.text.x = if (show_x_labels) {element_text(color = "#000000", family = "Optima", size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1)}
          else element_blank(),
          axis.text.y = element_text(color = "#000000", family = "Optima", size = 9, face = "bold"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = "#000000", linewidth = 0.15),
          strip.text.x = element_text(colour = "#000000", size = 22, face = "bold", family = "Optima"),
          strip.text.y = element_text(colour = "#000000", size = 18, face = "bold", family = "Optima"),
          strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = 0.15),
          axis.line = element_line(colour = "#000000", linewidth = 0.15)) +
    guides(fill = guide_colourbar(title = "", title.theme = element_text(size = 16, face = "bold"),
                                  label.theme = element_text(size = 10, face = "bold"), label.position = "right",
                                  barwidth = 1.25, barheight = 18, order = 1, frame.linetype = 1, frame.colour = NA,
                                  ticks.colour = NA, direction = "vertical", even.steps = TRUE,
                                  draw.ulim = TRUE, draw.llim = TRUE))}


# Generate the facets for each K ~
unique_K_values <- rev(unique(fulldfUp$K))
plots <- lapply(seq_along(unique_K_values), function(i) {
  facet_type <- if (i == 1) "CHRType" else "default"
  show_x_labels <- (i == length(unique_K_values))
  is_last_plot <- (i == 2)
  plot_for_K(unique_K_values[i], facet_type = facet_type, show_x_labels = show_x_labels, is_last_plot = is_last_plot)})


# Creates panel ~
combined_plot <- wrap_plots(plots, ncol = 1)


# Saves panel ~
ggsave(combined_plot, file = "Y150239Genomics--evalAdmix.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 14, height = 50, dpi = 600)
ggsave(combined_plot, file = "Y150239Genomics--evalAdmix.jpeg",
       limitsize = FALSE, scale = 1, width = 14, height = 20, dpi = 600)


#
##
### The END ~~~~~


#TEST <-
#  fulldfUp |>
#  mutate(Sample_ID_2 = reorder_within(Sample_ID_2, as.numeric(Value), K), linewidth = .15, colour = "#000000") |>
#  ggplot(aes(Sample_ID_1, Sample_ID_2)) +
#  facet_grid(K ~ ., scales = "free", space = "free") +
#  geom_tile(aes(fill = as.numeric(Value)), linewidth = .15, colour = "#000000") +
#  scale_y_reordered() +
#  theme(panel.background = element_rect(fill = "#ffffff"),
#        panel.border = element_blank(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.spacing = unit(1, "lines"),
#        legend.position = "right",
#        legend.key = element_blank(),
#        legend.background = element_blank(),
#        legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
#        legend.box = "vertical",
#        legend.box.margin = margin(t = 20, b = 30, r = 0, l = 0),
#        axis.title = element_blank(),
#        axis.text.x = element_text(color = "#000000", family = "Optima", size = 9, face = "bold"),
#        axis.text.y = element_text(color = "#000000", family = "Optima", size = 9, face = "bold"),
#        axis.ticks.x = element_blank(),
#        axis.ticks.y = element_line(color = "#000000", linewidth = 0.15),
#        strip.text.x = element_text(colour = "#000000", size = 22, face = "bold", family = "Optima"),
#        strip.text.y = element_text(colour = "#000000", size = 18, face = "bold", family = "Optima"),
#        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = 0.15),
#        axis.line = element_line(colour = "#000000", linewidth = 0.15)) +
#  guides(fill = guide_colourbar(title = "", title.theme = element_text(size = 16, face = "bold"),
#                                label.theme = element_text(size = 10, face = "bold"), label.position = "right",
#                               barwidth = 1.25, barheight = 18, order = 1, frame.linetype = 1, frame.colour = NA,
#                                ticks.colour = NA, direction = "vertical", even.steps = TRUE,
#                                draw.ulim = TRUE, draw.llim = TRUE))


# Saves panel ~
#ggsave(TEST, file = "TEST.png",
#       limitsize = FALSE, scale = 1, width = 10, height = 10, dpi = 600)