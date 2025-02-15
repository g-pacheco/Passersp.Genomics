### The BEGINNING ~~~~~
##
# ~ Plots Y150239Genomics--evalAdmix | Written by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, ggnewscale, data.table, tidytext, patchwork)


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


# Imports and process data
corres <- list()
annot <- list()
ord_list <- list()
corres_files <- dir(pattern = ".corres")
annot_files <- dir(pattern = ".labels")
qopt_files <- dir(pattern = ".qopt")


for (k in seq_along(annot_files)) {
  annot[[k]] <- read.table(annot_files[k], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(annot[[k]]) <- c("Annot")
  annot[[k]]$Population <- ifelse(grepl("FR0", annot[[k]]$Annot), "Sales",
                           ifelse(grepl("KAZ", annot[[k]]$Annot), "Chokpak",
                           ifelse(grepl("Lesina", annot[[k]]$Annot), "Lesina",
                           ifelse(grepl("Crotone", annot[[k]]$Annot), "Crotone",
                           ifelse(grepl("Guglionesi", annot[[k]]$Annot), "Guglionesi",
                           ifelse(grepl("PI22NLD0001M", annot[[k]]$Annot), "Y150239",
                           ifelse(grepl("PDOM2022NLD0077M", annot[[k]]$Annot), "Meerkerk",
                           ifelse(grepl("PDOM2022NLD0", annot[[k]]$Annot), "Utrecht", "Error"))))))))
  annot[[k]]$Ind <- with(annot[[k]], ave(Population, Population, FUN = function(x) sprintf("%s_%02d", x, seq_along(x))))
  
  qopt_df <- NULL
  if (length(qopt_files) >= k && file.exists(qopt_files[k])) {
    qopt_df <- as.matrix(read.table(qopt_files[k], header = FALSE))}
  corres_df <- as.data.frame(read.table(corres_files[k]))
  labels <- annot[[k]]$Annot
  pop <- annot[[k]]$Population
  ord <- orderInds(q = qopt_df, pop = pop)
  ord_list[[k]] <- ord
  corres_df <- corres_df[ord, ord]
  ordered_labels <- labels[ord]
  rownames(corres_df) <- ordered_labels
  colnames(corres_df) <- ordered_labels
  corres_df$Sample_ID_1 <- rownames(corres_df)
  corres_df$Population_1 <- annot[[k]]$Population
  corres_df$CHRType <- str_extract(corres_files[k], "(Allosome|Autosomes)")
  corres_df$K <- str_extract(corres_files[k], "(K2|K3|K4|K5|K6|K7)")
  corres_df$K <- ifelse(grepl("K2", corres_df$K), "K = 2",
                 ifelse(grepl("K3", corres_df$K), "K = 3",
                 ifelse(grepl("K4", corres_df$K), "K = 4",
                 ifelse(grepl("K5", corres_df$K), "K = 5",
                 ifelse(grepl("K6", corres_df$K), "K = 6",
                 ifelse(grepl("K7", corres_df$K), "K = 7", "Error"))))))
  corres[[k]] <- corres_df}


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


# Applies the compute_mean_correlations function ~
final_list <- compute_mean_correlations(corres, ord_list, pop)


# Combines all matrices for different Ks into data frame ~
final_combined <- bind_rows(final_list, .id = "K_Value")


# Converts data frame into long ~ 
final_long_format <- final_combined %>%
  pivot_longer(cols = -c(K_Value, Sample_ID_1, Population_1, CHRType, K),
               names_to = "Sample_ID_2",
               values_to = "Value") %>%
  select(Sample_ID_1, Sample_ID_2, Population_1, CHRType, K, Value)


# Gets Population_2 ~
fulldf <- final_long_format %>%
  mutate(Population_2 = ifelse(grepl("FR0", Sample_ID_2), "Sales",
                        ifelse(grepl("KAZ", Sample_ID_2), "Chokpak",
                        ifelse(grepl("Lesina", Sample_ID_2), "Lesina",
                        ifelse(grepl("Crotone", Sample_ID_2), "Crotone",
                        ifelse(grepl("Guglionesi", Sample_ID_2), "Guglionesi",
                        ifelse(grepl("PI22NLD0001M", Sample_ID_2), "Y150239",
                        ifelse(grepl("PDOM2022NLD0077M", Sample_ID_2), "Meerkerk",
                        ifelse(grepl("PDOM2022NLD0", Sample_ID_2), "Utrecht", "Error"))))))))) %>%
  select(1:3, Population_2, everything())


# Defines the generate_ordered_permutations function ~ 
generate_ordered_permutations <- function(individuals, k) {
  perm <- do.call(rbind, lapply(individuals, function(id1) {
    data.frame(Sample_ID_1 = id1, Sample_ID_2 = individuals, K = k)}))
  return(perm)}

all_permutations <- do.call(rbind, lapply(seq_along(corres), function(i) {
  individuals <- corres[[i]]$Sample_ID_1
  k <- corres[[i]]$K[1]
  generate_ordered_permutations(individuals, k)}))


# Sets the Order column per K ~
all_permutations <- all_permutations %>%
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
split_fulldf <- split(fulldf, fulldf$K)
split_permutations <- split(all_permutations, all_permutations$K)


# Applies the reordering function to each subset of fulldf ~
fulldfUp <- do.call(rbind, lapply(names(split_fulldf), function(k) {
  reordered_df <- reorder_fulldf(split_fulldf[[k]], split_permutations[[k]])
  reordered_df$Order <- split_permutations[[k]]$Order
  return(reordered_df)}))


# Sets factor levels for Sample_ID_1 & Sample_ID_2 per K ~ 
fulldfUp <- fulldfUp %>%
  group_by(K) %>%
  mutate(
    Sample_ID_1 = factor(Sample_ID_1, levels = unique(Sample_ID_1)),
    Sample_ID_2 = factor(Sample_ID_2, levels = unique(Sample_ID_1))) %>%
  ungroup()


# Defines color palette and breaks ~
color_palette <- c("#023858", "#ffffff", "#a50f15")
nHalf <- 10
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


# Safely convert to numeric factor
fulldfUp <- fulldfUp %>% mutate(Sample_ID_Factor = as.numeric(Sample_ID_1))


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
  group_by(K) %>%
  do({df <- assign_triangle(.)
  df$K <- unique(df$K)
  df}) %>%
  ungroup()


TEST <-
  fulldfUp |>
  mutate(Sample_ID_2 = reorder_within(Sample_ID_2, as.numeric(Value), K), linewidth = .15, colour = "#000000") |>
  ggplot(aes(Sample_ID_1, Sample_ID_2)) +
  facet_grid(K ~ ., scales = "free", space = "free") +
  geom_tile(aes(fill = as.numeric(Value)), linewidth = .15, colour = "#000000") +
  scale_y_reordered() +
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
        axis.text.x = element_text(color = "#000000", family = "Optima", size = 9, face = "bold"),
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
                                draw.ulim = TRUE, draw.llim = TRUE))


# Saves panel ~
ggsave(TEST, file = "TEST.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 10, height = 20, dpi = 600)
  
  

ggplot(fulldfUp, aes(x = Sample_ID_1, y = Sample_ID_2, fill = as.numeric(Value))) +
  geom_tile(data = subset(subset_data, Triangle == "Population"), linewidth = 0) +
  geom_tile(data = subset(subset_data, Triangle != "Population"), linewidth = .15, colour = "#000000") +
  scale_fill_gradientn(colors = rampcols, na.value = "#d6d6d6", breaks = c(-.3, 0, .3), rampbreaks, limits = c(-.3, .3)) +
  scale_x_discrete(labels = if (show_x_labels) {
    function(x) {labels <- rep("", length(x))
      for (i in seq_along(population_positions$center)) {
        labels[population_positions$center[i]] <- population_positions$Population_1[i]}
      labels}} else NULL, expand = c(0, 0), drop = FALSE) +
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
        axis.text.x = element_blank(),
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
                                draw.ulim = TRUE, draw.llim = TRUE))


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
    geom_tile(data = subset(subset_data, Triangle != "Population"), linewidth = 0.15, colour = "#000000") +
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
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 14, height = 20, dpi = 600)
ggsave(combined_plot, file = "Y150239Genomics--evalAdmix.png",
       limitsize = FALSE, scale = 1, width = 14, height = 20, dpi = 600)


#
##
### The END ~~~~~