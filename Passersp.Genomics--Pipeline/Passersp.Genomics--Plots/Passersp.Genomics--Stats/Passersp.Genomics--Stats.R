### The BEGINNING ~~~~~
##
# ~ Plots Passer sp. Genomics -- Stats | Written by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, extrafont, ggh4x)


# Loads datasets ~
Stats <- read.table("Passersp.Genomics--Stats.txt", sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
Adaptors <- read.table("Passersp.Genomics--Stats_Adaptors.txt", sep = "\t", header = FALSE,
                       stringsAsFactors = FALSE); colnames(Adaptors) <- c("Sample_ID", "reads_adaptors")


# Combines DFs ~
fulldf <- merge(Stats, Adaptors, by = "Sample_ID")


# Gets total_reads ~
fulldf$total_reads <-
       fulldf$seq_reads_pairs * 2


# Gets percentages ~
fulldf$percentage_retained_reads <-
       fulldf$seq_retained_reads * 100 / fulldf$total_reads
fulldf$reads_adaptors <-
  fulldf$reads_adaptors * 100 / fulldf$total_reads


# Fixes percentages ~
fulldf$hits_raw_frac <- fulldf$hits_raw_frac * 100
fulldf$hits_clonality <- fulldf$hits_clonality  * 100
fulldf$hits_unique_frac <- fulldf$hits_unique_frac * 100


# Expands PCA_Annot by adding Population ~
fulldf$Population <- ifelse(grepl("FR0", fulldf$Sample_ID), "Sales",
                     ifelse(grepl("KAZ", fulldf$Sample_ID), "Chokpak",
                     ifelse(grepl("Lesina", fulldf$Sample_ID), "Lesina",
                     ifelse(grepl("Crotone", fulldf$Sample_ID), "Crotone",
                     ifelse(grepl("Guglionesi", fulldf$Sample_ID), "Guglionesi",
                     ifelse(grepl("PI22NLD0001M", fulldf$Sample_ID), "Focal Ind.",
                     ifelse(grepl("PD22NLD0146F", fulldf$Sample_ID), "Garderen",
                     ifelse(grepl("PD22NLD0147F", fulldf$Sample_ID), "Garderen",
                     ifelse(grepl("PDOM2022NLD0077M", fulldf$Sample_ID), "Meerkerk",
                     ifelse(grepl("PDOM2022NLD0", fulldf$Sample_ID), "Utrecht",
                     ifelse(grepl("TreeSparrow", fulldf$Sample_ID), "Tree Sparrow", "Error")))))))))))


# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                            levels = c("Utrecht",
                                       "Sales",
                                       "Crotone",
                                       "Guglionesi",
                                       "Lesina",
                                       "Chokpak",
                                       "Focal Ind.",
                                       "Garderen",
                                       "Meerkerk",
                                       "Tree Sparrow"))

# Cleans DF ~
fulldf <- fulldf %>% select(Population, total_reads, percentage_retained_reads,
                            hits_unique_frac, hits_coverage)


# Converts DF from wide into long ~
fulldfUp <- gather(fulldf, Stat, Value, "total_reads", "percentage_retained_reads",
                   "hits_unique_frac", "hits_coverage")


# Reorders Stat ~
fulldfUp$Stat <- factor(fulldfUp$Stat, ordered = T,
                        levels = c("hits_coverage",
                                   "hits_unique_frac",
                                   #"hits_clonality",
                                   #"hits_raw_frac",
                                   "percentage_retained_reads",
                                   #"reads_adaptors",
                                   "total_reads"))


fulldfUp %>% filter(Stat == "hits_coverage") %>%
             summarise(mean(Value))


# Corrects facet labels ~
ylabels <- c("total_reads" = "# of Reads",
             #"reads_adaptors" = "% of Reads With Adaptors",
             "percentage_retained_reads" = "% of Reads Retained",
             #"hits_raw_frac" = "% of Mapped Reads",
             #"hits_clonality" = "% of Clonality",
             "hits_unique_frac" = "% of Uniquely Mapped Reads",
             "hits_coverage" = "Mean Depth")


# Custom y-axis limits ~
limits_fun <- function(x){
  limitVal <- max(x)
  print(x)
  if (limitVal < 60){
    c(0, 60)}
  else if (limitVal < 100){
    c(90, 100)}
  else { 
    c(70000000, 227746600)}}


# Custom y-axis breaks ~
breaks_fun <- function(y){
  caseVal <- min(y)
  print(y)
  if (caseVal < 100 & caseVal > 50){
    seq(20, 60, by = 10)}
  else if (caseVal < 100){
    seq(92, 100, by = 2)}}


# Custom y-axis labels ~
plot_index_labels <- 0
labels_fun <- function(z) {
  plot_index_labels <<- plot_index_labels + 1L
  switch(plot_index_labels,
         scales::label_number(accuracy = 1, suffix = "X")(z),
         #scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         #scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         #scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         scales::label_number(accuracy = 1, scale = 1/1000000, big.mark = "", suffix = "M")(z))}


# Creates the panel ~
Y150239Genomics_Stat <- 
 ggplot() +
  geom_boxplot(data = fulldfUp, aes(x = Population, y = Value),
               outlier.shape = NA, width = .5, lwd = .25, colour = "#000000", fill = "#C19EBE", alpha = .7) +
  scale_y_continuous(labels = labels_fun) +
  facet_grid2(Stat ~ ., scales = "free", axes = "all", remove_labels = "x", labeller = labeller(Stat = ylabels)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#E5E7E9", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.spacing = unit(.2, "cm"),
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title = element_blank(),
        axis.text.x = element_text(family = "Optima", colour = "#000000", size = 11, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(family = "Optima", color = "#000000", size = 11, face = "bold"),
        axis.ticks.x = element_line(color = "#000000", linewidth = .3),
        axis.ticks.y = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 13, face = "bold"),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank())


# Saves the panel ~
ggsave(Y150239Genomics_Stat, file = "Passersp.Genomics--Stats.pdf",
       device = cairo_pdf, width = 12, height = 12, scale = 1, dpi = 600)
ggsave(Y150239Genomics_Stat, file = "Passersp.Genomics--Stats.jpeg",
       width = 12, height = 12, scale = 1, dpi = 600)


#
## 
### The END ~~~~~