### The BEGINNING ~~~~~
##
# ~ Plots Passer sp. Genomics -- GlobalDepth | Written by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, extrafont, ggh4x)


# Loads data ~
fulldf <- read.table("AllSamples_bcftools.raw.vcf.Filtered.MeanDepth.ldepth.mean", header = TRUE)
fulldf$Type <- ""


# Expands fulldf by adding chrtype ~
fulldf$chrtype <- ifelse(grepl("chrZ", fulldf$CHROM), "Chromosome Z",
                  ifelse(grepl("mtDNA", fulldf$CHROM), "mtGenome", "Autosomes"))


# Reorders chrtype ~
fulldf$chrtype <- factor(fulldf$chrtype, ordered = T,
                         levels = c("Autosomes",
                                    "Chromosome Z",
                                    "mtGenome"))


# Get quantile values ~
quantiles_df <- fulldf %>%
  group_by(chrtype) %>%
  summarize(x_quantile = quantile(MEAN_DEPTH, .95))


# Expands quantiles_df by adding annotation ~
quantiles_df <- quantiles_df %>%
  add_column(group = c("Autosomes", "Chromosome Z", "mtGenome")) %>%
  add_column(label = c(sprintf("Quantile 95%%: %.0fX", quantiles_df$x_quantile))) %>%
  add_column(value_x = c(62.5, 62.5, 62.5)) %>%
  add_column(v_just = c(2.14, 2.3, 2.95)) %>%
  add_column(Type = "")


# Custom y-axis breaks ~
breaks_fun <- function(y){
  caseVal <- max(y)
  if (caseVal > 1.6){
    c(.4, .8, 1.2, 1.6)}
  else if (caseVal < .1){
    c(.01, .02, .03, .04)}
  else {
    c(.05, .1, .15, .2, .25)}}


# Custom y-axis limits ~
limits_fun <- function(x){
  limitVal <- max(x)
  if (limitVal > 1.6){
    c(0, 1.8)}
  else if (limitVal < .1){
    c(0, .046)}
  else {
    c(0, .259)}}


# Creates the plot ~
GlobalCoverage <-
 ggplot(fulldf, aes(x = MEAN_DEPTH, fill = Type, colour = Type), alpha = .1) +
  geom_density(alpha = .15, linewidth = .3, adjust = 1) +
  facet_grid2(chrtype ~., scales = "free", axes = "all", remove_labels = "x") +
  scale_fill_manual(values = c("#fa9fb5")) +
  scale_colour_manual(values = c("#000000")) +
  geom_label(data = quantiles_df, aes(x = value_x, y = Inf, vjust = v_just, label = label), alpha = .1,
             size = 3.7, fontface = "bold", label.padding = unit(.5, "lines"), family = "Optima", show.legend = FALSE) +
  scale_x_continuous("Mean Depth (X)",
                     breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                     labels = c("10", "20", "30", "40", "50", "60", "70", "80", "90", "100"),
                     limits = c(0, 102),
                     expand = c(0, 0)) +
  scale_y_continuous("Density",
                     breaks = breaks_fun,
                     limits = limits_fun,
                     expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#d9d9d9", linetype = "dashed", linewidth = .05),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.spacing.y = unit(.2, "cm"),
        axis.title.x = element_text(family = "Optima", size = 16, face = "bold", color = "#000000", margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Optima", size = 16, face = "bold", color = "#000000", margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.ticks = element_line(linewidth = .3, color = "#000000"),
        axis.text = element_text(family = "Optima", size = 11, color = "#000000", face = "bold"),
        strip.text = element_text(colour = "#000000", size = 13, face = "bold", family = "Optima"),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3),
        legend.position = "none")

  
# Saves plot ~
ggsave(GlobalCoverage, file = "Passersp.Genomics--MeanDepth.pdf",
       width = 12, height = 12, scale = 1, device = cairo_pdf, dpi = 600)
ggsave(GlobalCoverage, file = "Passersp.Genomics--MeanDepth.jpeg",
       width = 12, height = 12, scale = 1, dpi = 600)


#
##
### The END ~~~~~