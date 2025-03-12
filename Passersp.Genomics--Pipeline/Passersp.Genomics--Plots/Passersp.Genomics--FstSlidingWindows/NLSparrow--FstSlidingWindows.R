### The BEGINNING ~~~~~
##
# ~ Plots NLSparrow--FstSlidingWindows | By George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, extrafont, lemon)


# Imports extra fonts ~
loadfonts(device = "win", quiet = TRUE)


# Loads datasets ~
Target_Guglionesi <- read.table("NLSparrow_MinInd90_10K-10K_Guglionesi.Target--FstWindows.tsv", header = FALSE)
Target_Sales <- read.table("NLSparrow_MinInd90_10K-10K_Sales.Target--FstWindows.tsv", header = FALSE)
Companion_Guglionesi <- read.table("NLSparrow_MinInd90_10K-10K_Guglionesi.Companion--FstWindows.tsv", header = FALSE)
Companion_Sales <- read.table("NLSparrow_MinInd90_10K-10K_Sales.Companion--FstWindows.tsv", header = FALSE)
Guglionesi_Sales <- read.table("NLSparrow_MinInd90_10K-10K_Guglionesi.Sales--FstWindows.tsv", header = FALSE)


# Adds column names to MissingData ~
colnames(Target_Guglionesi) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(Target_Sales) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(Companion_Guglionesi) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(Companion_Sales) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
colnames(Guglionesi_Sales) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")


# Adds Pops column to each DF ~
Target_Guglionesi$Pops <- factor(paste("Target Vs Guglionesi"))
Target_Sales$Pops <- factor(paste("Target Vs Sales"))
Companion_Guglionesi$Pops <- factor(paste("Companion Vs Guglionesi"))
Companion_Sales$Pops <- factor(paste("Companion Vs Sales"))
Guglionesi_Sales$Pops <- factor(paste("Guglionesi Vs Sales"))


# Gets column names ~
Fst_Window_ColumnNames <- colnames(Target_Guglionesi)


# Merges DFs ~
fulldfT <- full_join(Target_Guglionesi, Target_Sales, by = Fst_Window_ColumnNames)
fulldfC <- full_join(Companion_Guglionesi, Companion_Sales, by = Fst_Window_ColumnNames)
fulldfS <- full_join(Target_Guglionesi, Guglionesi_Sales, by = Fst_Window_ColumnNames)
#fulldf <- full_join(fulldf0, Guglionesi_Sales, by = Fst_Window_ColumnNames)


# Reorders Species ~
fulldfT$Pops <- factor(fulldfT$Pops, ordered = T,
                       levels = c("Target Vs Guglionesi",
                                  "Target Vs Sales"))
                                  #"Guglionesi Vs Sales"

fulldfC$Pops <- factor(fulldfC$Pops, ordered = T,
                       levels = c("Companion Vs Guglionesi",
                                  "Companion Vs Sales"))

fulldfS$Pops <- factor(fulldfS$Pops, ordered = T,
                       levels = c("Target Vs Guglionesi",
                                  "Guglionesi Vs Sales"))
#"Guglionesi Vs Sales"


# Reorders Species ~
fulldfT$CHR <- factor(fulldfT$CHR, ordered = T,
                      levels = c("chr1", "chr1A", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                                 "chr18", "chr19", "chr20", "chr21", "chr22", "chrLGE22", "chr23", "chr24", "chr25",
                                 "chr26", "chr27", "chr28", "chrZ"))

fulldfS$CHR <- factor(fulldfS$CHR, ordered = T,
                      levels = c("chr1", "chr1A", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                                 "chr18", "chr19", "chr20", "chr21", "chr22", "chrLGE22", "chr23", "chr24", "chr25",
                                 "chr26", "chr27", "chr28", "chrZ"))
#"Guglionesi Vs Sales"


fulldfTUp <- fulldfT %>% filter(Fst > 0)
fulldfCUp <- fulldfC %>% filter(Fst > 0)
fulldfSUp <- fulldfS %>% filter(Fst > 0)


# Corrects the y-strip facet labels ~
y_strip_labels <- c("chr1" = "CHR 01", "chr1A" = "CHR 01A", "chr2" = "CHR 02", "chr3" = "CHR 03", "chr4" = "CHR 04",
                    "chr5" = "CHR 05", "chr6" = "CHR 06", "chr7" = "CHR 07", "chr8" = "CHR 08",
                    "chr9" = "CHR 09", "chr10" = "CHR 10", "chr11" = "CHR 11", "chr12" = "CHR 12",
                    "chr13" = "CHR 13", "chr14" = "CHR 14", "chr15" = "CHR 15", "chr16" = "CHR 16",
                    "chr17" = "CHR 17", "chr18" = "CHR 18", "chr19" = "CHR 19", "chr20" = "CHR 20",
                    "chr21" = "CHR 21", "chr22" = "CHR 22", "chrLGE22" = "CHR LEG22", "chr23" = "CHR 23",
                    "chr24" = "CHR 24", "chr25" = "CHR 25", "chr26" = "CHR 26", "chr27" = "CHR 27",
                    "chr28" = "CHR 28", "chrZ" = "CHR Z")



#CHR_09 <- fulldf %>% filter(CHR == "NC_045160.1" & Pops == "North Sea Vs Baltic Sea" & gPoint > 15000000 & gPoint < 20000000)


# Creates Manhattan panel ~
Fst_Window <-
  ggplot() +
  geom_line(data = fulldfSUp, aes(x = gPoint, y = Fst, colour = Pops), linetype = 1, linewidth = .4) +
  facet_rep_grid(CHR ~. , scales = "free", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000), 
                     labels = c("10Mb", "20Mb", "30Mb", "40Mb", "50Mb", "60Mb", "70Mb", "80Mb"),
                     limits = c(0, 80500000),
                     expand = c(.005, .005)) +
  scale_y_continuous("Fst Across Chrmosomes",
                     breaks = c(.25, .5, .75), 
                     labels = c(".25", ".5", ".75"),
                     limits = c(0, 1),
                     expand = c(.01, .01)) +
  scale_colour_manual(values = c("#d95f02", "#1b9e77", "#4daf4a")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text.x = element_text(colour = "#000000", size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "#000000", size = 12, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(colour = guide_legend(title = "Fst Comparisons", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = "none")


# Saves Manhattan plot ~
ggsave(Fst_Window, file = "NLSparrow_MinInd90_10K-10K--FstSlidingWindows.Sales.pdf",
       device = cairo_pdf, scale = 1, width = 28, height = 42, dpi = 100)


# Creates Manhattan panel ~
CHR_09_Plot <-
  ggplot() +
  geom_point(data = CHR_09, aes(x = gPoint, y = log(Fst), colour = Pops), size = .35) +
  facet_rep_grid(CHR ~. , scales = "free", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(17000000, 17100000, 17200000, 17300000, 17400000, 17500000), 
                     labels = c("17Mb", "17.1Mb", "17.2Mb", "17.3Mb", "17.4Mb", "17.5Mb"),
                     limits = c(17000000, 17500000),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous("Fst Across Chrmosomes",
                     breaks = c(.25, .5, .75), 
                     labels = c(".25", ".50", ".75"),
                     limits = c(0, .8),
                     expand = c(0.01, 0.01)) +
  scale_colour_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 15),
        axis.ticks = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(colour = guide_legend(title = "Fst Comparisons:", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(size = 1.4)),
         fill = "none")


# Saves Manhattan plot ~
ggsave(CHR_09_Plot, file = "BSG_AtlanticHerring_CHR_09--Fst-Windows.pdf",
       device = cairo_pdf, scale = 1, width = 26, height = 16, dpi = 600)


# Converts DF from wide into long ~
fulldfUp <- gather(fulldf, Estimate, Value, "Fst")
fulldfUp <- fulldfUp %>% filter(Value >= "0")


# Creates Histogram plot ~
Histo <-  
  ggplot() +
  geom_histogram(data = fulldfUp, aes(x = Value, fill = Pops), colour = "#000000", alpha = .5, bins = 200) +
  coord_trans(y = "log1p") +
  scale_x_continuous("5K-Window Fst",
                     breaks = c(.2, .4, .6, .8, 1), 
                     labels = c("0.2", "0.4", "0.6", "0.8", "1.0"),
                     limits = c(0, 1),
                     expand = c(0.005, 0.005)) +
  scale_y_continuous("Frequency (log1p)",
                     breaks = c(10, 100, 1000, 10000), 
                     labels = c("10", "100", "1000", "10000"),
                     limits = c(0, 40000),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "#000000"),
        legend.position = "top",
        legend.background = element_rect(fill = FALSE),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold", color = "#000000", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, face = "bold", color = "#000000", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = 11, color = "#000000", face = "bold"),
        axis.ticks = element_line(size = .3, color = "#000000"),
        axis.line = element_line(colour = "#000000", size = .3, color = "#000000")) +
  guides(fill = guide_legend(title = "Comparisons", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))



# Saves Histogram plot ~
ggsave(Histo, file = "BSG_AtlanticHerring_Ind66_5K_BalticSea--Fst_Histogram_Log.pdf",
       device = cairo_pdf, width = 4, height = 2, scale = 4, dpi = 600)
ggsave(Histo, file = "BSG_AtlanticHerring_Ind66_5K_BalticSea--Fst_Histogram_Log.jpg",
       width = 4, height = 2, scale = 4, dpi = 600)


#
##
### The END ~~~~~