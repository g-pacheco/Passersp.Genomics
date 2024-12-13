### The BEGINNING ~~~~~
##
# Y150239Genomics--Kinship by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(optparse, tidyverse, plyr, RColorBrewer, extrafont, ggforce, ggstar, ggrepel, RcppCNPy, reshape2, lemon,
               gridExtra, grid, ggpubr, rphylopic, viridis)


# Loads data while loading accompanying annotations ~
annot <- list()
rab <- list()
annotL <- dir(pattern = ".labels")
rabL <- dir(pattern = ".res")

for (k in 1:length(annotL)) {
  annot[[k]] <- read.table(annotL[k], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(annot[[k]]) <- c("Ind1")
  annot[[k]]$Ind2 <- c(annot[[k]]$Ind1[-c(1)], NA)
  annot[[k]]$Population <- gsub("[A-z._]*(Autosomes|Allosome).", "", annotL[k])
  annot[[k]]$Population <- gsub(".labels", "", annot[[k]]$Population)
  annot[[k]]$SubPopulation <- annot[[k]]$Population
  
  ### Only needed because the population Utrecht has subpopulation ~
  annot[[k]] <- annot[[k]] %>%
    mutate(SubPopulation = case_when(Ind1 %in% c("PD22NLD0146F", "PD22NLD0147F") ~ "Garderen",
                                     Ind1 == "PDOM2022NLD0077M" ~ "Meerkerk",
                                     Ind1 == "PI22NLD0001M" ~ "Y150239", TRUE ~ Population))
  annot[[k]]$CHRType <- str_extract(annotL[k], "(Allosome|Autosomes)")
  annot[[k]] <- annot[[k]] %>%
    group_by(SubPopulation) %>% mutate(Ind1b = paste(SubPopulation, sprintf("%02d", row_number()), sep = "_")) %>%
    ungroup() %>%
    mutate(Ind2b = lead(Ind1b, default = paste0(SubPopulation[1], "_", sprintf("%02d", n() + 1))))
  rab[[k]] <- read.table(rabL[k], header = TRUE, stringsAsFactors = FALSE)
  rab[[k]]$a <- rab[[k]]$a + 1
  rab[[k]]$Population <- annot[[k]]$Population[match(rab[[k]]$a, seq_along(annot[[k]]$Ind1))]
  rab[[k]]$CHRType <- str_extract(rabL[k], "(Allosome|Autosomes)")
  rab[[k]]$Ind1 <- annot[[k]]$Ind1b[match(rab[[k]]$a, seq_along(annot[[k]]$Ind1))]
  rab[[k]]$Ind2 <- annot[[k]]$Ind2b[match(rab[[k]]$b, seq_along(annot[[k]]$Ind1))]
  rab[[k]]$a <- annot[[k]]$Ind1[match(rab[[k]]$a, seq_along(annot[[k]]$Ind1))]
  rab[[k]]$b <- annot[[k]]$Ind2[match(rab[[k]]$b, seq_along(annot[[k]]$Ind1))]
  rab[[k]]$Pair <- paste(rab[[k]]$Ind1, "Vs", rab[[k]]$Ind2)
  
  ### Only needed because of Meerkerk_01 is amidst the Utrecht subpopulation ~
  rab[[k]] <- rab[[k]] %>% mutate(Invert = str_detect(Ind2, "Meerkerk") & !Ind1 %in% c("Garderen_01", "Garderen_02"),
                                  Temp_Ind1 = Ind1,
                                  Temp_Ind2 = Ind2,
                                  Ind1 = ifelse(Invert, Temp_Ind2, Temp_Ind1),
                                  Ind2 = ifelse(Invert, Temp_Ind1, Temp_Ind2),
                                  Pair = paste(Ind1, "Vs", Ind2)) %>%
                                  select(-Invert, -Temp_Ind1, -Temp_Ind2)
  tryCatch({rab[[k]] <- rab[[k]] %>%
            select(Ind1, Ind2, a, b, rab, Pair, CHRType, Population)}, error = function(e) {
            print(paste("Error in rab[[", k, "]] column selection:", e$message))})
  rab[[k]] <- rab[[k]] %>%
              arrange(Ind1, Ind2)}


# Expands rab list ~
fulldf <- bind_rows(rab)


# Fixes individual Y150239 ~
fulldf <- fulldf %>% mutate(across(c(Ind1, Ind2, Pair), ~ str_replace_all(., "Y150239_01", "Y150239")))


# Corrects Population names ~
levels(fulldf$Population <- sub("TreeSparrow", "Tree Sparrow", fulldf$Population))


# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                            levels = c("Utrecht",
                                       "Sales",
                                       "Crotone",
                                       "Guglionesi", 
                                       "Lesina",
                                       "Chokpak",
                                       "Tree Sparrow"))


# Creates plot (Heatmap) ~
Kinship_Plot_Heatmap <-
  ggplot(subset(fulldf, CHRType == "Autosomes"), aes(Ind1, y = Ind2, fill = rab)) + 
  geom_tile(colour = "#000000") +
  scale_fill_continuous(low = "#ffffff", high = "#f768a1") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  facet_wrap(Population ~., scales = "free", ncol = 2) +
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
        axis.text.x = element_text(family = "Optima", color = "#000000", size = 8.25, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(family = "Optima", color = "#000000", size = 8.25, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .2),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 14, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .2),
        axis.line = element_line(colour = "#000000", linewidth = .2)) +
  guides(fill = guide_legend(title = "Rab", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                             label.theme = element_text(family = "Optima", size = 15), reverse = TRUE))


# Saves plot (Heatmap) ~
ggsave(Kinship_Plot_Heatmap, file = "Y150239Genomics--Kinship.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 12, height = 14, dpi = 600)
ggsave(Kinship_Plot_Heatmap, file = "Y150239Genomics--Kinship.png",
       limitsize = FALSE, scale = 1, width = 12, height = 14, dpi = 600)


# Creates plot (Boxplot) ~
Kinship_Plot_Boxplot <-
  ggplot(fulldf, aes(x = Population, y = rab)) +
  geom_boxplot(fill = "#ffffff", colour = "#000000", show.legend = FALSE, linewidth = .285, width = .15, fatten = 1,
               outlier.shape = 21, outlier.fill = "#005824", outlier.colour = "#000000", outlier.stroke = .2, outlier.alpha = .9) +
  geom_point(data = subset(fulldf, rab >= .125), shape = 21, fill = "#df65b0", colour = "#000000", stroke = .2, alpha = .9) +
  geom_hline(yintercept = .125, linetype = "twodash", color = "#df65b0", linewidth = .3, alpha = .9) +
  geom_label_repel(data = subset(fulldf, rab >= .125), aes(label = Pair),
                   size = 3, nudge_x = .1, family = "Optima", fontface = "bold") +
  scale_x_discrete(expand = c(.05, .05)) +
  scale_y_continuous ("Rab", 
                      breaks = c(.1, .2, .3, .4, .5), 
                      labels = c("0.1", "0.2", "0.3", "0.4", "0.5"),
                      limits = c(0, .52),
                      expand = c(.005, .005)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
        legend.box = "vertical",
        legend.box.margin = margin(t = 20, b = 30, r = 0, l = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Optima", color = "#000000", size = 12, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(family = "Optima", color = "#000000", size = 10, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(family = "Optima", color = "#000000", size = 9, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(family = "Optima", colour = "#000000", size = 14, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Rab", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                             label.theme = element_text(family = "Optima", size = 15)))


# Saves plot (Boxplot) ~
ggsave(Kinship_Plot_Boxplot, file = "Y150239Genomics--Kinship_Boxplot.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 8, height = 5, dpi = 600)
ggsave(Kinship_Plot_Boxplot, file = "Y150239Genomics--Kinship_Boxplot.png",
       limitsize = FALSE, scale = 1, width = 8, height = 5, dpi = 600)

#
##
### The END ~~~~~