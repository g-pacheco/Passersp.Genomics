### The BEGINNING ~~~~~
##
# Plots Passer sp. Genomics -- Map | Written by George Pacheco ~


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(rnaturalearth, rnaturalearthdata, ggtext, sf, ggspatial, tidyverse, ggforce, ggstar,
               ggrepel, gridExtra, grid, cowplot, patchwork, ggpubr)


# Imports .shp files ~
Global <- ne_countries(scale = 'large', returnclass = 'sf')
PDnr <- read_sf(dsn = ".", layer = "Passer_domesticus")
PInr <- read_sf(dsn = ".", layer = "Passer_italiae")
PHnr <- read_sf(dsn = ".", layer = "Passer_hispaniolensis")


# Loads coordinates ~
SamplingEffort <- read.csv2("Passersp.Genomics--Locations.txt", sep = "\t", header = TRUE, encoding = "UTF-8")
SamplingEffort$Longitude <- as.numeric(SamplingEffort$Longitude)
SamplingEffort$Latitude <- as.numeric(SamplingEffort$Latitude)


# Transforms coordinates ~
SamplingEffort_sf <- st_as_sf(SamplingEffort, coords = c("Longitude", "Latitude"), crs = 4326)


# Reorganises the data ~
SamplingEffort_sf$Species <- factor(SamplingEffort_sf$Species,
                                    levels = c("House", "Italian", "Spanish", NA))


# Reorders Population ~
SamplingEffort$Location <- factor(SamplingEffort$Location, ordered = T,
                                  levels = c("Utrecht",
                                             "Sales",
                                             "Crotone",
                                             "Guglionesi",
                                             "Lesina",
                                             "Chokpak",
                                             NA))


# Defines the shapes to be used for each Group ~
shapes.legend <- as.vector(c(1, 9, 12, 28, 11, 23))


# Creates base map ~
MyLegend_Plot <-
 ggplot() +
  geom_star(data = SamplingEffort, aes(x = Longitude, y = Latitude, starshape = Location, fill = Species),
            size = 4.5, alpha = .9, starstroke = .15) +
  scale_starshape_manual(values = shapes.legend, na.translate = FALSE) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#c994c7"), na.translate = FALSE,  drop = FALSE) +
  scale_colour_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#c994c7"), na.translate = FALSE, drop = FALSE) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#E5E7E9", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0, b = 0, r = .2, l = .2, unit = "cm"),
        panel.border = element_rect(colour = "#000000", linewidth = .5, fill = NA),
        panel.spacing = unit(.2, "cm"),
        legend.position = "top",
        legend.box = "vertical",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 0, b = -15, r = 0, l = 0),
        axis.title = element_blank(),
        axis.text = element_text(size = 11, color = "#000000", face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Species", title.theme = element_text(size = 18, face = "bold"),
                             label.theme = element_text(size = 16,),
                             override.aes = list(starshape = 15, size = 5, starstroke = .15), nrow = 1, order = 1),
         starshape = guide_legend(title = "Locations", title.theme = element_text(size = 16, face = "bold"),
                                  label.theme = element_text(size = 14,),
                                  override.aes = list(size = 5, starstroke = .15), nrow = 1, order = 2),
         colour = "none")


# Saves map ~
ggsave(MyLegend_Plot, file = "Passersp.Genomics--Map.pdf", device = cairo_pdf,
       width = 14, height = 8, scale = 1, limitsize = FALSE, dpi = 600)


# Defines the shapes to be used for each Group ~
shapes.auto <- as.vector(c(1, 9, 12, 28, 11, 23, 15))


# Creates base map ~
MapBody <-
  ggplot() +
  #geom_sf(data = PDnr, fill = "#1E90FF", alpha = .2, color = NA) +
  
  #geom_sf(data = PHnr, fill = "#ee0000", alpha = .2, color = NA) +
  geom_sf(data = Global, fill = "#f2f2f2", color = "#000000") +
  geom_sf(data = PInr, fill = "#FFD700", alpha = .5, color = NA) +
  #geom_sf(data = SamplingEffort_sf, aes(fill = Species), size = 4.5, alpha = .9, show.legend = "point", shape = 21, stroke = .3, colour = "#000000") +
  geom_star(data = SamplingEffort, aes(x = Longitude, y = Latitude, starshape = Location, fill = Species),
            size = 4.5, alpha = .9, starstroke = .15) +
  scale_starshape_manual(values = shapes.legend, na.translate = FALSE) +
  coord_sf(xlim = c(-13.6, 75), ylim = c(35, 61), expand = FALSE) +
  geom_label_repel(data = subset(SamplingEffort, Location == "Utrecht"), aes(x = Longitude, y = Latitude, label = Labels, fill = Species), alpha = .9, show.legend = FALSE,
                   size = 4.5, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = -1.5, nudge_y = 2.5,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#c994c7"), na.translate = FALSE,  drop = FALSE) +
  scale_colour_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#c994c7"), na.translate = FALSE, drop = FALSE) +
  scale_x_continuous(breaks = seq(-120, 125, by = 10)) +
  scale_y_continuous(breaks = seq(-20, 70, by = 10)) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering,
                         height = unit(1.5, "cm"), width = unit(1.5, "cm"),) +
  annotation_scale(location = 'tr', line_width = 1.25, text_cex = 1.2, style = "ticks", pad_y = unit(2, "cm")) +
  theme(panel.background = element_rect(fill = "#b3cde3"),
        panel.grid.major = element_line(color = "#ffffff", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = .2, b = .2, r = .2, l = .2, unit = "cm"),
        panel.border = element_rect(colour = "black", linewidth = .5, fill = NA),
        panel.spacing = unit(.2, "cm"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 11, color = "#000000", face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


# Isolates legend ~
MyLegendBlog <- get_legend(MyLegend_Plot)


# Gets final plot ~
MapFull <- ggarrange(MapBody, legend.grob = MyLegendBlog)


# Saves map ~
ggsave(MapFull, file = "Passersp.Genomics--Map.pdf", device = cairo_pdf,
       width = 14, height = 8, scale = 1, limitsize = FALSE, dpi = 600)
#ggsave(Map, file = "Passersp.Genomics--Map.jpeg",
#       width = 12, height = 12, scale = 1, limitsize = FALSE, dpi = 600)


#geom_label_repel(data = subset(SamplingEffort, Location == "Sales"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
#                 size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = -2.5, nudge_y = 1.5,
#                 arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
#geom_label_repel(data = subset(SamplingEffort, Location == "Crotone"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
#                 size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = 1, nudge_y = -2.5,
#                arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
#geom_label_repel(data = subset(SamplingEffort, Location == "Guglionesi"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
#                 size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = -2.25, nudge_y = -2.5,
#                 arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
#geom_label_repel(data = subset(SamplingEffort, Location == "Lesina"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
#                 size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = 4, nudge_y = 2,
#                 arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
#geom_label_repel(data = subset(SamplingEffort, Location == "Chokpak"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
#                 size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = 0, nudge_y = 3,
#                 arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +


#
##
### The END ~~~~~