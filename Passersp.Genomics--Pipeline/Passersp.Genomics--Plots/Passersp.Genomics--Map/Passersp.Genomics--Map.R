### The BEGINNING ~~~~~
##
# Plots Passer sp. Genomics -- Map | Written by George Pacheco ~


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(rnaturalearth, rnaturalearthdata, sf, ggspatial, tidyverse,
               ggrepel, extrafont, cowplot, gridExtra, patchwork, jpeg, png, grid, magick)


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
SamplingEffort_sf$Species <- factor(SamplingEffort_sf$Species, levels = c("House", "Italian", "Spanish"))


# Creates base map ~
Map <-
 ggplot() +
  geom_sf(data = PDnr, fill = "#1E90FF", alpha = .15, color = NA) +
  geom_sf(data = PInr, fill = "#FFD700", alpha = .15, color = NA) +
  geom_sf(data = PHnr, fill = "#ee0000", alpha = .15, color = NA) +
  geom_sf(data = Global, fill = "#ffffff", alpha = .25, color = "#000000") +
  geom_sf(data = SamplingEffort_sf, aes(fill = Species), size = 4.5, alpha = .75, show.legend = "point", shape = 21, colour = "#000000") +
  coord_sf(xlim = c(-13.6, 75), ylim = c(35, 61), expand = FALSE) +
  geom_label_repel(data = subset(SamplingEffort, Location == "Utrecht"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
                   family = "Optima", size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = -1.5, nudge_y = 2.5,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  geom_label_repel(data = subset(SamplingEffort, Location == "Sales"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
                   family = "Optima", size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = -2.5, nudge_y = 1.5,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  geom_label_repel(data = subset(SamplingEffort, Location == "Crotone"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
                   family = "Optima", size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = 1, nudge_y = -2.5,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  geom_label_repel(data = subset(SamplingEffort, Location == "Guglionesi"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
                   family = "Optima", size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = -2.25, nudge_y = -2.5,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  geom_label_repel(data = subset(SamplingEffort, Location == "Lesina"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
                   family = "Optima", size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = 4, nudge_y = 2,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  geom_label_repel(data = subset(SamplingEffort, Location == "Chokpak"), aes(x = Longitude, y = Latitude, label = Location, fill = Species), alpha = .9, show.legend = FALSE,
                   family = "Optima", size = 4, fontface = "bold", colour = "#252525", point.padding = 2, segment.size = .3, nudge_x = 0, nudge_y = 3,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#c994c7"), drop = FALSE) +
  scale_colour_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#c994c7"), drop = FALSE) +
  scale_x_continuous(breaks = seq(-120, 125, by = 10)) +
  scale_y_continuous(breaks = seq(-20, 70, by = 10)) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering,
                         height = unit(2, "cm"), width = unit(2, "cm"),) +
  annotation_scale(location = 'tr', line_width = 1.25, text_cex = 1.2, style = "ticks", pad_y = unit(2, "cm"),) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#E5E7E9", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = .5, fill = NA),
        panel.spacing.y = unit(.2, "cm"),
        legend.position = "top",
        legend.box = "vertical",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 10, r = 0, l = 0),
        axis.title = element_blank(),
        axis.text = element_text(family = "Optima", size = 11, color = "#000000", face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Species", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                             label.theme = element_text(size = 14, family = "Optima"),
                             override.aes = list(starshape = 15, size = 5, starstroke = .15), nrow = 1, order = 1),
         colour = "none")


# Saves map ~
ggsave(Map, file = "Passersp.Genomics--Map.pdf", device = cairo_pdf,
       width = 14, height = 8, scale = 1, limitsize = FALSE, dpi = 600)
ggsave(Map, file = "Passersp.Genomics--Map.jpeg",
       width = 12, height = 12, scale = 1, limitsize = FALSE, dpi = 600)


# adding image to graph 
Oi <- Map +                  
             inset_element(p = PDimg,
             left = .5,
             bottom = .5,
             right = .5,
             top = .5)


ggdraw() +
  draw_image(file.path("ItalianSparrow.pdf")) +
  draw_plot(Map)


ggdraw() +
  draw_plot(Map)
  draw_image("Passer.png")


# Map2 - Faroe Islands ~
Map2 <-
 ggplot() + 
  geom_sf(data = FRO, fill = "#ffffff", color = "black", alpha = .5) +
  geom_sf(data = NR[NR$origin == "1",], fill = "#d4b9da", alpha = .35, color = NA) +
  geom_sf(data = Coords_FRO_sf, aes(fill = Class_Article), size = 5, alpha = 1, show.legend = FALSE, shape = 21, colour = "black") +
  #geom_label_repel(data = Coords_FRO, size = 4.5, seed = 10, min.segment.length = 0, force = 30, segment.curvature = 1, nudge_x = 0, nudge_y = 0, max.overlaps = Inf,
  #                 fontface = "bold", colour = "black", aes(x = Longitude, y = Latitude, label = LocationOnly,
  #                 fill = Class_Article), alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = "Faroe Islands"), x = -7.65, y = 62.335, size = 7.5, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("#44AA99"), drop = FALSE) +
  scale_colour_manual(values = c("#44AA99"), drop = FALSE) +
  scale_x_continuous(breaks = seq(-8, -4, by = .5)) +
  scale_y_continuous(breaks = seq(59, 63.0, by = .5)) +
  coord_sf(xlim = c(-8.1, -6.185), ylim = c(61.35, 62.42), expand = FALSE) +
  annotation_north_arrow(location = "bl", which_north = "false", style = north_arrow_fancy_orienteering,
                         pad_x = unit(.05, "in"), pad_y = unit(.15, "in")) +
  annotation_scale(location = 'bl', line_width = 1.25, text_cex = 1.2, style = "ticks") +
  theme(panel.background = element_rect(fill = "#f7fbff"),
        panel.border = element_rect(colour = "#a50026", linewidth = .5, linetype = "dotdash", fill = NA),
        panel.grid.major = element_line(color = "#d9d9d9", linetype = "dashed", linewidth = 0.00005),
        plot.margin =  margin(t = 0, b = 0, r = .2, l = .2, unit = "cm"),
        axis.text.x = element_text(color = "black", size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = .5))


#ggsave(Map2, file = "FO.pdf", device = cairo_pdf, width = 13, height = 13, scale = 0.65, limitsize = FALSE, dpi = 1000)


# Map3 - British Isles ~
Map3 <-
  ggplot() + 
  geom_sf(data = GBR, fill = "#ffffff", color = "black", alpha = .5) +
  geom_sf(data = IRL, fill = "#ffffff", color = "black", alpha = .5) +
  geom_sf(data = IMN, fill = "#ffffff", color = "black", alpha = .5) +
  geom_sf(data = NR[NR$origin == "1",], fill = "#d4b9da", alpha = .35, color = NA) +
  geom_sf(data = Coords_GBR_sf, aes(fill = Class_Article), size = 5, alpha = 1, show.legend = FALSE, shape = 21, colour = "black") +
  geom_label_repel(data = Coords_GBR, size = 2.5, seed = 10, min.segment.length = 0, force = 15, segment.curvature = 1, nudge_x = 0, nudge_y = 0, max.overlaps = Inf,
                   fontface = "bold", colour = "black", aes(x = Longitude, y = Latitude, label = LocationOnly,
                   fill = Class_Article), alpha = .9, show.legend = FALSE) +
  geom_text(aes(label = "British Isles"), x = -9, y = 60.1, size = 7.5, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("#44AA99", "#F0E442"), drop = FALSE) +
  scale_colour_manual(values = c("#44AA99", "#F0E442"), drop = FALSE) +
  scale_x_continuous(breaks = seq(-12, 2, by = 2)) +
  scale_y_continuous(breaks = seq(50, 61, by = 2)) +
  coord_sf(xlim = c(-12.6, 2.3), ylim = c(49.75, 60.98), expand = FALSE) +
  annotation_north_arrow(location = "bl", which_north = "false", style = north_arrow_fancy_orienteering,
                         pad_x = unit(.05, "in"), pad_y = unit(.15, "in")) +
  annotation_scale(location = 'bl', line_width = 1.25, text_cex = 1.2, style = "ticks") +
  theme(panel.background = element_rect(fill = "#f7fbff"),
        panel.border = element_rect(colour = "#a50026", linewidth = .5, linetype = "dotdash", fill = NA),
        panel.grid.major = element_line(color = "#d9d9d9", linetype = "dashed", linewidth = 0.00005),
        plot.margin =  margin(t = 0, b = 0, r = .2, l = .2, unit = "cm"),
        axis.text.x = element_text(color = "black", size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = .5))


# Map4 - Sri Lanka ~
Map4 <-
 ggplot() + 
  geom_sf(data = SLK, fill = "#ffffff", color = "black", alpha = .5) +
  geom_sf(data = NR_slk[NR_slk$origin == "1",], fill = "#d4b9da", alpha = .35, color = NA) +
  geom_sf(data = Coords_SLK_sf, aes(fill = Class_Article), size = 5, alpha = 1, show.legend = FALSE, shape = 21, colour = "black") +
  #geom_label_repel(data = Coords_SLK, size = 4.5, seed = 10, min.segment.length = 0, force = 50, segment.curvature = 1,
  #                 nudge_x = c(0, 0, -0.35), nudge_y = c(-0.15, 0.15, 0), max.overlaps = Inf,
  #                 fontface = "bold", colour = "black", aes(x = Longitude, y = Latitude, label = LocationOnly,
  #                 fill = Class_Article), alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = "Sri Lanka"), x = 79, y = 9.6, size = 7.5, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("#44AA99", "#F0E442"), drop = FALSE) +
  scale_colour_manual(values = c("#44AA99", "#F0E442"), drop = FALSE) +
  scale_x_continuous(breaks = seq(79, 82, by = 1)) +
  scale_y_continuous(breaks = seq(6, 11, by = 1)) +
  coord_sf(xlim = c(78.4, 82), ylim = c(5.825, 9.9), expand = FALSE) +
  annotation_north_arrow(location = "bl", which_north = "false", style = north_arrow_fancy_orienteering,
                         pad_x = unit(.05, "in"), pad_y = unit(.15, "in")) +
  annotation_scale(location = "bl", line_width = 1.25, text_cex = 1.2, style = "ticks") +
  theme(panel.background = element_rect(fill = "#f7fbff"),
        panel.border = element_rect(colour = "#a50026", size = .5, linetype = "dotdash", fill = NA),
        panel.grid.major = element_line(color = "#d9d9d9", linetype = "dashed", linewidth = 0.00005),
        plot.margin = margin(t = 0, b = 0, r = .2, l = .2, "cm"),
        axis.text.x = element_text(color = "black", size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = .5))


#ggsave(Map4, file = "SLK.pdf", device = cairo_pdf, width = 13, height = 13, scale = .65, limitsize = FALSE, dpi = 600)


# Creates final panel ~
MapPanel <- Map1 / (Map2 | Map3 | Map4) + plot_layout(widths = c(1))


# Saves panel ~
ggsave(MapPanel, file = "FPG--MapWaldir.pdf", device = cairo_pdf,
       width = 20, height = 18, scale = .8, limitsize = FALSE, dpi = 600)
ggsave(MapPanel, file = "FPG--Map.jpg",
       width = 20, height = 18, scale = .8, limitsize = FALSE, dpi = 300)


#
##
### The END ~~~~~