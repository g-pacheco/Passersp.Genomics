### The BEGINNING ~~~~~
##
# NLSparrow--Map | By George Pacheco ~


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(rnaturalearth, rnaturalearthdata, rgeos, sf, ggspatial, tidyverse,
               ggrepel, extrafont, cowplot, gridExtra, patchwork, jpeg, png, grid, magick)


# Imports .shp files ~
Global <- ne_countries(scale = 'large', returnclass = 'sf')
PDnr <- read_sf(dsn = ".", layer = "Passer_domesticus")
PInr <- read_sf(dsn = ".", layer = "Passer_italiae")
PHnr <- read_sf(dsn = ".", layer = "Passer_hispaniolensis")
PMnr <- read_sf(dsn = ".", layer = "Passer_montanus")


# Creates base map ~
Map <-
 ggplot() +
  #geom_sf(data = PDnr, fill = "#1E90FF", alpha = .9, color = NA) +
  geom_sf(data = PMnr, fill = "#fed976", alpha = .9, color = NA) +
  #geom_sf(data = PInr, fill = "#FFD700", alpha = .9, color = NA) +
  #geom_sf(data = PHnr, fill = "#ee0000", alpha = .9, color = NA) +
  geom_sf(data = Global, fill = "#ffffff", alpha = .25, color = "#000000") +
  #geom_sf(data = Coords_Global_sf, aes(fill = Class_Article), size = 3, alpha = 1, show.legend = "point", shape = 21, colour = "black") +
  #coord_sf(xlim = c(-13.6, 30), ylim = c(35, 61), expand = FALSE) +
  #geom_label_repel(data = Coords_Global, size = 3.75, seed = 10, min.segment.length = 0, force = 25, segment.curvature = 1,
  #                 nudge_x = 0, nudge_y = 0, max.overlaps = Inf, fontface = "bold", colour = "black",
  #                 aes(x = Longitude, y = Latitude, label = LocationOnly,
  #                 fill = Class_Article), alpha = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = c("#44AA99", "#F0E442", "#E69F00"), drop = FALSE) +
  scale_colour_manual(values = c("#44AA99", "#F0E442", "#E69F00"), drop = FALSE) +
  scale_x_continuous(breaks = seq(-120, 125, by = 10)) +
  scale_y_continuous(breaks = seq(-20, 70, by = 10)) +
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_fancy_orienteering,
                         height = unit(2, "cm"), width = unit(2, "cm"), pad_y = unit(3.25, "cm")) +
  annotation_scale(location = 'bl', line_width = 1.25, text_cex = 1.2, style = "ticks", pad_y = unit(3, "cm")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_rect(colour = "black", linewidth = .5, fill = NA),
        panel.grid.major = element_line(color = "#d9d9d9", linetype = "dashed", linewidth = .00005),
        plot.margin = margin(t = .005, b = .005, r = .2, l = .2, unit = "cm"),
        legend.background = element_rect(fill = "#c6dbef", linewidth = .15, color = "#5e5e5e", linetype = "dotted"),
        legend.key = element_rect(fill = "#c6dbef"),
        legend.position = "top",
        axis.text.x = element_text(color = "black", size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black", size = 13),
        axis.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = .5)) +
  guides(fill = guide_legend(title = "Biological Status", title.theme = element_text(size = 10.5, face = "bold"),
                             label.theme = element_text(size = 8, face = "italic"),
                             override.aes = list(size = 3, alpha = .9)), colour = "none")


#ggsave(Map, file = "XXX.pdf", device = cairo_pdf,
#       width = 16, height = 16, scale = 1, limitsize = FALSE, dpi = 600)
ggsave(Map, file = "YYY.jpeg",
       width = 16, height = 16, scale = 1, limitsize = FALSE, dpi = 600)








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