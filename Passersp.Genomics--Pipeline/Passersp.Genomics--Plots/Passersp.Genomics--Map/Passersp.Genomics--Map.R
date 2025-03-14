### The BEGINNING ~~~~~
##
# Plots Passer sp. Genomics -- Map | Written by George Pacheco ~


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(rnaturalearth, rnaturalearthdata, ggtext, sf, ggspatial, tidyverse, ggforce, ggstar,
               ggrepel, sysfonts, extrafont, gridExtra, grid, cowplot, patchwork, ggpubr, showtext)


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
shapes.legend <- as.vector(c(9, 1, 28, 12, 11, 23))


# Creates base map ~
MyLegend_Plot <-
 ggplot() +
  geom_star(data = SamplingEffort, aes(x = Longitude, y = Latitude, starshape = Location, fill = Species),
            size = 4.85, alpha = .9, starstroke = .15) +
  scale_starshape_manual(values = shapes.legend, na.translate = FALSE) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000"), na.translate = FALSE,  drop = FALSE) +
  scale_colour_manual(values = c("#1E90FF", "#FFD700", "#ee0000"), na.translate = FALSE, drop = FALSE) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(.1, "cm")) +
  guides(fill = guide_legend(title = "Species", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                             label.theme = element_text(family = "Optima", size = 15),
                             override.aes = list(size = 5, starstroke = .15), nrow = 1, order = 1),
         starshape = guide_legend(title = "Locations", title.theme = element_text(family = "Optima", size = 16, face = "bold"),
                                  label.theme = element_text(family = "Optima", size = 15),
                                  override.aes = list(size = 5, starstroke = .15), nrow = 1, order = 2),
         colour = "none")



# Expands fulldf by adding Labels ~
SamplingEffort$Location <- ifelse(SamplingEffort$Labels %in% c("Focal Ind.", "Garderen", "Meerkerk"), "Focal Area", SamplingEffort$Location)


# Defines the shapes to be used for each Group ~
shapes.auto <- as.vector(c(9, 1, 28, 12, 11, 23, 15))


# Creates base map ~
MapBody <-
  ggplot() +
  geom_sf(data = Global, fill = "#fff7f3", color = "#bdbdbd") +
  geom_sf(data = PInr, fill = "#FFD700", alpha = .35, color = NA) +
  geom_star(data = subset(SamplingEffort, Location != "Focal Area"), aes(x = Longitude, y = Latitude, starshape = Location, fill = Species),
            size = 4.85, colour = "#000000", alpha = .9, starstroke = .15) +
  geom_star(data = subset(SamplingEffort, Location == "Focal Area"), aes(x = Longitude, y = Latitude, starshape = Location),
            size = 2.85, fill = "#d9d9d9", colour = "#000000", alpha = .85, starstroke = .15) +
  scale_starshape_manual(values = shapes.auto) +
  coord_sf(xlim = c(-13.6, 75), ylim = c(35, 61), expand = FALSE) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000")) +
  scale_colour_manual(values = c("#1E90FF", "#FFD700", "#ee0000")) +
  scale_x_continuous(breaks = seq(-120, 125, by = 10)) +
  scale_y_continuous(breaks = seq(-20, 70, by = 10)) +
  geom_label_repel(data = subset(SamplingEffort, Labels == "Garderen"), aes(x = Longitude, y = Latitude, label = Labels), show.legend = FALSE,
                   size = 3.5, fontface = "bold", colour = "#000000", fill = "#d9d9d9", alpha = .85, point.padding = .75, segment.size = .3, nudge_x = 5, nudge_y = -.5,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  geom_label_repel(data = subset(SamplingEffort, Labels == "Meerkerk"), aes(x = Longitude, y = Latitude, label = Labels), show.legend = FALSE,
                   size = 3.5, fontface = "bold", colour = "#000000", fill = "#d9d9d9", alpha = .85, point.padding = 1.5, segment.size = .3, nudge_x = -2, nudge_y = -3.1,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"), ends = "last", type = "open")) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering,
                         pad_x = unit(.125, "in"), pad_y = unit(.125, "in")) +
  annotation_scale(text_family = "Optima", location = "tr", line_width = 1, text_cex = 1, style = "ticks",
                   pad_x = unit(.2, "in"), pad_y = unit(.65, "in")) +
  theme(panel.background = element_rect(fill = "#deebf7"),
        panel.grid.major = element_line(color = "#bdbdbd", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "#000000", linewidth = .3, fill = NA),
        plot.margin = margin(t = -50, b = 0, r = 10, l = 10),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 11, color = "#000000"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


# Isolates legend ~
MyLegendBlog <- get_legend(MyLegend_Plot)


# Gets the final map ~
MapFull <- ggarrange(MapBody, nrow = 1, legend.grob = MyLegendBlog)


# Saves map ~
ggsave(MapFull, file = "Passersp.Genomics--Map.pdf", device = cairo_pdf,
       width = 14, height = 8, scale = 1, limitsize = FALSE, dpi = 600)
ggsave(MapFull, file = "Passersp.Genomics--Map.jpeg",
       width = 14, height = 8, scale = 1, limitsize = FALSE, dpi = 600)


#
##
### The END ~~~~~