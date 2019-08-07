# Final clean analyses for weevils 
# John Godlee (johngodlee@gmail.com)
# 2019_07_12

# Preamble ----
# Get R version
print(version)
citation(package = "base", lib.loc = NULL, auto = NULL)
citation("glmmTMB")
citation("MuMIn")
citation("emmeans")

# Remove old crap
rm(list=ls())
#dev.off()

# Set working directory to the location of the source file
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("~/google_drive/postgrad/extra_projects/weevils/analysis")

# Packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rgdal)
library(gstat)
library(igraph)
library(ggraph)
library(viridis)
library(MASS)
library(AER)
library(lme4)
library(RColorBrewer)
library(nlme)
library(MuMIn)
library(stargazer)
library(ggeffects)
library(glmmTMB)
library(emmeans)

# Import data 
damage <- read.csv("data/damage_new.csv")

site_loc <- read.csv("data/site_loc.csv")

# Spatial data
uk_outline <- readOGR("data/GBR_adm", "GBR_adm1")

europe_outline <-readOGR("data/NUTS_RG_10M_2016_4326_LEVL_0/", 
  "NUTS_RG_10M_2016_4326_LEVL_0")

cg_loc <- data.frame(long = -3.21, lat = 55.86)

geog_zone_pal <- brewer.pal(n = length(unique(site_loc$geog_zone)), name = "Set2")
geog_zone_pal[6] <- "#C42D2D"

big_region_pal <- brewer.pal(n = length(unique(site_loc$big_region)), name = "Set1")


# Merge datasets
damage_full <- left_join(site_loc, damage, by = c("site_code" = "population"))
names(damage_full) <- c("site_name", "seed_zone", "geog_zone", "site_code", 
  "big_region", "dec_latitude", "dec_longitude", 
  "site_area_ha", "gsl", "growing_deg_days_c", 
  "feb_mean_temp_c", "jul_mean_temp_c", "loc", 
  "family", "individual", "field_code", 
  "x_coord", "y_coord", "xy_coord", 
  "block", "curr_damage", "mm2_damage")

damage_full$geog_zone <- as.character(damage_full$geog_zone)
damage_full$family <- as.character(damage_full$family)
damage_full$site_name <- as.character(damage_full$site_name)
damage_full$big_region <- as.character(damage_full$big_region)

site_loc$geog_zone <-   as.character(site_loc$geog_zone)


# Create a dataset with no zero values for current damage
damage_nozero_dam <- damage_full %>%
  filter(curr_damage > 0)

# Create column on whether damage is present
damage_full$has_damage <- ifelse(damage_full$curr_damage > 0, 1, 0)

# Data description ----

# How many saplings had weevil damage?
round(
  (1 - (length(damage_full$site_name) - length(damage_nozero_dam$site_name)) / 
      length(damage_full$seed_zone)) * 100, 
  digits = 1)

length(damage_full$site_name)

length(damage_nozero_dam$site_name)

# Exploring covariance of damaged area per site

cov_summ <- damage_nozero_dam %>%
  group_by(site_name) %>%
  summarise(mean_mm2_damage = mean(mm2_damage), 
    sd_mm2_damage = sd(mm2_damage)) %>%
  mutate(cov_mm2_damage = sd_mm2_damage / mean_mm2_damage * 100) %>%
  arrange(cov_mm2_damage)

print(cov_summ, n = length(cov_summ$cov_mm2_damage))


# Which site had the highest number of affected saplings?

num_summ <- damage_nozero_dam %>%
  group_by(site_name) %>%
  summarise(n = n()) %>%
  arrange(n)

print(num_summ, n = length(num_summ$n))

area_summ <- damage_nozero_dam %>%
  group_by(site_name) %>%
  summarise(cum_mm2_damage = sum(mm2_damage)) %>%
  arrange(cum_mm2_damage)

print(area_summ, n = length(area_summ$cum_mm2_damage))

# Compare site names with their codes

name_summ <- damage_full %>%
  group_by(site_name, site_code) %>%
  summarise() %>%
  arrange(site_name)

print(name_summ, n = length(name_summ$site_name))

# Fortify shapefiles for ggplot2
scotland_outline <- uk_outline[uk_outline$NAME_1 == "Scotland",]
scotland_fort <- fortify(scotland_outline, region = "NAME_1")

europe_fort <- fortify(europe_outline, region = "NUTS_ID")

# Plot map of Europe with study region
region_map <- ggplotGrob(
  ggplot() + 
    geom_polygon(aes(x = long, y = lat, group = group), 
      colour = "black", fill = "#CCCCCC",
      data = europe_fort, alpha = 1) + 
    geom_polygon(aes(x = long, y = lat, group = group), 
      colour = NA, fill = "#6B6B6B",
      data = scotland_fort, alpha = 1) + 
    geom_rect(aes(xmin = -8, xmax = 0, ymin = 55, ymax = 60), 
      colour = "red", fill = NA) +
    theme_void() + 
    theme(legend.position = "right") + 
    coord_map() + 
    xlim(-12, 25) + 
    ylim(35, 70) + 
    xlab("Longitude") + 
    ylab("Latitude") + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
)

site_map <- ggplot() + 
  annotation_custom(grob = region_map,
    xmin = -2, xmax = +0.5,
    ymin = 55, ymax = 57.0) + 
  geom_polygon(aes(x = long, y = lat, group = group), 
    colour = "black", fill = "white",
    data = scotland_fort, alpha = 1) + 
  stat_ellipse(aes(x = dec_longitude, y = dec_latitude, 
    colour = big_region),
    geom = "polygon", fill = NA, segments = 200, 
    type = "t", level = 0.96,
    data = site_loc) + 
  geom_point(aes(x = dec_longitude, y = dec_latitude, 
    fill = geog_zone, shape = big_region),
    size = 3, colour = "black",
    data = site_loc) +
  geom_label_repel(aes(x = dec_longitude, y = dec_latitude,
    fill = geog_zone, label = site_code),
    label.padding = unit(0.13, "lines"), point.padding = 0.2, min.segment.length = 0.3,
    data = site_loc) + 
theme_classic() + 
  coord_quickmap() +  
  xlim(-7.5, 0) + 
  ylim(54.5, 59.5) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  scale_fill_manual(values = geog_zone_pal) + 
  scale_shape_manual(values = c(21, 22, 24), guide = FALSE) + 
  scale_colour_manual(values = big_region_pal) + 
  guides(
    fill = guide_legend(title="Geog. Zone", 
      override.aes = list(size = 5, label = "", shape = 15)),
    color = guide_legend(title="Region")) + 
  geom_point(aes(x = long, y = lat), 
    shape = 23, fill = "black", size = 3,
    data = cg_loc) + 
  geom_text(aes(x = long, y = lat),
    label = "Common\nGarden",
    hjust = -0.15, vjust = 0.7,
    data = cg_loc) + 
  theme(legend.position=c(.9,.75))

ggsave(file="../paper/img/site_map.pdf", plot=site_map, width=8, height=8)

# Plot map of sites with points sized by mm^2 damage found
site_summ <- damage_full %>%
  group_by(site_code) %>%
  summarise(mm2_damage = sum(mm2_damage)) %>%
  right_join(.,site_loc, by = c("site_code" = "site_code")) %>%
  mutate(geog_zone = as.character(geog_zone))

bubble_map <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), 
    colour = "black", fill = NA,
    data = scotland_fort, alpha = 1) + 
  geom_point(aes(x = dec_longitude, y = dec_latitude, 
    fill = geog_zone, size = mm2_damage),
    shape = 21, colour = "black",
    data = site_summ) + 
  geom_point(aes(x = long, y = lat), shape = 23, fill = "black",
    data = cg_loc) + 
  theme_classic() + 
  theme(legend.position=c(.9,.5)) + 
  coord_map() + 
  xlim(-7.5, 0) + 
  ylim(54.5, 59.5) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  guides(fill=guide_legend(title="Geog. Zone", override.aes = list(size = 5))) + 
  guides(size=guide_legend(title=expression(paste("Total bark\ndamage (mm"^2,")")))) + 
  scale_fill_manual(values = geog_zone_pal) + 
  geom_text(aes(x = long, y = lat),
    label = "Common\nGarden",
    hjust = -0.15, vjust = 0.7,
    data = cg_loc)

ggsave(file="../paper/img/bubble_map.pdf", plot=bubble_map, width=8, height=8)

# Plot map of sites with points sized by how many trees were attacked
site_summ_binom <- damage_full %>%
  group_by(site_code) %>%
  summarise(has_damage = sum(has_damage)) %>%
  right_join(.,site_loc, by = c("site_code" = "site_code")) %>%
  mutate(geog_zone = as.character(geog_zone))

bubble_freq_map <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), 
    colour = "black", fill = NA,
    data = scotland_fort, alpha = 1) + 
  geom_point(aes(x = dec_longitude, y = dec_latitude, 
    fill = geog_zone, size = has_damage),
    shape = 21, colour = "black",
    data = site_summ_binom) + 
  geom_point(aes(x = long, y = lat), shape = 23, fill = "black",
    data = cg_loc) + 
  theme_classic() + 
  theme(legend.position=c(.9,.5)) + 
  coord_map() + 
  xlim(-7.5, 0) + 
  ylim(54.5, 59.5) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  guides(fill=guide_legend(title="Geog. Zone", override.aes = list(size = 5))) + 
  guides(size=guide_legend(title="Saplings damaged")) + 
  scale_fill_manual(values = geog_zone_pal)

ggsave(file="../paper/img/bubble_freq_map.pdf", plot=bubble_freq_map, width=10, height=8)


# Latitudinal and longitudinal variation in damage

latitude <- ggplot(damage_nozero_dam, aes(x = dec_latitude, y = log(mm2_damage))) + 
  geom_point(aes(fill = geog_zone), 
    position = "dodge", shape = 21, colour = "black", size = 2) + 
  stat_smooth(method = "lm", colour = "black") + 
  theme_classic() +
  labs(x = "Latitude", y = expression(paste("Bark damage (mm"^2,")"))) + 
  guides(fill=guide_legend(title="Geog. Zone", override.aes = list(size = 5))) + 
  scale_fill_manual(values = geog_zone_pal)

ggsave(file="../paper/img/latitude.pdf", plot=latitude, width=10, height=5)

longitude <- ggplot(damage_full, aes(x = dec_longitude, y = mm2_damage)) + 
  geom_point(aes(colour = geog_zone)) + 
  stat_smooth(method = "lm", colour = "black") + 
  theme_classic() +
  labs(x = "Longitude", y = expression(paste("Bark damage (mm"^2,")"))) + 
  guides(colour=guide_legend(title="Geog. Zone", override.aes = list(size = 5))) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../paper/img/longitude.pdf", plot=longitude, width=10, height=8)


# Plot distribution of Growing Degree Days across seed collection sites
deg_days <- site_loc %>%
  group_by(seed_zone) %>%
  summarise(mean_growing_deg_days_c = mean(growing_deg_days_c),
    sd_growing_deg_days_c = sd(growing_deg_days_c)) %>%
  mutate(seed_zone = factor(seed_zone, 
    levels = seed_zone[order(.$mean_growing_deg_days_c)])) %>%
  ggplot(aes(x = seed_zone, y = mean_growing_deg_days_c)) + 
  geom_bar(stat = "identity", aes(fill = seed_zone)) + 
  geom_errorbar(aes(x = seed_zone, 
    ymin = (mean_growing_deg_days_c - sd_growing_deg_days_c),
    ymax = (mean_growing_deg_days_c + sd_growing_deg_days_c)),
    width = 0.5) + 
  labs(x = "Seed zone",
    y = "Mean growing degree days (C)") + 
  theme_classic() +
  theme(legend.position = "none")

ggsave(file="../paper/img/deg_days.pdf", plot=deg_days, width=10, height=8)


# Minimum and Maximum Growing Degree Days across seed collection sites
min(site_loc$growing_deg_days_c)
max(site_loc$growing_deg_days_c)

# Plot a dendrogram of nested structure of data with current damage
# Create base node
damage_full$country <- "Scotland"

# Create simple df
heirarchy <- damage_full %>% 
  dplyr::select(country, big_region, geog_zone, site_code) %>%
  distinct()

# transform it to an edge list for each level
edges_level1_2 = heirarchy %>% 
  dplyr::select(country, big_region) %>% 
  unique %>% 
  rename(from=country, to=big_region)

edges_level2_3 = heirarchy %>% 
  dplyr::select(big_region, geog_zone) %>% 
  unique %>% 
  rename(from=big_region, to=geog_zone)

edges_level3_4 = heirarchy %>% 
  dplyr::select(geog_zone, site_code) %>% 
  unique %>% 
  rename(from=geog_zone, to=site_code)

edge_list=rbind(edges_level1_2, edges_level2_3, edges_level3_4)

# Create graph object
dam_graph <- graph_from_data_frame( edge_list )

# Add line weights by damage
region_weights <- damage_full %>% 
  group_by(big_region) %>%
  summarise(damage = sum(mm2_damage / n()))

geog_zone_weights <- damage_full %>%
  group_by(geog_zone) %>%
  summarise(damage = sum(mm2_damage / n())) %>% 
  mutate(geog_zone = factor(as.character(geog_zone), levels = c("1", "2", "3", "4", "5", "6")))
geog_zone_weights <- geog_zone_weights[order(geog_zone_weights$geog_zone), ]

site_code_weights <- damage_full %>%
  group_by(site_code) %>%
  summarise(damage = sum(mm2_damage / n())) %>%
  mutate(site_code = factor(as.character(site_code), 
    levels = c("AM",
      "GE", "RD", "SO", "BE", "LC", "SD",
      "GA", "GC", "AB", "AC", "BB", "GD",
      "GT", "RM", "CG", "GL", "BW", "CC", "CR", "MG")))
site_code_weights <- site_code_weights[order(site_code_weights$site_code), ]

E(dam_graph)$weight <- c(region_weights$damage, geog_zone_weights$damage, site_code_weights$damage)

# Plot graph
dendro <- ggraph(dam_graph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal(aes(width = weight)) +
  geom_node_point() +
  geom_node_label(position = "identity", aes(label = name), label.padding = unit(0.22, "lines")) + 
  theme_void() + 
  annotate("rect", xmin = 0, xmax = 22, ymin = 1.6, ymax = 2.4, fill = "red", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = 22, ymin = 0.6, ymax = 1.4, fill = "blue", alpha = 0.1) + 
  annotate("rect", xmin = 0, xmax = 22, ymin = -0.2, ymax = 0.4, fill = "green", alpha = 0.1) +
  geom_text(aes(x = 24, y = 2), label = "Region", hjust = "center") + 
  geom_text(aes(x = 24, y = 1), label = "Geog.\nzone", hjust = "center") +
  geom_text(aes(x = 24, y = 0), label = "Site", hjust = "center") + 
  theme(legend.position = "none")

ggsave(file="../paper/img/dendro.pdf", plot=dendro, width=10, height=5)


# Bar chart number with number of saplings damaged

damage_nozero_dam_order_bar <- damage_nozero_dam %>%
  group_by(site_code, geog_zone) %>%
  summarise(n = n())

barchart <- ggplot(damage_nozero_dam_order_bar, 
  aes(x = site_code, y = n, fill = geog_zone)) + 
  geom_bar(stat = "identity", width = 1, 
    colour = "black") + 
  labs(x = "Population and Geog. Zone", 
    y = "Number of saplings damaged") + 
  scale_fill_manual(values = geog_zone_pal) + 
  facet_grid(~geog_zone, 
    scales = "free_x", 
    space = "free_x",
    switch = "x") + 
  theme_classic() + 
  theme(panel.spacing = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    strip.text = element_text(size = 12),
    legend.position = "none")

ggsave(file="../paper/img/barchart.pdf", plot=barchart, width=10, height=5)


# mm2_damage boxplot for saplings with lesions

damage_nozero_dam <- damage_nozero_dam %>%
  mutate(site_code = factor(as.character(site_code), 
    levels = c("AM",
      "GE", "RD", "SO", "BE", "LC", "SD",
      "GA", "GC", "AB", "AC", "BB", "GD",
      "GT", "RM", "CG", "GL", "BW", "CC", "CR", "MG")))

boxplot <- ggplot(damage_nozero_dam, 
  aes(x = site_code, y = log(mm2_damage), fill = geog_zone)) + 
  geom_boxplot() + 
  labs(x = "Population and Geog. Zone",
    y = expression(paste("log(Bark area damaged) (mm"^2,")"))) +
  facet_grid(~geog_zone, 
    scales = "free_x", 
    space = "free_x",
    switch = "x") + 
  theme_classic() + 
  scale_fill_manual(values = geog_zone_pal) + 
  guides(fill=guide_legend(title="Geog. Zone")) + 
  theme(legend.position = "bottom") + 
  theme(panel.spacing = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    strip.text = element_text(size = 12),
    legend.position = "none")

ggsave(file="../paper/img/boxplot.pdf", plot=boxplot, width=10, height=5)


# Plot a grid map of saplings and lesions to visualise spatial autocorrelation
damage_zero <- damage_full %>%
  filter(mm2_damage == 0)

sapling_map <- ggplot(damage_full, aes(x = x_coord, y = y_coord)) +  
  geom_point(data = damage_zero, aes(x = x_coord, y = y_coord), shape = 4, size = 1) + 
  geom_point(data = damage_nozero_dam, aes(x = x_coord, y = y_coord, colour = mm2_damage, size = mm2_damage)) + 
  scale_colour_viridis(name = expression(paste("Bark area damaged (mm"^2,")"))) + 
  scale_size_continuous(guide = FALSE) + 
  coord_equal() + 
  theme_void() + 
  theme(legend.position = "bottom") +
  ylim(-0.2, 8.8)

ggsave(file="../paper/img/sapling_map.pdf", plot=sapling_map, width=10, height=2.5)