# Final clean analyses for weevils 
# John Godlee (johngodlee@gmail.com)
# 2019_07_12

# Preamble ----
# Get R version
print(version)
citation(package = "base", lib.loc = NULL, auto = NULL)
citation("glmmTMB")
citation("MuMIn")

# Remove old crap
rm(list=ls())
#dev.off()

# Set working directory to the location of the source file
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("~/google_drive/postgrad/extra_projects/weevils/analysis")

# Packages
library(dplyr)
library(ggplot2)
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
library(svglite)


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
    ymin = 55.5, ymax = 57.5) + 
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
  theme_classic() + 
  theme(legend.position = "right") + 
  coord_quickmap() +  
  xlim(-7.5, 0) + 
  ylim(54.5, 59.5) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  scale_fill_manual(values = geog_zone_pal, guide = FALSE) + 
  scale_shape_manual(values = c(21, 22, 24), guide = FALSE) + 
  scale_colour_manual(values = big_region_pal, guide = FALSE) + 
  geom_point(aes(x = long, y = lat), 
    shape = 23, fill = "black", size = 3,
    data = cg_loc) + 
  geom_text(aes(x = long, y = lat),
    label = "CG",
    hjust = -0.8,
    data = cg_loc)

ggsave(file="../paper/img/site_map.pdf", plot=site_map, width=10, height=8)

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
  theme(legend.position = "right") + 
  coord_map() + 
  xlim(-7.5, 0) + 
  ylim(54.5, 59.5) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  guides(fill=guide_legend(title="Geog. zone", override.aes = list(size = 5))) + 
  guides(size=guide_legend(title=expression(paste("Total bark\ndamage (mm"^2,")")))) + 
  scale_fill_manual(values = geog_zone_pal)

ggsave(file="../paper/img/bubble_map.pdf", plot=bubble_map, width=10, height=8)

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
  theme(legend.position = "right") + 
  coord_map() + 
  xlim(-7.5, 0) + 
  ylim(54.5, 59.5) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  guides(fill=guide_legend(title="Geog. zone", override.aes = list(size = 5))) + 
  guides(size=guide_legend(title="Saplings damaged")) + 
  scale_fill_manual(values = geog_zone_pal)

ggsave(file="../paper/img/bubble_freq_map.pdf", plot=bubble_freq_map, width=10, height=8)


# Latitudinal and longitudinal variation in damage

latitude <- ggplot(damage_nozero_dam, aes(x = dec_latitude, y = log(mm2_damage))) + 
  geom_point(aes(colour = geog_zone)) + 
  stat_smooth(method = "lm", colour = "black") + 
  theme_classic() +
  labs(x = "Latitude", y = expression(paste("Bark damage (mm"^2,")"))) + 
  guides(colour=guide_legend(title="Geog. zone", override.aes = list(size = 5))) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../paper/img/latitude.pdf", plot=latitude, width=10, height=8)

longitude <- ggplot(damage_full, aes(x = dec_longitude, y = mm2_damage)) + 
  geom_point(aes(colour = geog_zone)) + 
  stat_smooth(method = "lm", colour = "black") + 
  theme_classic() +
  labs(x = "Longitude", y = expression(paste("Bark damage (mm"^2,")"))) + 
  guides(colour=guide_legend(title="Geog. zone", override.aes = list(size = 5))) + 
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

ggsave(file="../paper/img/dendro.pdf", plot=dendro, width=10, height=8)


# Bar chart number with number of saplings damaged

damage_nozero_dam_order_bar <- damage_nozero_dam %>%
  group_by(site_code, geog_zone) %>%
  summarise(n = n())

barchart <- ggplot(damage_nozero_dam_order_bar, 
  aes(x = site_code, y = n, fill = geog_zone)) + 
  geom_bar(stat = "identity", width = 1, 
    colour = "black") + 
  labs(x = "Population and Geog. zone", 
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

ggsave(file="../paper/img/barchart.pdf", plot=barchart, width=10, height=8)

  
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
  labs(x = "Population and Geog. zone",
    y = expression(paste("log(Bark area damaged) (mm"^2,")"))) +
  facet_grid(~geog_zone, 
    scales = "free_x", 
    space = "free_x",
    switch = "x") + 
  theme_classic() + 
  scale_fill_manual(values = geog_zone_pal) + 
  guides(fill=guide_legend(title="Geog. zone")) + 
  theme(legend.position = "bottom") + 
  theme(panel.spacing = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    strip.text = element_text(size = 12),
    legend.position = "none")

ggsave(file="../paper/img/boxplot.pdf", plot=boxplot, width=10, height=8)


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

ggsave(file="../paper/img/sapling_map.pdf", plot=sapling_map, width=10, height=8)


# Test for spatial auto-correlation with semi-variogram

# GLS model comparison with different auto-correlation structures
m1 <- gls(log(mm2_damage) ~ 1, correlation = corExp(form = ~x_coord + 
    y_coord, nugget = TRUE), data = damage_nozero_dam)

m2 <- gls(log(mm2_damage) ~ 1, correlation = corGaus(form = ~x_coord + 
    y_coord, nugget = TRUE), data = damage_nozero_dam)

m3 <- gls(log(mm2_damage) ~ 1, correlation = corSpher(form = ~x_coord + 
    y_coord, nugget = TRUE), data = damage_nozero_dam)

m4 <- gls(log(mm2_damage) ~ 1, correlation = corLin(form = ~x_coord + 
    y_coord, nugget = TRUE), data = damage_nozero_dam)

m5 <- gls(log(mm2_damage) ~ 1, correlation = corRatio(form = ~x_coord + 
    y_coord, nugget = TRUE), data = damage_nozero_dam)

m_null <- gls(log(mm2_damage) ~ 1, data = damage_nozero_dam)

gls_aic_comp <- data.frame(cor_struc = c("Exponential", "Gaussian", "Spherical", "Linear", "Rational quadratic", "Null"),
  AIC(m1, m2, m3, m4, m5, m_null),
  logLik = c(logLik(m1)[[1]], logLik(m2)[[1]], logLik(m3)[[1]], 
    logLik(m4)[[1]], logLik(m5)[[1]], logLik(m_null)[[1]]),
  r2 = c(r.squaredLR(m1), r.squaredLR(m2), r.squaredLR(m3), r.squaredLR(m4), r.squaredLR(m5), r.squaredLR(m_null))
  ) %>%
  mutate(model_name = rownames(.)) %>%
  dplyr::select(-df) %>%
  arrange(AIC) %>%
  dplyr::select(cor_struc, AIC, logLik, r2)

names(gls_aic_comp) <- c("Cor. Struct.", "AIC", "logLik", "r2m")

fileConn<-file("../paper/include/gls_aic_comp.tex")
writeLines(stargazer(gls_aic_comp, 
  summary = FALSE, rownames = FALSE, label = "cor_table", 
  table.placement = "H", digit.separate = 0), 
  fileConn)
close(fileConn)

# Make variogram
m2_semivar <- variogram(m2$residuals~1, data = damage_nozero_dam, 
  locations = ~x_coord + y_coord)

# Fit model
m2_semivar_fit = fit.variogram(m2_semivar, 
  model = vgm(psill = 1, model = "Exp", range = 15))

# Fortify model for ggplot2
m2_semivar_fit_fort <- variogramLine(m2_semivar_fit, 
  maxdist = max(m2_semivar$dist))

# Create semivariogram  of raw data
mm2_dam_semivar <- variogram(log(mm2_damage)~1, data = damage_nozero_dam, 
  locations = ~x_coord+y_coord)

# Fit model
mm2_dam_semivar_fit = fit.variogram(mm2_dam_semivar, 
  model = vgm(psill = 1, model = "Exp", range = 15))

# Fortify model for ggplot2
mm2_dam_semivar_fit_fort <- variogramLine(mm2_dam_semivar_fit, 
  maxdist = max(mm2_dam_semivar$dist))

# Make semivariogram plot
semivariogram <- ggplot() + 
  geom_point(data = mm2_dam_semivar, aes(x = dist, y = gamma)) + 
  geom_line(data = mm2_dam_semivar_fit_fort, aes(x = dist, y = gamma)) + 
  theme_classic() + 
  labs(x = "Distance (m)", y = "Semivariance (\u03B3)") + 
  geom_vline(aes(xintercept = 5), linetype = 2) + 
  geom_hline(aes(yintercept = max(mm2_dam_semivar_fit_fort$gamma)), linetype = 2) + 
  ylim(0.8, 1.3)

ggsave(file="../paper/img/variogram.pdf", plot=semivariogram, width=10, height=8)


##' Semivariogram indicates that there isn't 
##' any great amount of spatial autocorrelation

# Statistical modelling ----

# Mixed effects model - logistic regression
geog_rsite_rfamily_binom <- glmmTMB(has_damage ~ geog_zone + (1|site_name/family), 
  data = damage_full, family = binomial)
geog_rsite_binom <- glmmTMB(has_damage ~ geog_zone + (1|site_name), 
  data = damage_full, family = binomial)
site_rfamily_binom <- glmmTMB(has_damage ~ site_name + (1|family), 
  data = damage_full, family = binomial)
site_rgeog_binom <- glmmTMB(has_damage ~ site_name + (1|geog_zone), 
  data = damage_full, family = binomial)
site_rgeog_rfamily_binom <- glmmTMB(has_damage ~ site_name + (1|geog_zone) + (1|family), 
  data = damage_full, family = binomial)
family_rgeog_rsite_binom <- glmmTMB(has_damage ~ family + (1|geog_zone/site_name), 
  data = damage_full, family = "binomial")
family_rgeog_binom <- glmmTMB(has_damage ~ family + (1|geog_zone),
  data = damage_full, family = "binomial")

null_binom <- glmmTMB(has_damage ~ 1, 
  data = damage_full, family = "binomial")
rand_rsite_binom <- glmmTMB(has_damage ~ (1|site_name), 
  data = damage_full, family = "binomial")
rand_rsite_rfamily_binom <- glmmTMB(has_damage ~ (1|site_name/family), 
  data = damage_full, family = "binomial")
rand_rgeog_rsite_rfamily_binom <- glmmTMB(has_damage ~ (1|geog_zone/site_name/family), 
  data = damage_full, family = "binomial")

binom_list <- list(geog_rsite_rfamily_binom, geog_rsite_binom, 
  site_rfamily_binom, site_rgeog_binom,
  site_rgeog_rfamily_binom,  family_rgeog_rsite_binom, family_rgeog_binom, 
  null_binom, rand_rsite_binom, rand_rsite_rfamily_binom, rand_rgeog_rsite_rfamily_binom)

r2m_mod <- function(x){
  tryCatch(
    r.squaredGLMM(x)[1,1], 
    error = function(err) 0)
}

r2c_mod <- function(x){
  tryCatch(
    r.squaredGLMM(x)[2,1], 
    error = function(err) 0)
}

r2c_mod_lmer <- function(x){
  tryCatch(
    r.squaredGLMM(x)[2], 
    error = function(err) 0)
}

binom_lmer_aic_comp <- data.frame(
  fixed_effects = c("Geog. zone", "Geog. zone", "Site", "Site", "Site", "Family", 
    "Family", "NA", "NA", "NA", "NA"),
  random_effects = c("Site/Family", "Site", "Family", "Geog. zone", "Geog. zone + Family", 
    "Geog. zone / Site", "Geog. zone", 
    "NA", "Site", "Site / Family", "Geog. zone / Site / Family"),
  AIC = as.vector(unlist(lapply(binom_list, AIC))),
  logLik = as.vector(unlist(lapply(binom_list, logLik))),
  r2m = as.vector(unlist(lapply(binom_list, r2m_mod))),
  r2c = as.vector(unlist(lapply(binom_list, r2c_mod)))) %>%
  arrange(AIC)

binom_lmer_aic_comp

fileConn<-file("../paper/include/binom_lmer_aic_comp.tex")
writeLines(stargazer(binom_lmer_aic_comp, 
  summary = FALSE, rownames = FALSE, label = "binom_comp", digit.separate = 0), fileConn)
close(fileConn)

# Plot of predicted probabilities of `has_damage`
site_predict_binom <- ggpredict(site_rfamily_binom, c("site_name"))
site_predict_binom$x <- as.factor(site_predict_binom$x)
site_predict_binom <- data.frame(site_predict_binom)
site_predict_binom$x <- sort(unique(damage_full$site_code))
site_predict_binom$geog_zone <- damage_full %>% 
  group_by(site_code) %>% 
  summarise(geog_zone = first(geog_zone)) %>% 
  pull(geog_zone)


pred_binom_site <- ggplot() + 
  geom_point(data = site_predict_binom, 
    aes(x = x, y = predicted, colour = geog_zone),
    size = 3) + 
  geom_errorbar(data = site_predict_binom, 
    aes(x = x, ymin = conf.low, ymax = conf.high, colour = geog_zone),
    width = 0.5) + 
  facet_grid(~geog_zone, 
    scales = "free_x", 
    space = "free_x",
    switch = "x") + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(x = "Site", y = "Predicted probability of damage") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
    breaks = c(seq(from = 0, to = 1, by = 0.1)), 
    limits = c(0,0.73)) + 
  theme(panel.spacing = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    strip.text = element_text(size = 12),
    legend.position = "none") + 
  scale_fill_manual(values = geog_zone_pal) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../paper/img/pred_binom_site.pdf", plot=pred_binom_site, width=10, height=8)


geog_zone_predict_binom <- ggpredict(geog_rsite_binom, c("geog_zone"))
geog_zone_predict_binom$x <- as.character(geog_zone_predict_binom$x)

pred_binom_geog <- ggplot() + 
  geom_point(data = geog_zone_predict_binom, 
    aes(x = x, y = predicted, colour = x),
    size = 3) + 
  geom_errorbar(data = geog_zone_predict_binom, 
    aes(x = x, ymin = conf.low, ymax = conf.high, colour = x),
    width = 0.5) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(x = "Geog. zone", y = "Predicted probability of damage") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
    breaks = c(seq(from = 0, to = 0.6, by = 0.1)), 
    limits = c(0,0.6)) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../paper/img/pred_binom_geog.pdf", plot=pred_binom_geog, width=10, height=8)


# Tukeys's HSD test to see which groups are different.

# Geog. zone
tukey_geog_zone_no_family_lmer <- emmeans(geog_rsite_binom, "geog_zone")

tukey_lmer_geog <- pwpp(tukey_geog_zone_no_family_lmer, values = TRUE) + 
  theme_classic() + 
  scale_colour_manual(values = geog_zone_pal) + 
  ylab("Geog. zone") + 
  theme(panel.grid.major.y = element_line(colour="#E0E0E0"))

ggsave(file="../paper/img/tukey_lmer_geog.pdf", plot=tukey_lmer_geog, width=10, height=8)


# Site
tukey_site_lmer <- emmeans(site_rfamily_binom, "site_name")

tukey_lmer_site <- pwpp(tukey_site_lmer, values = TRUE) + 
  theme_classic() +
  ylab("Geog. zone") + 
  theme(panel.grid.major.y = element_line(colour="#E0E0E0")) + 
  scale_x_continuous(breaks = c(0.05, 0.1, 0.5, 1), limits = c(0,1))

  
  
ggsave(file="../paper/img/tukey_lmer_site.pdf", plot=tukey_lmer_site, width=10, height=8)

# Mixed effects model - Non-zero glm

geog_rsite_rfamily_lmer <- glmmTMB(log(mm2_damage) ~ geog_zone + (1|site_name/family), 
  data = damage_nozero_dam)
geog_rsite_lmer <- glmmTMB(log(mm2_damage) ~ geog_zone + (1|site_name), 
  data = damage_nozero_dam)
site_rfamily_lmer <- glmmTMB(log(mm2_damage) ~ site_name + (1|family), 
  data = damage_nozero_dam)
site_rgeog_lmer <- glmmTMB(log(mm2_damage) ~ site_name + (1|geog_zone), 
  data = damage_nozero_dam)
site_rgeog_rfamily_lmer <- glmmTMB(log(mm2_damage) ~ site_name + (1|geog_zone) + (1|family), 
  data = damage_nozero_dam)
family_rgeog_rsite_lmer <- glmmTMB(log(mm2_damage) ~ family + (1|geog_zone/site_name),
  data = damage_nozero_dam)
family_rgeog_lmer <- glmmTMB(log(mm2_damage) ~ family + (1|geog_zone),
  data = damage_nozero_dam)
null_lmer <- glmmTMB(log(mm2_damage) ~ 1, 
  data = damage_nozero_dam)
rand_rsite_rfamily_lmer <- glmmTMB(log(mm2_damage) ~ (1|site_name/family), 
  data = damage_nozero_dam)
rand_rgeog_rsite_lmer <- glmmTMB(log(mm2_damage) ~ (1|geog_zone/site_name), 
  data = damage_nozero_dam)

lmer_list <- list(
  geog_rsite_rfamily_lmer, geog_rsite_lmer, site_rfamily_lmer, 
  site_rgeog_lmer, site_rgeog_rfamily_lmer, family_rgeog_rsite_lmer, 
  family_rgeog_lmer, null_lmer, rand_rsite_rfamily_lmer, rand_rgeog_rsite_lmer)

lmer_aic_comp <- data.frame(
  fixed_effects = c("Geog. zone", "Geog. zone", "Site", "Site", "Site",
    "Family", "Family", "NA", "NA", "NA"),
  random_effects = c("Site / Family", "Site", "Family", "Geog. zone", 
    "Geog. zone + Family", "Geog. zone / Site", "Geog. zone", 
    "NA", "Site / Family", "Geog. zone / Site"),
  AIC = as.vector(unlist(lapply(lmer_list, AIC))),
  logLik = as.vector(unlist(lapply(lmer_list, logLik))),
  r2m = as.vector(unlist(lapply(lmer_list, r2m_mod))),
  r2c = as.vector(unlist(lapply(lmer_list, r2c_mod_lmer)))) %>%
  arrange(AIC)

fileConn<-file("../paper/include/lmer_aic_comp.tex")
writeLines(stargazer(lmer_aic_comp, 
  summary = FALSE, rownames = FALSE,
  label = "lmer_comp", digit.separate = 0), fileConn)
close(fileConn)

# Plots of predicted values
site_predict_lmer <- ggpredict(site_rfamily_lmer, c("site_name"))
site_predict_lmer$x <- as.factor(site_predict_lmer$x)
site_predict_lmer <- data.frame(site_predict_lmer)
site_predict_lmer$x <- sort(unique(damage_full$site_code))
site_predict_lmer$geog_zone <- damage_full %>% 
  group_by(site_code) %>% 
  summarise(geog_zone = first(geog_zone)) %>% 
  pull(geog_zone)

pred_lmer_site <- ggplot() + 
  geom_point(data = site_predict_lmer, 
    aes(x = x, y = predicted, colour = geog_zone),
    size = 3) + 
  geom_errorbar(data = site_predict_lmer, 
    aes(x = x, ymin = conf.low, ymax = conf.high, colour = geog_zone),
    width = 0.5) + 
  facet_grid(~geog_zone, 
    scales = "free_x", 
    space = "free_x",
    switch = "x") + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(x = "Geog. zone and Seed population", y = expression(paste("Predicted bark damage (mm"^2,")"))) + 
  theme(panel.spacing = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    strip.text = element_text(size = 12),
    legend.position = "none") + 
  scale_fill_manual(values = geog_zone_pal) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../paper/img/pred_lmer_site.pdf", plot=pred_lmer_site, width=10, height=8)


geog_zone_predict_lmer <- ggpredict(geog_rsite_lmer, c("geog_zone"))
geog_zone_predict_lmer$x <- as.character(geog_zone_predict_lmer$x)


pred_lmer_geog <- ggplot() + 
  geom_point(data = geog_zone_predict_lmer, 
    aes(x = x, y = predicted, colour = x),
    size = 3) + 
  geom_errorbar(data = geog_zone_predict_lmer, 
    aes(x = x, ymin = conf.low, ymax = conf.high, colour = x),
    width = 0.5) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(x = "Geog. zone", y = expression(paste("Predicted bark damage (mm"^2,")"))) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../paper/img/pred_lmer_geog.pdf", plot=pred_lmer_geog, width=10, height=8)

# Tukeys's HSD test to see which groups are different.

# Family
tukey_geog_zone_no_family_lmer <- emmeans(geog_rsite_lmer, "geog_zone")

tukey_lmer_geog <- pwpp(tukey_geog_zone_no_family_lmer, values = TRUE) + 
  theme_classic() + 
  scale_colour_manual(values = geog_zone_pal) + 
  ylab("Geog. zone") + 
  theme(panel.grid.major.y = element_line(colour="#E0E0E0"))

ggsave(file="../paper/img/tukey_lmer_geog.pdf", plot=tukey_lmer_geog, width=10, height=8)


# Site
tukey_site_lmer <- emmeans(site_rfamily_lmer, "site_name")

tukey_lmer_site <- pwpp(tukey_site_lmer, values = TRUE) + 
  theme_classic() +
  ylab("Geog. zone") + 
  theme(panel.grid.major.y = element_line(colour="#E0E0E0"))

ggsave(file="../paper/img/tukey_lmer_site.pdf", plot=tukey_lmer_site, width=10, height=8)

# Linear mixed effects model of latitude and longitude against total damage

dec_lat_lmer <- glmmTMB(log(mm2_damage) ~ dec_latitude + (1|family), data = damage_nozero_dam)
dec_lon_lmer <- glmmTMB(log(mm2_damage) ~ dec_longitude + (1|site_name/family), data = damage_nozero_dam)

summary(dec_lat_lmer)
summary(dec_lon_lmer)

test <- Anova(dec_lat_lmer)

r.squaredGLMM(dec_lat_lmer)
r.squaredGLMM(dec_lon_lmer)


dec_lat_predict <- ggpredict(dec_lat_lmer, terms = c("dec_latitude"), type = "fe")

# Plot predicted values by fixed effect
pred_lat <- ggplot() + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
    data = dec_lat_predict, alpha = 0.3) + 
  geom_line(data = dec_lat_predict, 
    aes(x = x, y = predicted)) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  labs(x = "Latitude", y = expression(paste("Predicted bark damage (mm"^2,")"))) 

ggsave(file="../paper/img/pred_lat.pdf", plot=pred_lat, width=10, height=8)

