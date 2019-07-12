# Final clean analyses for weevils 
# John Godlee (johngodlee@gmail.com)
# 2019_07_12

# Preamble ----
# Get R version
print(version)
citation(package = "base", lib.loc = NULL, auto = NULL)
citation("nlme")

# Remove old crap
rm(list=ls())
dev.off()

# Set working directory to the location of the source file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library(dplyr)
library(ggplot2)
library(rgdal)
library(gstat)

# Import data 
damage <- read.csv("data/damage.csv")

site_loc <- read.csv("data/site_loc.csv")

# Spatial data
uk_outline <- readOGR("data/GBR_adm", "GBR_adm1")

europe_outline <-readOGR("data/NUTS_RG_10M_2016_4326_LEVL_0/", 
  "NUTS_RG_10M_2016_4326_LEVL_0")

cg_loc <- data.frame(long = -3.21, lat = 55.86)

# Merge datasets
damage_full <- left_join(site_loc, damage, by = c("site_code" = "population"))
names(damage_full) <- c("site_name", "seed_zone", "site_code", 
  "big_region", "dec_latitude", "dec_longitude", 
  "site_area_ha", "gsl", "growing_deg_days_c", 
  "feb_mean_temp_c", "jul_mean_temp_c", "loc", 
  "family", "individual", "field_code", 
  "x_coord", "y_coord", "xy_coord", 
  "block", "prev_damage", "curr_damage", "total_damage")

# Create a dataset with no zero values for current damage
damage_nozero_dam <- damage_full %>%
  filter(prev_damage > 0,
    curr_damage > 0)

# Create column on whether damage is present
damage_full$has_damage <- ifelse(damage_full$curr_damage > 0, 1, 0)

# Data description ----

# How many saplings had weevil damage?
round(
  (1 - (length(damage_full$site_name) - length(damage_nozero_dam$site_name)) / 
      length(damage_full$seed_zone)) * 100, 
  digits = 1)

# Fortify shapefiles for ggplot2
scotland_outline <- uk_outline[uk_outline$NAME_1 == "Scotland",]
scotland_fort <- fortify(scotland_outline, region = "NAME_1")

europe_fort <- fortify(europe_outline, region = "NUTS_ID")

# Plot map of scotland with sites
ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), 
    colour = "black", fill = NA,
    data = scotland_fort, alpha = 1) + 
  geom_point(aes(x = dec_longitude, y = dec_latitude, 
    colour = seed_zone, shape = big_region),
    data = site_loc) + 
  geom_point(aes(x = long, y = lat), shape = 23, fill = "black",
    data = cg_loc) + 
  theme_classic() + 
  theme(legend.position = "right") + 
  coord_map() + 
  xlim(-7.5, 0) + 
  ylim(54.5, 59.5) + 
  xlab("Longitude") + 
  ylab("Latitude")

# Plot map of Europe with study region
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
  ylab("Latitude")


# Plot map of sites with points sized by how many lesions found
site_summ <- damage_full %>%
  group_by(site_code) %>%
  summarise(curr_damage = sum(curr_damage)) %>%
  right_join(.,site_loc, by = c("site_code" = "site_code"))

ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), 
    colour = "black", fill = NA,
    data = scotland_fort, alpha = 1) + 
  geom_point(aes(x = dec_longitude, y = dec_latitude, 
    fill = seed_zone, size = curr_damage),
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
  ylab("Latitude")

# Plot distribution of Growing Degree Days across seed collection sites
site_loc %>%
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

# Minimum and Maximum Growing Degree Days across seed collection sites
min(site_loc$growing_deg_days_c)
max(site_loc$growing_deg_days_c)

# Plot a dendrogram of nested structure of data with current damage
# Create base node
damage_full$country <- "scotland"

# Create simple df
heirarchy <- damage_full %>% 
  dplyr::select(country, big_region, seed_zone, site_code) %>%
  distinct()

# transform it to an edge list for each level
edges_level1_2 = heirarchy %>% 
  dplyr::select(country, big_region) %>% 
  unique %>% 
  rename(from=country, to=big_region)

edges_level2_3 = heirarchy %>% 
  dplyr::select(big_region, seed_zone) %>% 
  unique %>% 
  rename(from=big_region, to=seed_zone)

edges_level3_4 = heirarchy %>% 
  dplyr::select(seed_zone, site_code) %>% 
  unique %>% 
  rename(from=seed_zone, to=site_code)

edge_list=rbind(edges_level1_2, edges_level2_3, edges_level3_4)

# Create graph object
dam_graph <- graph_from_data_frame( edge_list )

# Add line weights by damage
region_weights <- damage_full %>% 
  group_by(big_region) %>%
  summarise(damage = sum(curr_damage))

seed_zone_weights <- damage_full %>%
  group_by(seed_zone) %>%
  summarise(damage = sum(curr_damage)) %>% 
  mutate(seed_zone = factor(as.character(seed_zone), levels = c("N", "NW", "NC", "EC", "NE", "SW", "SC")))
seed_zone_weights <- seed_zone_weights[order(seed_zone_weights$seed_zone), ]

site_code_weights <- damage_full %>%
  group_by(site_code) %>%
  summarise(damage = sum(curr_damage)) %>%
  mutate(site_code = factor(as.character(site_code), levels = c( "GE", "RD", "SO",
    "BE", "LC", "SD", "AM", "GA", "GC", "AB", "GD", "RM", "AC",
    "BB", "GT", "CG", "CR", "GL", "BW", "CC", "MG")))
site_code_weights <- site_code_weights[order(site_code_weights$site_code), ]

E(dam_graph)$weight <- c(region_weights$damage, seed_zone_weights$damage, site_code_weights$damage)

# Plot graph
ggraph(dam_graph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal(aes(width = weight)) +
  geom_node_point() +
  geom_node_label(position = "identity", aes(label = name)) + 
  theme_void() + 
  annotate("rect", xmin = 0, xmax = 22, ymin = 1.6, ymax = 2.4, fill = "red", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = 22, ymin = 0.6, ymax = 1.4, fill = "blue", alpha = 0.1) + 
  annotate("rect", xmin = 0, xmax = 22, ymin = -0.2, ymax = 0.4, fill = "green", alpha = 0.1) +
  geom_text(aes(x = 24, y = 2), label = "Region", hjust = "right") + 
  geom_text(aes(x = 24, y = 1), label = "Seed\nzone", hjust = "right") +
  geom_text(aes(x = 24, y = 0), label = "Site", hjust = "right")


# Bar chart number with lesions
ggplot(damage_nozero_dam, aes(x = seed_zone)) + 
  geom_bar(colour = "black" , aes(group = site_code, fill = seed_zone), position = "dodge") + 
  theme_classic() + 
  labs(x = "Seed zone", 
    y = "Number of saplings with lesions") + 
  theme(legend.position = "none")

# Number of lesions boxplot for saplings with lesions
ggplot(damage_nozero_dam, aes(x = site_code, y = log(curr_damage))) + 
  geom_boxplot(aes(fill = seed_zone)) + 
  labs(x = "Population",
    y = "log(Total number of lesions)") +
  theme_classic()

# Plot a grid map of saplings and lesions to visualise spatial autocorrelation
ggplot(damage_full, aes(x = y_coord, y = x_coord)) + 
  geom_point(aes(colour = total_damage, size = curr_damage)) + 
  scale_colour_viridis(name = "No. lesions") + 
  scale_size_continuous(guide = FALSE) + 
  coord_equal() + 
  theme_classic() + 
  theme(legend.position = "bottom")

# Test for spatial auto-correlation with semi-variogram

# Create semivariogram 
curr_dam_semivar <- variogram(curr_damage~1, data = damage_full, locations = ~x_coord+y_coord)

# Fit model
curr_dam_semivar_fit = fit.variogram(curr_dam_semivar, 
  model = vgm(psill = 1, model = "Exp", range = 0.5))

# Fortify model for ggplot2
curr_dam_semivar_fit_fort <- variogramLine(curr_dam_semivar_fit, maxdist = max(curr_dam_semivar$dist))

# Make semivariogram plot
ggplot() + 
  geom_point(data = curr_dam_semivar, aes(x = dist, y = gamma)) + 
  geom_line(data = curr_dam_semivar_fit_fort, aes(x = dist, y = gamma)) + 
  ylim(0, 52) + 
  theme_classic() + 
  labs(x = "Distance (m)", y = "Semivariance (\u03B3)") + 
  geom_vline(aes(xintercept = 1), linetype = 2) + 
  geom_hline(aes(yintercept = max(curr_dam_semivar_fit_fort$gamma)), linetype = 2)

##' Semivariogram indicates that there isn't 
##' any great amount of spatial autocorrelation

# Statistical modelling ----

# Logistic Regression
glm_binom_site <- glm(has_damage ~ site_name, data = damage_full, family = "binomial")
glm_binom_seed_zone <- glm(has_damage ~ seed_zone, data = damage_full, family = "binomial")
glm_binom_big_region <- glm(has_damage ~ big_region, data = damage_full, family = "binomial")
glm_binom_null <- glm(has_damage ~ 1, data = damage_full, family = "binomial")

AIC(glm_binom_site, glm_binom_seed_zone, glm_binom_big_region, glm_binom_null)


##' None of these models are better than a null model at explaining 
##' probability of having damage

# Test for overdispersion on none zero data
glm_pois <- glm(curr_damage ~ 1, data = damage_nozero_dam, family = "poisson")
dispersiontest(glm_pois, trafo=1)

##' Yes, Overdispersion present, very high alpha

# Construct Negative Binomial and poisson models
glm_negbi_site <- glm.nb(curr_damage ~ site_name, data = damage_nozero_dam, )
glm_negbi_seed_zone <- glm.nb(curr_damage ~ seed_zone, data = damage_nozero_dam)
glm_negbi_big_region <- glm.nb(curr_damage ~ big_region, data = damage_nozero_dam)
glm_negbi_null <- glm.nb(curr_damage ~ 1, data = damage_nozero_dam)

glm_pois_site <- glm(curr_damage ~ site_name, data = damage_nozero_dam, family = "poisson")
glm_pois_seed_zone <- glm(curr_damage ~ seed_zone, data = damage_nozero_dam, family = "poisson")
glm_pois_big_region <- glm(curr_damage ~ big_region, data = damage_nozero_dam, family = "poisson")
glm_pois_null <- glm(curr_damage ~ 1, data = damage_nozero_dam, family = "poisson")

# Create list of models
mod_list <- list(glm_negbi_site, glm_negbi_seed_zone, glm_negbi_big_region, glm_negbi_null,
  glm_pois_site, glm_pois_seed_zone, glm_pois_big_region, glm_pois_null)

# Log likelihood model estimates
loglik_comp <- data.frame(model_name = aic_comp$model_name, loglik = unlist(lapply(mod_list, logLik)))

# AIC values for models
aic_comp <- AIC(glm_negbi_site, glm_negbi_seed_zone, glm_negbi_big_region, glm_negbi_null,
  glm_pois_site, glm_pois_seed_zone, glm_pois_big_region, glm_pois_null) %>%
  mutate(model_name = rownames(.)) %>%
  dplyr::select(-df) %>%
  arrange(AIC) %>%
  dplyr::select(model_name, AIC)

# ANOVA of Neg-binomial models
anova(glm_negbi_null, glm_negbi_seed_zone, glm_negbi_big_region, glm_negbi_site,
  test="Chisq") 

# ANOVA of poisson models
anova(glm_pois_null, glm_pois_seed_zone, glm_pois_big_region, glm_pois_site, test = "Chisq")
