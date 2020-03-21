# Statistical modelling final version - updated from Kyle's comments
# John Godlee (johngodlee@gmail.com)
# 2019_10_15

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
#setwd("~/git_proj/weevils/analysis")

# Packages 
library(dplyr)
library(ggplot2)
library(ggeffects)
library(ggrepel)
library(viridis)
library(glmmTMB)
library(RColorBrewer)
library(stargazer)
library(emmeans)
library(MuMIn)
library(car)
library(nlme)
library(automap)
library(gridExtra)
library(gstat)
library(tidyr)
library(vegan)

# Import data ----
damage <- read.csv("data/damage_new.csv")
site_loc <- read.csv("data/site_loc.csv")
chem <- read.csv("data/chem.csv")

geog_zone_pal <- brewer.pal(n = length(unique(site_loc$geog_zone)), name = "Set2")
geog_zone_pal[6] <- "#C42D2D"

big_region_pal <- brewer.pal(n = length(unique(site_loc$big_region)), name = "Set1")

# Clean data ----
# Merge datasets
damage_full <- left_join(site_loc, damage, by = c("site_code" = "population"))
names(damage_full) <- c("site_name", "seed_zone", "geog_zone", "site_code", 
  "big_region", "dec_latitude", "dec_longitude", 
  "site_area_ha", "gsl", "growing_deg_days_c", 
  "feb_mean_temp_c", "jul_mean_temp_c", "loc", 
  "family", "individual", "field_code", 
  "x_coord", "y_coord", "xy_coord", 
  "block", "curr_damage", "mm2_damage")
damage_full <- left_join(damage_full, chem, by = c("field_code" = "tree"))

damage_full$geog_zone <- as.character(damage_full$geog_zone)
damage_full$family <- as.character(damage_full$family)
damage_full$site_name <- as.character(damage_full$site_name)
damage_full$big_region <- as.character(damage_full$big_region)

site_loc$geog_zone <- as.character(site_loc$geog_zone)

# damage_full$mm2_damage <- damage_full$mm2_damage / damage_full$dry_mass_per

# Create a dataset with no zero values for current damage
damage_nozero_dam <- damage_full %>%
  filter(curr_damage > 0)

# Create column on whether damage is present
damage_full$has_damage <- ifelse(damage_full$curr_damage > 0, 1, 0)

# Create spatial dataset for when I'm plotting semivariograms
damage_full_sp <- damage_full
coordinates(damage_full_sp) = ~x_coord + y_coord

damage_nozero_dam_sp <- damage_nozero_dam
coordinates(damage_nozero_dam_sp) = ~x_coord + y_coord

# Mixed effects model - logistic regression ----
geog_rfamily_binom <- glmmTMB(has_damage ~ geog_zone + (1|family),
  data = damage_full, family = binomial)
geog_rsite_rfamily_binom <- glmmTMB(has_damage ~ geog_zone + (1|site_name/family), 
  data = damage_full, family = binomial)
site_rfamily_binom <- glmmTMB(has_damage ~ site_code + (1|family), 
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
rand_rfamily_binom <- glmmTMB(has_damage ~ (1|family), 
  data = damage_full, family = "binomial")
rand_rsite_binom <- glmmTMB(has_damage ~ (1|site_name), 
  data = damage_full, family = "binomial")
rand_rsite_rfamily_binom <- glmmTMB(has_damage ~ (1|site_name/family), 
  data = damage_full, family = "binomial")
rand_rgeog_rsite_rfamily_binom <- glmmTMB(has_damage ~ (1|geog_zone/site_name/family), 
  data = damage_full, family = "binomial")

binom_list <- list(geog_rfamily_binom, geog_rsite_rfamily_binom, 
  site_rfamily_binom, site_rgeog_binom,
  site_rgeog_rfamily_binom, null_binom, 
  rand_rfamily_binom, rand_rsite_binom, 
  rand_rsite_rfamily_binom, rand_rgeog_rsite_rfamily_binom)

r2c_mod <- function(x){
  tryCatch(
    r.squaredGLMM(x)[1,1], 
    error = function(err) 0)
}

r2m_mod <- function(x){
  tryCatch(
    r.squaredGLMM(x)[2,1], 
    error = function(err) 0)
}

r2m_mod_lmer <- function(x){
  tryCatch(
    r.squaredGLMM(x)[2], 
    error = function(err) 0)
}

binom_lmer_aic_comp <- data.frame(
  fixed_effects = c("Geog. Zone", "Geog. Zone", "Site", "Site", "Site",
    "NA", "NA", "NA", "NA", "NA"),
  random_effects = c("Parent", "Site / Parent", "Parent", "Geog. Zone", "Geog. Zone + Parent", 
    "NA", "Parent", "Site", "Site / Parent", "Geog. Zone / Site / Parent"),
  AIC = as.vector(unlist(lapply(binom_list, AIC))),
  logLik = as.vector(unlist(lapply(binom_list, logLik))),
  r2m = as.vector(unlist(lapply(binom_list, r2m_mod))),
  r2c = as.vector(unlist(lapply(binom_list, r2c_mod)))) %>%
  arrange(AIC)

binom_lmer_aic_comp

fileConn<-file("../manuscript/include/binom_lmer_aic_comp.tex")
writeLines(stargazer(binom_lmer_aic_comp, 
  summary = FALSE, rownames = FALSE, label = "binom_comp", digit.separate = 0), fileConn)
close(fileConn)

# Plot of predicted probabilities of `has_damage`

## Best Geog. Zone model
geog_zone_predict_binom <- ggpredict(geog_rfamily_binom, c("geog_zone"))
geog_zone_predict_binom$x <- as.character(geog_zone_predict_binom$x)

pred_binom_geog <- ggplot() + 
  geom_point(data = geog_zone_predict_binom, 
    aes(x = x, y = predicted, colour = x),
    size = 3) + 
  geom_errorbar(data = geog_zone_predict_binom, 
    aes(x = x, ymin = conf.low, ymax = conf.high, colour = x),
    width = 0.5) + 
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)) + 
  labs(x = "Geog. Zone", y = "Predicted probability of damage") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
    breaks = c(seq(from = 0, to = 0.6, by = 0.1)), 
    limits = c(0,0.6)) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../manuscript/img/pred_binom_geog.pdf", plot=pred_binom_geog, width=10, height=8)

## Is there spatial autocorrelation amongst residuals in the best model?
damage_full$geog_rfamily_binom_resid <- resid(geog_rfamily_binom)

geog_rfamily_binom_autofit = autofitVariogram(damage_full$geog_rfamily_binom_resid~1, 
  damage_full_sp, model = "Exp")

geog_rfamily_binom_vario <- variogram(damage_full$geog_rfamily_binom_resid~1, 
  locations = ~x_coord + y_coord, data = damage_full)

geog_rfamily_binom_vgm <- fit.variogram(geog_rfamily_binom_vario, 
  model = vgm(model = "Exp", 
    psill = geog_rfamily_binom_autofit$var_model$psill[2], 
    range = geog_rfamily_binom_autofit$var_model$range[2],
    kappa = geog_rfamily_binom_autofit$var_model$kappa[2],
    nugget = geog_rfamily_binom_autofit$var_model$psill[1]))

geog_rfamily_binom_vgm_fort <- variogramLine(geog_rfamily_binom_vgm, 
  maxdist = max(geog_rfamily_binom_vario$dist))

# Make semivariogram plot
semivariogram_geog_zone_binom <- ggplot() + 
  geom_point(data = geog_rfamily_binom_vario, 
    aes(x = dist, y = gamma)) + 
  geom_line(data = geog_rfamily_binom_vgm_fort, 
    aes(x = dist, y = gamma)) + 
  theme_classic() + 
  labs(x = "Distance (m)", y = "Semivariance (\u03B3)") + 
  geom_vline(aes(xintercept = geog_rfamily_binom_autofit$var_model$range[2]), linetype = 2) + 
  ylim(0.15, 0.26)

ggsave(file="../manuscript/img/semivariogram_geog_zone_binom.pdf", 
  plot=semivariogram_geog_zone_binom, width=10, height=5)

## Best site model 
site_predict_binom <- ggpredict(site_rfamily_binom, c("site_code"))
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
    legend.position = "none", 
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)) + 
  scale_fill_manual(values = geog_zone_pal) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../manuscript/img/pred_binom_site.pdf", plot=pred_binom_site, width=10, height=8)

## Is there spatial autocorrelation amongst residuals in the best model?
damage_full$site_rfamily_binom_resid <- resid(site_rfamily_binom)

site_rfamily_binom_autofit = autofitVariogram(damage_full$site_rfamily_binom_resid~1, 
  damage_full_sp, model = "Exp")

site_rfamily_binom_vario <- variogram(damage_full$site_rfamily_binom_resid~1, 
  locations = ~x_coord + y_coord, data = damage_full)

site_rfamily_binom_vgm <- fit.variogram(site_rfamily_binom_vario, 
  model = vgm(model = "Exp", 
    psill = site_rfamily_binom_autofit$var_model$psill[2], 
    range = site_rfamily_binom_autofit$var_model$range[2],
    kappa = site_rfamily_binom_autofit$var_model$kappa[2],
    nugget = site_rfamily_binom_autofit$var_model$psill[1]))

site_rfamily_binom_vgm_fort <- variogramLine(site_rfamily_binom_vgm, 
  maxdist = max(site_rfamily_binom_vario$dist))

# Make semivariogram plot
semivariogram_site_binom <- ggplot() + 
  geom_point(data = site_rfamily_binom_vario, 
    aes(x = dist, y = gamma)) + 
  geom_line(data = site_rfamily_binom_vgm_fort, 
    aes(x = dist, y = gamma)) + 
  theme_classic() + 
  labs(x = "Distance (m)", y = "Semivariance (\u03B3)") + 
  geom_vline(aes(xintercept = site_rfamily_binom_autofit$var_model$range[2]), linetype = 2) + 
  ylim(0.15, 0.26)

ggsave(file="../manuscript/img/semivariogram_site_binom.pdf", 
  plot=semivariogram_site_binom, width=10, height=5)

# Tukeys's HSD test to see which groups are different

## Geog. Zone
tukey_geog_zone_no_family_binom <- emmeans(geog_rfamily_binom, "geog_zone")

tukey_binom_geog <- pwpp(tukey_geog_zone_no_family_binom, values = TRUE, sort = FALSE)

tukey_binom_geog_marg <- data.frame(
  y = tukey_binom_geog$layers[[3]]$data$geog_zone , 
  label = tukey_binom_geog$layers[[3]]$data$fmtval)

tukey_binom_geog_comp <- data.frame(
  x = tukey_binom_geog$data$p.value, 
  plus = tukey_binom_geog$data$plus, 
  minus = tukey_binom_geog$data$minus, 
  midpt = tukey_binom_geog$data$midpt)

tukey_binom_geog_ggplot <- ggplot(tukey_binom_geog_comp) + 
  geom_segment(
    aes(x = x, xend = x, y = plus, yend = midpt, colour = minus)) +
  geom_segment(
    aes(x = x, xend = x, y = minus, yend = midpt, colour = plus)) +
  geom_point(
    aes(x = x, y = plus, colour = minus), size = 3) + 
  geom_label(data = tukey_binom_geog_marg, 
    aes(x = 0.1, y = y, label = label),
    label.padding = unit(0.15, "lines")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour="#E0E0E0"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 12, colour = geog_zone_pal),
    legend.position = "none") + 
  ylab("Geog. Zone") +
  xlab("Tukey-adjusted P value") +
  scale_x_continuous(breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1)) +
  scale_colour_manual(values = geog_zone_pal) +
  coord_trans(x = "log10") 

ggsave(file="../manuscript/img/margin_binom_geog.pdf", plot=tukey_binom_geog_ggplot, width=10, height=8)

## Site
tukey_site_binom <- emmeans(site_rfamily_binom, "site_code")

tukey_binom_site <- pwpp(tukey_site_binom, values = TRUE) 

tukey_binom_site_marg <- data.frame(
  y = tukey_binom_site$layers[[3]]$data$site_code , 
  label = tukey_binom_site$layers[[3]]$data$fmtval)

tukey_binom_site_comp <- data.frame(
  x = tukey_binom_site$data$p.value, 
  plus = tukey_binom_site$data$plus, 
  minus = tukey_binom_site$data$minus, 
  midpt = tukey_binom_site$data$midpt)

tukey_binom_site_comp <- tukey_binom_site_comp %>%
  mutate(plus_colour = case_when(
    plus == "AM"| plus == "GE"| plus == "RD"| plus == "SO"| plus == "BE" ~ geog_zone_pal[1],
    plus == "LC"| plus == "SD"| plus == "GA"~ geog_zone_pal[2],
    plus == "GC"| plus == "AB" ~ geog_zone_pal[3],
    plus == "AC"| plus == "BB"| plus == "GD"| plus == "GT"| plus == "RM"| plus == "CG" ~ geog_zone_pal[4],
    plus == "GL"| plus == "BW" ~ geog_zone_pal[5],
    plus == "CC"| plus == "CR"| plus == "MG" ~ geog_zone_pal[6]
  ),
    minus_colour = case_when(
      minus == "AM"| minus == "GE"| minus == "RD"| minus == "SO"| minus == "BE" ~ geog_zone_pal[1],
      minus == "LC"| minus == "SD"| minus == "GA"~ geog_zone_pal[2],
      minus == "GC"| minus == "AB" ~ geog_zone_pal[3],
      minus == "AC"| minus == "BB"| minus == "GD"| minus == "GT"| minus == "RM"| minus == "CG" ~ geog_zone_pal[4],
      minus == "GL"| minus == "BW" ~ geog_zone_pal[5],
      minus == "CC"| minus == "CR"| minus == "MG" ~ geog_zone_pal[6]
    ))

tukey_binom_site_ggplot <- ggplot(tukey_binom_site_comp) + 
  geom_segment(
    aes(x = x, xend = x, y = plus, yend = minus), 
    colour = tukey_binom_site_comp$minus_colour, alpha = 0.5) +
  geom_segment(
    aes(x = x, xend = x, y = minus, yend = minus), 
    colour = tukey_binom_site_comp$plus_colour, alpha = 0.5) +
  geom_point(
    aes(x = x, y = plus), colour = tukey_binom_site_comp$plus_colour, size = 3) + 
  geom_label(data = tukey_binom_site_marg, 
    aes(x = 0.1, y = y, label = label),
    label.padding = unit(0.15, "lines")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour="#E0E0E0"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 12, colour = c(rep(geog_zone_pal[1], times = 5),
      rep(geog_zone_pal[2], times = 3), rep(geog_zone_pal[3], times = 2),
      rep(geog_zone_pal[4], times = 6), rep(geog_zone_pal[5], times = 2), rep(geog_zone_pal[6], times = 3))),
    legend.position = "none") + 
  ylab("Site") +
  xlab("Tukey-adjusted P value") +
  scale_x_continuous(breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1)) +
  scale_y_discrete(limits =  c("AM",
    "GE", "RD", "SO", "BE", "LC", "SD",
    "GA", "GC", "AB", "AC", "BB", "GD",
    "GT", "RM", "CG", "GL", "BW", "CC", "CR", "MG"), drop = FALSE) + 
  coord_trans(x = "log10") 

ggsave(file="../manuscript/img/margin_binom_site.pdf",
  plot=tukey_binom_site_ggplot, width=10, height=8)

# Mixed effects model - Non-zero glm ----

geog_rsite_rfamily_lmer <- glmmTMB(log(mm2_damage) ~ geog_zone + (1|site_name/family), 
  data = damage_nozero_dam)
geog_rsite_lmer <- glmmTMB(log(mm2_damage) ~ geog_zone + (1|site_name), 
  data = damage_nozero_dam)
site_rfamily_lmer <- glmmTMB(log(mm2_damage) ~ site_code + (1|family), 
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
  site_rgeog_lmer, site_rgeog_rfamily_lmer, null_lmer, 
  rand_rsite_rfamily_lmer, rand_rgeog_rsite_lmer)

lmer_aic_comp <- data.frame(
  fixed_effects = c("Geog. Zone", "Geog. Zone", "Site", "Site", "Site", "NA", "NA", "NA"),
  random_effects = c("Site / Parent", "Site", "Parent", "Geog. Zone", 
    "Geog. Zone + Parent", 
    "NA", "Site / Parent", "Geog. Zone / Site"),
  AIC = as.vector(unlist(lapply(lmer_list, AIC))),
  logLik = as.vector(unlist(lapply(lmer_list, logLik))),
  r2m = as.vector(unlist(lapply(lmer_list, r2m_mod_lmer))),
  r2c = as.vector(unlist(lapply(lmer_list, r2c_mod)))) %>%
  arrange(AIC)

fileConn<-file("../manuscript/include//lmer_aic_comp.tex")
writeLines(stargazer(lmer_aic_comp, 
  summary = FALSE, rownames = FALSE,
  label = "lmer_comp", digit.separate = 0), fileConn)
close(fileConn)

# Plots of predicted values

## Geog. zone
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
  theme(legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)) + 
  labs(x = "Geog. Zone", y = expression(paste("Predicted bark damage (mm"^2,")"))) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../manuscript/img/pred_lmer_geog.pdf", plot=pred_lmer_geog, width=10, height=8)

## Is there spatial autocorrelation amongst residuals in the best model?
damage_nozero_dam$geog_rsite_lmer_resid <- resid(geog_rsite_lmer)

geog_rsite_lmer_autofit = autofitVariogram(damage_nozero_dam$geog_rsite_lmer_resid~1, 
  damage_nozero_dam_sp, model = "Exp")

geog_rsite_lmer_vario <- variogram(damage_nozero_dam$geog_rsite_lmer_resid~1, 
  locations = ~x_coord + y_coord, data = damage_nozero_dam)

geog_rsite_lmer_vgm <- fit.variogram(geog_rsite_lmer_vario, 
  model = vgm(model = "Exp", 
    psill = geog_rsite_lmer_autofit$var_model$psill[2], 
    range = geog_rsite_lmer_autofit$var_model$range[2],
    kappa = geog_rsite_lmer_autofit$var_model$kappa[2],
    nugget = geog_rsite_lmer_autofit$var_model$psill[1]))

geog_rsite_lmer_vgm_fort <- variogramLine(geog_rsite_lmer_vgm, 
  maxdist = max(geog_rsite_lmer_vario$dist))

# Make semivariogram plot
semivariogram_geog_lmer <- ggplot() + 
  geom_point(data = geog_rsite_lmer_vario, 
    aes(x = dist, y = gamma)) + 
  geom_line(data = geog_rsite_lmer_vgm_fort, 
    aes(x = dist, y = gamma)) + 
  theme_classic() + 
  labs(x = "Distance (m)", y = "Semivariance (\u03B3)") + 
  geom_vline(aes(xintercept = geog_rsite_lmer_autofit$var_model$range[2]), linetype = 2)

ggsave(file="../manuscript/img/semivariogram_geog_lmer.pdf", 
  plot=semivariogram_geog_lmer, width=10, height=5)

## Site
site_predict_lmer <- ggpredict(site_rfamily_lmer, c("site_code"))
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
  theme(legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)) + 
  labs(x = "Geog. Zone and Seed population", 
    y = expression(paste("Predicted bark damage (mm"^2,")"))) + 
  theme(panel.spacing = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    strip.text = element_text(size = 12),
    legend.position = "none") + 
  scale_fill_manual(values = geog_zone_pal) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../manuscript/img/pred_lmer_site.pdf", 
  plot=pred_lmer_site, width=10, height=8)

## Is there spatial autocorrelation amongst residuals in the best model?
damage_nozero_dam$site_rfamily_lmer_resid <- resid(site_rfamily_lmer)

site_rfamily_lmer_autofit = autofitVariogram(damage_nozero_dam$site_rfamily_lmer_resid~1, 
  damage_nozero_dam_sp, model = "Sph")

site_rfamily_lmer_vario <- variogram(damage_nozero_dam$site_rfamily_lmer_resid~1, 
  locations = ~x_coord + y_coord, data = damage_nozero_dam)

site_rfamily_lmer_vgm <- fit.variogram(site_rfamily_lmer_vario, 
  model = vgm(model = "Sph", 
    psill = site_rfamily_lmer_autofit$var_model$psill[2], 
    range = site_rfamily_lmer_autofit$var_model$range[2],
    kappa = site_rfamily_lmer_autofit$var_model$kappa[2],
    nugget = site_rfamily_lmer_autofit$var_model$psill[1]))

site_rfamily_lmer_vgm_fort <- variogramLine(site_rfamily_lmer_vgm, 
  maxdist = max(site_rfamily_lmer_vario$dist))

# Make semivariogram plot
semivariogram_site_lmer <- ggplot() + 
  geom_point(data = site_rfamily_lmer_vario, 
    aes(x = dist, y = gamma)) + 
  geom_line(data = site_rfamily_lmer_vgm_fort, 
    aes(x = dist, y = gamma)) + 
  theme_classic() + 
  labs(x = "Distance (m)", y = "Semivariance (\u03B3)") + 
  geom_vline(aes(xintercept = site_rfamily_lmer_autofit$var_model$range[2]), linetype = 2)

ggsave(file="../manuscript/img/semivariogram_site_lmer.pdf", 
  plot=semivariogram_site_lmer, width=10, height=5)

# Tukeys's HSD test to see which groups are different.

## Geog. Zone
tukey_geog_zone_no_family_lmer <- emmeans(geog_rsite_lmer, "geog_zone")

tukey_lmer_geog <- pwpp(tukey_geog_zone_no_family_lmer, values = TRUE, sort = FALSE)

tukey_lmer_geog_marg <- data.frame(
  y = tukey_lmer_geog$layers[[3]]$data$geog_zone , 
  label = tukey_lmer_geog$layers[[3]]$data$fmtval)

tukey_lmer_geog_comp <- data.frame(
  x = tukey_lmer_geog$data$p.value, 
  plus = tukey_lmer_geog$data$plus, 
  minus = tukey_lmer_geog$data$minus, 
  midpt = tukey_lmer_geog$data$midpt)

tukey_lmer_geog_ggplot <- ggplot(tukey_lmer_geog_comp) + 
  geom_segment(
    aes(x = x, xend = x, y = plus, yend = midpt, colour = minus)) +
  geom_segment(
    aes(x = x, xend = x, y = minus, yend = midpt, colour = plus)) +
  geom_point(
    aes(x = x, y = plus, colour = minus), size = 3) + 
  geom_label(data = tukey_lmer_geog_marg, 
    aes(x = 0.008, y = y, label = label),
    label.padding = unit(0.15, "lines")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour="#E0E0E0"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 12, colour = geog_zone_pal),
    legend.position = "none") + 
  ylab("Geog. Zone") +
  xlab("Tukey-adjusted P value") +
  scale_x_continuous(breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1)) + 
  scale_colour_manual(values = geog_zone_pal) +
  coord_trans(x = "log10") 

ggsave(file="../manuscript/img/margin_lmer_geog.pdf", 
  plot=tukey_lmer_geog_ggplot, width=10, height=8)

## Site

## Reorder factor levels
site_rfamily_lmer$frame$site_code <- factor(
  site_rfamily_lmer$frame$site_code,levels = c("AM",
    "GE", "RD", "SO", "BE", "LC", "SD",
    "GA", "GC", "AB", "AC", "BB", "GD",
    "GT", "RM", "CG", "GL", "BW", "CC", "CR", "MG"))

tukey_site_lmer <- emmeans(site_rfamily_lmer, "site_code")

tukey_lmer_site <- pwpp(tukey_site_lmer, values = TRUE, sort = FALSE)

tukey_lmer_site_marg <- data.frame(
  y = tukey_lmer_site$layers[[3]]$data$site_code , 
  label = tukey_lmer_site$layers[[3]]$data$fmtval)

tukey_lmer_site_comp <- data.frame(
  x = tukey_lmer_site$data$p.value, 
  plus = tukey_lmer_site$data$plus, 
  minus = tukey_lmer_site$data$minus, 
  midpt = tukey_lmer_site$data$midpt)

tukey_lmer_site_comp <- tukey_lmer_site_comp %>%
  mutate(plus_colour = case_when(
    plus == "AM"| plus == "GE"| plus == "RD"| plus == "SO"| plus == "BE" ~ geog_zone_pal[1],
    plus == "LC"| plus == "SD"| plus == "GA"~ geog_zone_pal[2],
    plus == "GC"| plus == "AB" ~ geog_zone_pal[3],
    plus == "AC"| plus == "BB"| plus == "GD"| plus == "GT"| plus == "RM"| plus == "CG" ~ geog_zone_pal[4],
    plus == "GL"| plus == "BW" ~ geog_zone_pal[5],
    plus == "CC"| plus == "CR"| plus == "MG" ~ geog_zone_pal[6]
  ),
    minus_colour = case_when(
      minus == "AM"| minus == "GE"| minus == "RD"| minus == "SO"| minus == "BE" ~ geog_zone_pal[1],
      minus == "LC"| minus == "SD"| minus == "GA"~ geog_zone_pal[2],
      minus == "GC"| minus == "AB" ~ geog_zone_pal[3],
      minus == "AC"| minus == "BB"| minus == "GD"| minus == "GT"| minus == "RM"| minus == "CG" ~ geog_zone_pal[4],
      minus == "GL"| minus == "BW" ~ geog_zone_pal[5],
      minus == "CC"| minus == "CR"| minus == "MG" ~ geog_zone_pal[6]
    ))

tukey_lmer_site_ggplot <- ggplot(tukey_lmer_site_comp) + 
  geom_segment(
    aes(x = x, xend = x, y = plus, yend = minus), 
    colour = tukey_lmer_site_comp$plus_colour, alpha = 0.5) +
  geom_segment(
    aes(x = x, xend = x, y = plus, yend = minus), 
    colour = tukey_lmer_site_comp$plus_colour, alpha = 0.5) +
  geom_point(
    aes(x = x, y = plus), 
    colour = tukey_lmer_site_comp$plus_colour, size = 3) + 
  geom_label(data = tukey_lmer_site_marg, 
    aes(x = 0.03, y = y, label = label),
    label.padding = unit(0.15, "lines")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour="#E0E0E0"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 12, colour = c(rep(geog_zone_pal[1], times = 5),
      rep(geog_zone_pal[2], times = 3), rep(geog_zone_pal[3], times = 2),
      rep(geog_zone_pal[4], times = 6), rep(geog_zone_pal[5], times = 2), 
      rep(geog_zone_pal[6], times = 3))),
    legend.position = "none") + 
  ylab("Site") +
  scale_y_discrete(limits =  c("AM",
    "GE", "RD", "SO", "BE", "LC", "SD",
    "GA", "GC", "AB", "AC", "BB", "GD",
    "GT", "RM", "CG", "GL", "BW", "CC", "CR", "MG"), drop = FALSE) + 
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.5, 1)) + 
  coord_trans(x = "log10") 

ggsave(file="../manuscript/img/margin_lmer_site.pdf", plot=tukey_lmer_site_ggplot, width=10, height=8)

## Parent
tukey_family_rgeog_lmer <- emmeans(family_rgeog_lmer, "family")

tukey_lmer_family <- pwpp(tukey_family_rgeog_lmer, values = TRUE, sort = FALSE)

tukey_lmer_family_comp <- data.frame(
  x = tukey_lmer_family$data$p.value, 
  plus = tukey_lmer_family$data$plus, 
  minus = tukey_lmer_family$data$minus, 
  midpt = tukey_lmer_family$data$midpt)

# COLOUR THIS BY GEOG. ZONE
tukey_lmer_family_ggplot <- ggplot(tukey_lmer_family_comp) + 
  geom_segment(
    aes(x = x, xend = x, y = plus, yend = midpt, colour = minus)) +
  geom_segment(
    aes(x = x, xend = x, y = minus, yend = midpt, colour = plus)) +
  geom_point(
    aes(x = x, y = plus, colour = minus), size = 3) + 
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour="#E0E0E0"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_text(size = 6),
    legend.position = "none") + 
  ylab("Parent") +
  xlab("Tukey-adjusted P value") +
  scale_x_continuous(breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1)) + 
  coord_trans(x = "log10") 

ggsave(file="../manuscript/img/margin_lmer_family.pdf", 
  plot=tukey_lmer_family_ggplot, width=5, height=12)

# Linear mixed effects model of latitude and longitude against total damage

dec_lat_lmer <- glmmTMB(log(mm2_damage) ~ dec_latitude + (1|site_name/family), 
  data = damage_nozero_dam)
dec_lon_lmer <- glmmTMB(log(mm2_damage) ~ dec_longitude + (1|site_name/family), 
  data = damage_nozero_dam)

summary(dec_lat_lmer)
summary(dec_lon_lmer)

test <- Anova(dec_lat_lmer)

r.squaredGLMM(dec_lat_lmer)
r.squaredGLMM(dec_lon_lmer)

dec_lat_predict <- ggpredict(dec_lat_lmer, terms = c("dec_latitude"), type = "fe")

# Plot predicted values by fixed effect
pred_lat <- ggplot() + 
  geom_ribbon(aes(x = x, ymin = log(conf.low), ymax = log(conf.high)),
    data = dec_lat_predict, alpha = 0.3) + 
  geom_point(data = damage_nozero_dam, 
    aes(x = dec_latitude, y = log(mm2_damage), fill = geog_zone), 
    position = "dodge", shape = 21, colour = "black", size = 2) + 
  geom_line(data = dec_lat_predict, 
    aes(x = x, y = log(predicted))) + 
  scale_fill_manual(values = geog_zone_pal, name = "Geog. zone") + 
  theme_classic() + 
  labs(x = "Latitude", y = expression(paste("Log bark damage (mm"^2,")"))) 

ggsave(file="../manuscript/img/pred_lat.pdf", plot=pred_lat, width=10, height=5)

# Grid arrange all semivariograms
pdf(file = "../manuscript/img/semivariogram_all.pdf", width=10, height=5)
grid.arrange(
  semivariogram_site_binom + ggtitle("Site best binomial model"),
  semivariogram_geog_zone_binom + ggtitle("Geog. zone best binomial model"),
  semivariogram_site_lmer + ggtitle("Site best damage model"),
  semivariogram_geog_lmer + ggtitle("Geog. zone best damage model"),
  nrow = 2)
dev.off()


# Models of chemical fingerprint vs. damage ----

# Create dataset of chem data only
damage_chem <- damage_nozero_dam[,c(23:87, 89:101)] %>% 
  gather(key = "chem", value = "value")

damage_chem$mm2_damage <- rep(damage_nozero_dam$mm2_damage, 
  times = length(damage_chem))

all_biochem <- ggplot(damage_chem, aes(x = value, y = log(mm2_damage))) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~chem, scales = "free_x")

pdf("all_biochem_regression.pdf", width = 15, height = 15)
all_biochem
dev.off()

biochem_single <- function(name){
  name <- enquo(name)
  ggplot(data = damage_nozero_dam, 
    aes(x = !!name, y = log(mm2_damage))) + 
    geom_point() + 
    geom_smooth(method = "lm")
}

biochem_single(P_all)
biochem_single(TU_all)
biochem_single(mois_per)
biochem_single(dry_mass_per)
biochem_single(a_humulene)
biochem_single(phen_mgg)
biochem_single(ox_phen_mgg)
biochem_single(fructose_mgg)
biochem_single(sugar_mgg)
biochem_single(a_pinene_m)
biochem_single(a_pinene_p)

chem_corr_df <- damage_nozero_dam %>%
  dplyr::select(starts_with("terp_", ignore.case = FALSE), 
    starts_with("TU", ignore.case = FALSE), 
    starts_with("P", ignore.case = FALSE), -TU_all, -P_all) %>%
  na.omit() 

pdf("/Users/johngodlee/Desktop/test.pdf", width = 25, height = 25)
ggcorrplot(cor(chem_corr_df, use = "complete.obs"), 
  type = "lower", 
  lab = TRUE,
  method = "square",
  colors = c("blue", "white", "red"),
  ggtheme = theme_classic, 
  show.legend = FALSE,
  outline.color = "black")
dev.off()

damage_nozero_dam_std <- damage_nozero_dam %>%
  dplyr::select(mm2_damage, 
    starts_with("P", ignore.case = FALSE),
    starts_with("terp_", ignore.case = FALSE), 
    starts_with("TU", ignore.case = FALSE), 
    -P_all, -TU_all) %>%
  mutate_at(vars(-one_of(c("mm2_damage"))),
    .funs = list(std = ~(as.vector(scale(.,))))) %>%
  dplyr::select(ends_with("_std"), mm2_damage)

all_chem_lm <- lm(log(damage_nozero_dam$mm2_damage) ~ ., damage_nozero_dam_std)

all_chem_lm_summ <- summary(all_chem_lm)

all_chem_lm_summ_df <- data.frame(
  var = row.names(all_chem_lm_summ$coefficients),
  slope = all_chem_lm_summ$coefficients[, 1], 
  std_err = all_chem_lm_summ$coefficients[, 2],
  t_value = all_chem_lm_summ$coefficients[, 3],
  p_value = all_chem_lm_summ$coefficients[, 4])
  
  
all_chem_lm_summ_df %>% 
  arrange(p_value) 

# Run a Principal component analysis on the phenolics and terpenes to get principal components
damage_nozero_dam_naomit <- damage_nozero_dam %>%
  dplyr::select(starts_with("terp_", ignore.case = FALSE), 
    starts_with("TU", ignore.case = FALSE),
    starts_with("P", ignore.case = FALSE), 
    site_name, seed_zone, geog_zone, site_code, big_region, 
    dec_latitude, dec_longitude, family, field_code, mm2_damage,
    sugar_mgg,
    -TU_all, -P_all) %>%
  na.omit()

damage_nozero_dam_naomit[rowSums(is.na(damage_nozero_dam_naomit)) > 0,]

damage_pca <- damage_nozero_dam_naomit %>%
  dplyr::select(starts_with("terp_", ignore.case = FALSE), 
    starts_with("TU", ignore.case = FALSE),
    starts_with("P", ignore.case = FALSE)) %>%
  rda(., scale = TRUE)

# Principal component axes variance explained
as.vector(damage_pca$CA$eig)/sum(damage_pca$CA$eig)
barplot(as.vector(damage_pca$CA$eig)/sum(damage_pca$CA$eig)) 

# How much explained by first two axes
sum((as.vector(damage_pca$CA$eig)/sum(damage_pca$CA$eig))[1:3])
##' 50% isn't great

# Plot distribution of sites and "species" (chemicals)
biplot(damage_pca)
plot(damage_pca, display = "sites", type = "points")
plot(damage_pca, display = "species", type = "text")

# Extract site scores
sitePCA <- damage_pca$CA$u 

damage_nozero_dam_naomit$chem_pca_1 <- sitePCA[,1]
damage_nozero_dam_naomit$chem_pca_2 <- sitePCA[,2]
damage_nozero_dam_naomit$chem_pca_3 <- sitePCA[,3]

# Run linear model of principal components against mm2 damage
chem_pca_lm <- lm(log(mm2_damage) ~ chem_pca_1 + chem_pca_2 + chem_pca_3, 
  data = damage_nozero_dam_naomit)

summary(chem_pca_lm)

# Linear model of alpha-pinene against mm2 damage
pinene_m_lm <- lm(log(mm2_damage) ~ terp_a_pinene_m, data = damage_nozero_dam_naomit)

summary(pinene_m_lm)

# Linear model of sugar content against mm2 damage
sugar_lm <- lm(log(mm2_damage) ~ sugar_mgg, data = damage_nozero_dam_naomit)

summary(sugar_lm)
