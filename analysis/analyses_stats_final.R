# Statistical modelling final version 
# John Godlee (johngodlee@gmail.com)
# 2019_08_07

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

ggsave(file="../paper/img/variogram.pdf", plot=semivariogram, width=10, height=5)


##' Semivariogram indicates that there isn't 
##' any great amount of spatial autocorrelation

# Statistical modelling ----

# Mixed effects model - logistic regression

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
  site_rgeog_rfamily_binom,  family_rgeog_rsite_binom, 
  family_rgeog_binom, null_binom, 
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
  fixed_effects = c("Geog. Zone", "Geog. Zone", "Site", "Site", "Site", "Parent", 
    "Parent", "NA", "NA", "NA", "NA", "NA"),
  random_effects = c("Parent", "Site / Parent", "Parent", "Geog. Zone", "Geog. Zone + Parent", 
    "Geog. Zone / Site", "Geog. Zone", 
    "NA", "Parent", "Site", "Site / Parent", "Geog. Zone / Site / Parent"),
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

ggsave(file="../paper/img/pred_binom_geog.pdf", plot=pred_binom_geog, width=10, height=8)

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

ggsave(file="../paper/img/pred_binom_site.pdf", plot=pred_binom_site, width=10, height=8)

# Tukeys's HSD test to see which groups are different.

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

ggsave(file="../paper/img/margin_binom_geog.pdf", plot=tukey_binom_geog_ggplot, width=10, height=8)


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

ggsave(file="../paper/img/margin_binom_site.pdf", plot=tukey_binom_site_ggplot, width=10, height=8)

# Mixed effects model - Non-zero glm

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
  site_rgeog_lmer, site_rgeog_rfamily_lmer, family_rgeog_rsite_lmer, 
  family_rgeog_lmer, null_lmer, rand_rsite_rfamily_lmer, rand_rgeog_rsite_lmer)

lmer_aic_comp <- data.frame(
  fixed_effects = c("Geog. Zone", "Geog. Zone", "Site", "Site", "Site",
    "Parent", "Parent", "NA", "NA", "NA"),
  random_effects = c("Site / Parent", "Site", "Parent", "Geog. Zone", 
    "Geog. Zone + Parent", "Geog. Zone / Site", "Geog. Zone", 
    "NA", "Site / Parent", "Geog. Zone / Site"),
  AIC = as.vector(unlist(lapply(lmer_list, AIC))),
  logLik = as.vector(unlist(lapply(lmer_list, logLik))),
  r2m = as.vector(unlist(lapply(lmer_list, r2m_mod_lmer))),
  r2c = as.vector(unlist(lapply(lmer_list, r2c_mod)))) %>%
  arrange(AIC)

fileConn<-file("../paper/include/lmer_aic_comp.tex")
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

ggsave(file="../paper/img/pred_lmer_geog.pdf", plot=pred_lmer_geog, width=10, height=8)

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
  labs(x = "Geog. Zone and Seed population", y = expression(paste("Predicted bark damage (mm"^2,")"))) + 
  theme(panel.spacing = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    strip.text = element_text(size = 12),
    legend.position = "none") + 
  scale_fill_manual(values = geog_zone_pal) + 
  scale_colour_manual(values = geog_zone_pal)

ggsave(file="../paper/img/pred_lmer_site.pdf", plot=pred_lmer_site, width=10, height=8)

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

ggsave(file="../paper/img/margin_lmer_geog.pdf", plot=tukey_lmer_geog_ggplot, width=10, height=8)

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
      rep(geog_zone_pal[4], times = 6), rep(geog_zone_pal[5], times = 2), rep(geog_zone_pal[6], times = 3))),
    legend.position = "none") + 
  ylab("Site") +
  scale_y_discrete(limits =  c("AM",
    "GE", "RD", "SO", "BE", "LC", "SD",
    "GA", "GC", "AB", "AC", "BB", "GD",
    "GT", "RM", "CG", "GL", "BW", "CC", "CR", "MG"), drop = FALSE) + 
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.5, 1)) + 
  coord_trans(x = "log10") 

ggsave(file="../paper/img/margin_lmer_site.pdf", plot=tukey_lmer_site_ggplot, width=10, height=8)

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

ggsave(file="../paper/img/margin_lmer_family.pdf", plot=tukey_lmer_family_ggplot, width=5, height=12)


# Linear mixed effects model of latitude and longitude against total damage

dec_lat_lmer <- glmmTMB(log(mm2_damage) ~ dec_latitude + (1|site_name/family), data = damage_nozero_dam)
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

ggsave(file="../paper/img/pred_lat.pdf", plot=pred_lat, width=10, height=5)

