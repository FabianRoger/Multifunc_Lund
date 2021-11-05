
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Illustrate the single function approach using the Jena data

# load relevant libraries
library(here)
library(readr)
library(dplyr)
library(ggplot2)

rm(list = ls())

# link to scripts with the relevant functions
source(here("Scripts/function_plotting_theme.R"))

# read in the Jena data
jena.dat <- read_csv(file = here("data/jena_data_Jochum_2020_clean.csv"))

# check the column names
col_names <- names(jena.dat)

# remove the individual species abundances
spp <- col_names[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", col_names)]
length(spp)
rm(col_names)

jena.dat <- 
  jena.dat %>%
  select(-all_of(spp))

# choose four ecosystem functions to examine in the single function approach
names(jena.dat)
func.names <- c("biomass", "soilorgC", "poll", "rootBM")

jena.dat <- 
  jena.dat %>%
  select(plotcode, block, sowndiv, realised_richness1, all_of(func.names) )

# check the summary statistics and distributions
summary(jena.dat)
length(unique(jena.dat$block))
lapply(jena.dat[, func.names], hist)

# scale the sowndiv as this is our predictor
jena.dat <- 
  jena.dat %>%
  mutate(sowndiv_scale = scale(sowndiv, center = TRUE, scale = TRUE)[,1])

# fit models with each variable as a response and sowndiv as a predictor

# biomass
hist(sqrt(jena.dat$biomass) )
lm.bm <- glm((biomass) ~ sowndiv_scale, data = jena.dat, family = "gaussian")
plot(lm.bm) # assumptions look good

summary(lm.bm)

# soilC
hist(jena.dat$soilorgC)
lm.sc <- glm(soilorgC ~ sowndiv_scale, data = jena.dat, family = "gaussian")
plot(lm.sc)

summary(lm.sc)

# pollination
hist(jena.dat$poll)
lm.po <- glm(round(poll, 0) ~ sowndiv_scale, data = jena.dat, family = "poisson")
plot(lm.po)

summary(lm.po)

# pollination
hist(jena.dat$rootBM)
lm.po <- glm(sqrt(rootBM) ~ sowndiv_scale, data = jena.dat, family = "gaussian")
plot(lm.po)

summary(lm.po)

# generate a correlation plot between these functions



