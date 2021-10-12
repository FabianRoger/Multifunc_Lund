
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Explore the simulated data and the Jena data for mistakes, realism etc.

# load relevant libraries
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# scripts to draw functions from
source(here("Scripts/function_plotting_theme.R"))
source(here("Scripts/MF_functions_collated.R"))


# simulated data

# load the simulation parameters
par.dat <- read_csv(file = here("data/parameters_sim.csv"))
head(par.dat)

# load the raw multifunctionality data
sim.raw <- read_csv(file = here("data/multifunctionality_data.csv"))
head(sim.raw)

# load the raw abundance data
abun.raw <- read_csv(file = here("data/raw_abundance_data.csv"))
head(abun.raw)


# choose a sim.id to examine
id.in <- 13

# plot the relative species abundance distributions of the plots
# plot rank-abundance distributions for four randomly chosen patches
p.df <- 
  abun.raw %>%
  filter(sim.id == id.in) %>%
  filter(abundance > 0) %>%
  mutate(patch = as.character(patch))

# get four random patches
ran.patch <- sample(x = unique(p.df$patch), size = 5)

p.df.ran <- 
  p.df %>%
  filter(patch %in% ran.patch) %>%
  group_by(patch) %>%
  mutate(abundance_rank = rank(-abundance)) %>%
  ungroup()

ggplot(data = p.df.ran,
         mapping = aes(x = (abundance_rank), y = abundance, colour = patch)) +
  stat_smooth(method='lm', formula = y~poly(x, 2), se = FALSE, size = 1) +
  geom_jitter(width = 0.1) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  scale_x_continuous(breaks = c(0:max(p.df.ran$abundance_rank))) +
  ylab("abundance") +
  xlab("rank") +
  theme_meta() +
  theme(legend.position = "none") 


# plot the relationship between species pool diversity and function
dat.in <- 
  sim.raw %>%
  filter(sim.id == id.in)

ggplot(data = dat.in,
       mapping = aes(x = local.sp.pool, y = total_abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_meta()

ggplot(data = dat.in,
       mapping = aes(x = richness, y = total_abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_meta()


# plot the relationship between richness and ecosystem function for each function
dat.in %>%
  pivot_longer(cols = starts_with("F_"),
               names_to = "function_id",
               values_to = "value") %>%
  ggplot(data = .,
         mapping = aes(x = richness, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_meta() +
  facet_wrap(~function_id, scales = "free")


# plot the relationship between richness and average multifunctionality
ggplot(data = dat.in,
       mapping = aes(x = richness, y = ave._MF)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_meta()


# plot the relationship between richness and summing MF
ggplot(data = dat.in,
       mapping = aes(x = richness, y = sum_MF)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_meta()

# choose two random functions and plot sum_MF again
fs <- sample(paste("F_", 1:9, sep = ""), 2)

dat.in %>%
  mutate(MF_sum2 = MF_sum(adf = dat.in, vars = fs) ) %>%
  ggplot(data = .,
         mapping = aes(x = richness, y = MF_sum2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_meta()








