
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Add multifunctionality to the abundance data from the stachova and leps (2010) model

# select scripts to call functions from
library(here)
source(here("Scripts/multispecies_lotka_volterra_model/stachova_leps_2010_LK_model_function.R"))

# simulate a set of communities
df.1 <- s_l_function(lsp = c(4, 8, 12, 16, 20),
                     reps = 12,
                     rsp = 30,
                     t_steps = 100,
                     n0 = 5,
                     ext.thresh = 0.2,
                     a_mean = 0.25, a_sd = 0.1, a_min = 0, a_max = 1.2, 
                     a_spp = 1, sim.comp = "sym",
                     k_min = 5, k_max = 150,
                     r_min = 0.1, r_max = 0.5)

# check basic relationships
library(ggplot2)
library(dplyr)

df.sum <- 
  df.1[[1]] %>%
  filter(total_abundance > 0)

ggplot(data = df.sum,
       mapping = aes(x = richness, y = total_abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

ggplot(data = df.sum,
       mapping = aes(x = local.sp.pool, y = total_abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# check the species parameters
df.1[[3]]

# get the species abundance data
df.spp <- df.1[[2]]

# make a species list
spp.list <- unique(df.spp$species)

# choose the number of functions
func.n <- 9

# define the distribution parameters

# select from either -0.5 and 0 so that some functions are never negatively associated with abundances
# this is realistic (e.g. productivity)
func.min <- sample(c(-0.5, 0), func.n, replace = TRUE)
func.max <- rep(1, func.n)

func.mat <- 
  lapply(spp.list, function(x){
  z <- runif(n = func.n, min = func.min, max = func.max)
  round(z, digits = 4)
})

func.mat <- data.frame(do.call(rbind, func.mat))
func.mat <- cbind(spp.list, func.mat)

names(func.mat) <- c("species", paste("F", 1:func.n, sep = "_"))

multi.func <- full_join(df.spp, func.mat, by = "species")


# multiply abundance by these function values
multi.func %>%
  mutate(across(.cols = starts_with(match = "F_"), ~(.*abundance) ))





