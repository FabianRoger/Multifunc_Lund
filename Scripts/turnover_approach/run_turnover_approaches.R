
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Run the turnover approaches on the Jena data

# load relevant libraries
library(here)
library(readr)
library(dplyr)
library(ggplot2)

rm(list = ls())

# link to scripts with the relevant functions
source(here("Scripts/turnover_approach/turnover_approach_functions.R"))
source(here("Scripts/function_plotting_theme.R"))
source(here("Scripts/process_model_data.R"))

# read in the Jena data
jena.dat <- read_csv(file = here("data/jena_data_Jochum_2020_clean.csv"))

# get a vector of species names
col_names <- names(jena.dat)

spp <- col_names[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", col_names)]
length(spp)
rm(col_names)

# check if all species are present in at least one plot
any(colSums(jena.dat[, spp]) == 0)
min(colSums(jena.dat[, spp]))

# convert the species abundances into presence-absence data
jena.dat <- 
  jena.dat %>%
  mutate(across(.cols = all_of(spp), ~if_else(. > 0, 1, 0)))

# get a vector of function names
names(jena.dat)
func.names <- c("biomass", "plantCN", "soilC", "soilorgC", "herbi",
                "micBMC", "phoact", "poll", "rootBM")

# get the random expectation for the jena data based on the AIC approach
jena.aic.ran <- 
  prop_species_pool_random(data = jena.dat, 
                           function_names = func.names, 
                           species_names = spp, 
                           method = "AIC", n_ran = 100, n = 1000
                           )

jena.aic.obs <- 
  prop_species_pool(data = jena.dat, 
                    function_names = func.names, 
                    species_names = spp, 
                    method = "SES", n_ran = 100
                    ) 


# plot the results

# check the max of the null data to see if it is below 1
max(jena.aic.ran$prop_species_pool)

# calculate the 97.5% and 2.5% quantiles for each number of functions
# plot this as well
jena.aic.ran.q <- 
  jena.aic.ran %>% 
  group_by(run, number_functions) %>%
  summarise(prop_species_pool_m = mean(prop_species_pool, na.rm = TRUE), .groups = "drop") %>%
  group_by(number_functions) %>%
  summarise(quant_97.5 = quantile(prop_species_pool_m, probs = c(0.975)),
            quant_2.5 = quantile(prop_species_pool_m, probs = c(0.025)),
            .groups = "drop")

# calculate mean and sd of the observed data
jena.aic.obs.m <- 
  jena.aic.obs %>%
  group_by(number_functions) %>%
  summarise(prop_species_pool_m = mean(prop_species_pool, na.rm = TRUE),
            prop_species_pool_sd = sd(prop_species_pool_m, na.rm = TRUE) )

# plot the null expectations versus the observed data
library(ggplot2)
g1 <- 
  ggplot() +
  geom_ribbon(data = jena.aic.ran.q,
              mapping = aes(x = number_functions, ymax = quant_97.5, ymin = quant_2.5),
              alpha = 0.25) +
  geom_jitter(data = jena.aic.obs,
              mapping = aes(x = number_functions, y = prop_species_pool), 
              width = 0.1, colour = "red", shape = 16, alpha = 0.1, size = 2) +
  geom_line(data = jena.aic.obs.m,
            mapping = aes(x = number_functions, y = prop_species_pool_m), size = 1, colour = "red") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  ylab("number of species") +
  xlab("number of functions") +
  theme_meta()

g1

# gsave(filename = here("Figures/aic_turnover_null.png"), plot = g1,
       # width = 11, height = 7.5, units = "cm", dpi = 450)

# randomly assigning functions to plots generates the same pattern as the empirical data
# why is this the case?

# maybe it is because the turnover approach is anti-conservative...



### END
