
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Run the turnover approaches on the Jena data

# load relevant libraries
library(here)
library(readr)
library(dplyr)
library(ggplot2)

rm(list = ls())

# link to scripts with the relevant functions
source(here("Scripts/05_turnover_approach/turnover_approach_functions.R"))
source(here("Scripts/01_general_functions/function_plotting_theme.R"))
source(here("Scripts/01_general_functions/process_model_data.R"))

# read in the Jena data
jena.dat <- read_csv(file = here("data/jena_data_Jochum_2020_clean.csv"))

# get a vector of species names
col_names <- names(jena.dat)

spp <- col_names[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", col_names)]
length(spp)
rm(col_names)

# check if all species are present in at least one plot
any(colSums(jena.dat[, spp]) == 0) # if FALSE then all species are present in at least one plot
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
                           method = "AIC", n_ran = 100, n = 100
                           )
jena.aic.ran$method <- "AIC"

jena.aic.obs <- 
  prop_species_pool(data = jena.dat, 
                    function_names = func.names, 
                    species_names = spp, 
                    method = "AIC", n_ran = 100
                    ) 
jena.aic.obs$method <- "AIC"

# get the random expectation for the jena data based on the SES approach
jena.ses.ran <- 
  prop_species_pool_random(data = jena.dat, 
                           function_names = func.names, 
                           species_names = spp, 
                           method = "SES", n_ran = 100, n = 100
  )
jena.ses.ran$method <- "SES"

jena.ses.obs <- 
  prop_species_pool(data = jena.dat, 
                    function_names = func.names, 
                    species_names = spp, 
                    method = "SES", n_ran = 100
  ) 
jena.ses.obs$method <- "SES"

# bind these datasets together
jena.ran <- bind_rows(jena.aic.ran, jena.ses.ran)
jena.obs <- bind_rows(jena.aic.obs, jena.ses.obs)

# plot the results

# check the max of the null data to see if it is below 1
max(jena.ran$proportion_species_pool)

# calculate the 97.5% and 2.5% quantiles for each number of functions
# plot this as well
jena.ran.q <- 
  jena.ran %>% 
  group_by(method, run, effect_direction, number_functions) %>%
  summarise(proportion_species_pool_m = mean(proportion_species_pool, na.rm = TRUE), .groups = "drop") %>%
  group_by(method, effect_direction, number_functions) %>%
  summarise(quant_97.5 = quantile(proportion_species_pool_m, probs = c(0.975)),
            quant_2.5 = quantile(proportion_species_pool_m, probs = c(0.025)),
            .groups = "drop")
head(jena.ran.q)


# calculate mean and sd of the observed data
jena.obs.m <- 
  jena.obs %>%
  group_by(method, effect_direction, number_functions) %>%
  summarise(proportion_species_pool_m = mean(proportion_species_pool, na.rm = TRUE),
            proportion_species_pool_sd = sd(proportion_species_pool, na.rm = TRUE), .groups = "drop" )

# plot the null expectations versus the observed data
library(ggplot2)

l.in <- list(jena.ran.q, jena.obs, jena.obs.m)

l.in <- 
  lapply(l.in, function(x) {
  
  x %>%
    mutate(par_id = paste(paste(method, ":", sep = ""), 
                          gsub(x = effect_direction, pattern = "_", replacement = " "), sep = " " ) )
  
})

g1 <- 
  ggplot() +
  geom_ribbon(data = l.in[[1]],
              mapping = aes(x = number_functions, ymax = quant_97.5, ymin = quant_2.5),
              alpha = 0.15) +
  geom_jitter(data = l.in[[2]],
              mapping = aes(x = number_functions, y = proportion_species_pool), 
              width = 0.1, colour = "red", shape = 1, alpha = 0.1, size = 2) +
  geom_line(data = l.in[[3]],
            mapping = aes(x = number_functions, y = proportion_species_pool_m), size = 1, colour = "red") +
  geom_ribbon(data = l.in[[3]],
              mapping = aes(x = number_functions, 
                            ymax = proportion_species_pool_m + proportion_species_pool_sd, 
                            ymin = proportion_species_pool_m - proportion_species_pool_sd),
              alpha = 0.2, fill = "red") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  ylab("Proportion species pool") +
  xlab("Number of functions") +
  facet_wrap(~par_id) +
  theme_meta() + 
  theme(axis.text = element_text(size = 9),
    strip.text = element_text(
      size = 11
    ))

g1

ggsave(filename = here("Figures/turnover_null.png"), plot = g1,
       width = 12, height = 11, units = "cm", dpi = 450)

### END
