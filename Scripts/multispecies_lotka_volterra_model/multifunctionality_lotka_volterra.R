
# Project: Mutifunctionality workshop (Lund 2019)

# Title: Adding multifunctionality into Lotka-Volterra models

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)
library(corrplot)

# choose scripts to draw functions from
source(here("Scripts/MF_functions_collated.R"))
source(here("Scripts/function_plotting_theme.R"))
source(here("Scripts/Multifunctionality-Simulations/Multifunc_simulations_functions.R"))

# load the Lotka-Volterra simulated data
lv_dat <- 
  read_csv(file = here("Scripts/multispecies_lotka_volterra_model/lv_mf_sim_dat.csv"))

# check dimensions of the data
nrow(lv_dat)

# make a vector of species names
sp_names <- unique(lv_dat$species)

# choose the number of functions to simulate
func_n <- 9

# set the function names
func_names <- paste0("F_", 1:func_n)

# for each species, assign a function value for x functions

# create an empty data matrix
Funcs <- 
  matrix(
  nrow = length(sp_names),
  ncol = func_n,
  dimnames = list(sp_names, func_names))


# we can fill this empty matrix in different ways:

# fill this with function values with different correlation levels among species

# number of species
specnum <- length(sp_names)

# choose pairwise correlation strength
COR <- 0

# make correlation matrix (strictly speaking a covariance matrix but for these simulations it does not matter)
Sigma <- matrix(COR, ncol = func_n, nrow = func_n)

# make three 'cluster' of correlated functions
Sigma[1:4,1:4] <- -2
Sigma[5:6,5:6] <- -2
Sigma[7:9,7:9] <- -2

diag(Sigma) <- rnorm(n = length(diag(Sigma)), mean = 10, sd = 2)

Sigma

# draw correlated functions (with mean 0)
corF <- mvrnorm(n = specnum, mu = rep(0, func_n), Sigma = Sigma)

# shift to positive
corF <- apply(corF, 2, function(x){ x + abs(min(x)) })

# fill the function matrix
Funcs[1:nrow(Funcs), 1:ncol(Funcs)] <- 
  as.vector(corF)

Funcs

Funcs %>%
  cor() %>%
  corrplot(method = "ellipse")


# convert the Funcs data to a dataframe
Funcs <- as.data.frame(Funcs)

# add a species column
Funcs$species <- row.names(Funcs)
row.names(Funcs) <- NULL


# extract the abundance of each species in the last time point
# multiply each species abundance by the function value
# calculate function scores for each species

bio_funcs <- 
  left_join(filter(lv_dat, time == last(time)),
            Funcs,
            by = "species") %>%
  mutate( across(.cols = all_of(func_names), ~(.*abundance) ) ) %>%
  group_by(replicate, species_pool) %>%
  summarise( across(.cols = c("abundance", all_of(func_names) ) , ~sum(.) ), .groups = "drop" )

# z-score standardise the different functions
# translate them to make sure they are positive
bio_funcs <- 
  bio_funcs %>%
  mutate(across(.cols = all_of(func_names), ~as.numeric(scale(.x, center = TRUE, scale = TRUE)) )) %>%
  mutate(across(.cols = all_of(func_names), ~(.x + abs(min(.x)))  ))

bio_funcs %>%
  select(abundance, contains("F")) %>%
  cor() %>%
  corrplot(method = "ellipse")


# write a function to output a slope between richness and function for different numbers of functions
est_n_func <- function(data, funcs, n_functions) {
  
  # set up combinations for each number of functions
  comb_f <- gtools::combinations(n = length(funcs), r = n_functions)
  
  # set up an output vector
  y <- vector("list", length = nrow(comb_f))
  
  for (i in 1:nrow(comb_f)) {
    
    v <- funcs[comb_f[i, ]]
    
    df <- data[, v]
    
    lm_df <- 
      data %>%
      mutate(species_pool = scale(species_pool, scale = TRUE, center = TRUE),
             ave_mf = MF_av(adf = df, vars = v),
             pasari_mf = MF_pasari(adf = df, vars = v),
             enf_mf = hill_multifunc(adf = df, vars = v, scale = 1, HILL = TRUE),
             thresh_30_mf = single_threshold_mf(adf = df, vars = v, thresh = 0.3)/length(v),
             thresh_50_mf = single_threshold_mf(adf = df, vars = v, thresh = 0.5)/length(v),
             thresh_70_mf = single_threshold_mf(adf = df, vars = v, thresh = 0.7)/length(v)) %>%
      mutate(sd_ave = sd(ave_mf),
             mean_ave = mean(ave_mf))
    
    mf_names <- names(select(lm_df, ends_with("mf")))
    
    s_cof <- vector(length = length(mf_names))
    
    for(s in 1:length(mf_names)) {
      
      lm1 <- lm(as.formula(paste0(mf_names[s], "~ species_pool")),
                data = lm_df)
      
      s_cof[s] <- lm1$coefficients[2]
      
    }
    
    s_cof <- c(s_cof, lm_df$sd_ave[1], lm_df$mean_ave[1])
    
    names(s_cof) <- c(mf_names, "sd_ave", "mean_ave")
    
    y[[i]] <- s_cof
    
  } 
  
  y <- bind_rows(y, .id = "run")
  
  y$number_functions <- n_functions
  
  y
  
}

# test this function
est_n_func(data = bio_funcs, funcs = func_names, n_functions = 2)


# loop over this function for all possible numbers of functions
slopes_n_func <- vector("list", length = length(func_names))

for(j in 2:length(func_names)) {
  
  slopes_n_func[[j]] <- est_n_func(data = bio_funcs, funcs = func_names, n_functions = j)
  
}

slopes_n_func <- bind_rows(slopes_n_func, .id = "n_func")

# plot the diversity-function relationship
bio_funcs %>%
  pivot_longer(cols = starts_with("F"),
               names_to = "eco_function",
               values_to = "value") %>%
  ggplot(data = .,
         mapping = aes(x = species_pool, y = value)) +
  geom_jitter(width = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~eco_function, scales = "free") +
  theme_meta()

# plot the standardised effect size between species pool diversity and function
slopes_n_func <- 
  slopes_n_func %>%
  pivot_longer(cols = ends_with("mf"),
               names_to = "mf_metric",
               values_to = "mf_value")

ggplot(data = slopes_n_func, 
       mapping = aes(x = number_functions, mf_value)) +
  geom_jitter(width = 0.5, alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~mf_metric, scales = "free") +
  theme_meta()

ggplot(data = slopes_n_func %>%
         pivot_longer(cols = c("sd_ave", "mean_ave"),
                      names_to = "met",
                      values_to = "val"),
       mapping = aes(x = number_functions, y = val)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  facet_wrap(~met, scales = "free") +
  theme_meta()




















