
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

# function to get the slope for different numbers of functions in a dataset
# data is a dataset with multiple functions as columns
# funcs is a vector of function names
est_n_func <- function(data, funcs) {
  
  s_est <- function(a, b, d) {
    
    # set up combinations for each number of functions
    comb_f <- gtools::combinations(n = length(b), r = d)
    
    # set up an output vector
    y <- vector("list", length = nrow(comb_f))
    
    for (i in 1:nrow(comb_f)) {
      
      v <- b[comb_f[i, ]]
      
      df <- a[, v]
      
      lm_df <- 
        a %>%
        mutate(species_pool = scale(species_pool, scale = TRUE, center = TRUE),
               ave_mf = MF_av(adf = df, vars = v),
               sum_mf = MF_sum(adf = df, vars = v),
               pasari_mf = MF_pasari(adf = df, vars = v),
               enf_mf = hill_multifunc(adf = df, vars = v, scale = 1, HILL = TRUE),
               thresh_30_mf = single_threshold_mf(adf = df, vars = v, thresh = 0.3),
               thresh_70_mf = single_threshold_mf(adf = df, vars = v, thresh = 0.7))
      
      mf_names <- names(select(lm_df, ends_with("mf")))
      
      s_cof <- vector(length = length(mf_names))
      
      for(s in 1:length(mf_names)) {
        
        lm1 <- lm(as.formula(paste0(mf_names[s], "~ species_pool")),
                  data = lm_df)
        
        s_cof[s] <- lm1$coefficients[2]
        
      }
      
      names(s_cof) <- c(mf_names)
      
      y[[i]] <- s_cof
      
    } 
    
    y <- bind_rows(y, .id = "run")
    
    y$number_functions <- d
    
    y
    
  }
  
  # loop over this function for all possible numbers of functions
  slopes_n_func <- vector("list", length = length(funcs))
  
  for(j in 2:length(funcs)) {
    
    slopes_n_func[[j]] <- 
      s_est(a = data, b = funcs, d = j)
    
  }
  
  bind_rows(slopes_n_func, .id = "n_func_rep")
  
}


# load the Lotka-Volterra simulated data
lv_mf_sims <- 
  read_csv(file = here("data/lv_mf_analysis_data.csv"))

# get a vector of function names
f_names <- names(select(lv_mf_sims, starts_with("F")))

# create a unique id column
lv_mf_sims$id <- paste(as.character(lv_mf_sims$run), as.character(lv_mf_sims$scenario), sep = ".")

# change the name of run to model
lv_mf_sims <- 
  lv_mf_sims %>%
  rename(model = run)

# split the data by this id column into a list
lv_mf_l <- 
  split(select(lv_mf_sims, -id), lv_mf_sims$id)


# apply the est_n_func to each of these unique replicates
mf_ests <- 
  lapply(lv_mf_l,
       
         function(x) {
           
         z <- est_n_func(data = x, funcs = f_names)
         
         z %>%
           mutate(scenario = x$scenario[1],
                  model = x$model[1])
         
       })


mf_ests_b <- bind_rows(mf_ests, .id = "id")


lv_mf_sims %>%
  pivot_longer(cols = all_of(f_names),
               names_to = "eco_function",
               values_to = "value") %>%
  ggplot(data = .,
       mapping = aes(x = species_pool, y = value, colour = id)) +
  geom_jitter(width = 0.3, alpha = 0.3) +
  geom_smooth(se = FALSE, method = "lm") +
  facet_wrap(~eco_function, scales = "free") +
  scale_colour_viridis_d() +
  theme_meta() +
  theme(legend.position = "none")

mf_ests_b %>%
  pivot_longer(cols = ends_with("mf"),
               names_to = "mf_metric",
               values_to = "mf_div_est") %>%
  ggplot(data = ., 
       mapping = aes(x = number_functions, y = mf_div_est, colour = id)) +
  geom_jitter(width = 0.3, alpha = 0.3) +
  geom_smooth(se = FALSE, method = "lm") +
  facet_wrap(~mf_metric, scales = "free") +
  scale_colour_viridis_d() +
  theme_meta() +
  theme(legend.position = "none")
















