
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code to pull model data into a list

# load relevant libraries
library(here)

# link to scripts with the relevant functions
source(here("Scripts/multispecies_lotka_volterra_model/ecological_drift_model.R"))
source(here("Scripts/multispecies_lotka_volterra_model/stachova_leps_2010_LK_model_function.R"))

# how many model reps for each parameter combination?
n_reps <- 5

# run the ecological drift model

# drift model parameters
p_change <- c(0.025, 0.05, 0.1)

drift.mod.list <- vector("list", length = length(p_change))
for (i in 1:length(p_change)) {
  
  drift.mod.list[[i]] <- 
    
    drift_model(lsp = c(2, 4, 6, 9),
                mono = "all",
                reps = 5,
                technical_reps = 2,
                rsp = 12,
                t_steps = 500,
                n0 = 500,
                prop_change = p_change[i],
                n_repeats = n_reps
                )
}

names(drift.mod.list) <- paste("drift_model_p_change_", p_change, sep = "" )


# run the lotka-volterra models under selection and complementarity assumptions

# lotka-volterra model parameters
alpha_mean <- c(0.7, 0.8, 0.9, 0.1, 0.2, 0.3)
sel_comp <- rep(c("selection", "complementarity"), each = length(alpha_mean)/2)

lk.mod.list <- vector("list", length = length(alpha_mean))
for (i in 1:length(alpha_mean)) {
  
  x <- 
    s_l_function(lsp = c(2, 4, 6, 9),
                 mono = "all",
                 reps = 5,
                 technical_reps = 2,
                 rsp = 12,
                 t_steps = 500,
                 n0 = 500,
                 a_mean = alpha_mean[i], a_sd = 0.2, a_min = -0.1,
                 a_max = 1.2, a_spp = 1, sim.comp = "sym",
                 k_min = 200, k_max = 800,
                 r_min = 0.1, r_max = 0.5,
                 NBE = sel_comp[i], random.seed = 1467,
                 n_repeats = n_reps)
  
  lk.mod.list[[i]] <- x$data.raw
  
}

names(lk.mod.list) <- paste("lk_model_alpha_mean_", alpha_mean, "_", sel_comp,  sep = "" )


# combine the drift model and lotka-volterra model outputs into a list
mod.list <- c(drift.mod.list, lk.mod.list)

# output this as a .rds file at some point

### END
