
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Prepare cluster of model replicates (i.e. 1000 simulations of one set of parameters)

# load relevant libraries
library(here)

# link to scripts with the relevant functions
source(here("Scripts/multispecies_lotka_volterra_model/process_model_data.R"))
source(here("Scripts/multispecies_lotka_volterra_model/per_capita_function_matrix_functions.R"))
source(here("Scripts/multispecies_lotka_volterra_model/turnover_approach_functions.R"))

# Laura will provide function matrices so that the identity effects in the DI-models are the same
# for now, I will use the func_matrix_generator function:

# choose the number of functions
n.f <- 5

# specialist
fm1 <- func_matrix_generator(species_list = unique(mod.list[[1]]$species),
                             func.n = n.f, func.spec = "specialist", prob.neg = 0.05)

# generalist
fm2 <- func_matrix_generator(species_list = unique(mod.list[[1]]$species),
                             func.n = n.f, func.spec = "generalist", prob.neg = 0.05)

func.list <- list(fm1, fm2)
func.list

# list of simulated data cluster (e.g. neutral model with 1000 runs with same parameters)
mod.list

# for each dataset in mod.list and for each function matrix in func.list, run the AIC and SES turnover approaches

mod.turnover.test <- vector("list", length = length(mod.list))
for (i in 1:length(mod.list)) {
  
  func.reps <- vector("list", length = length(func.list))
  for (j in 1:length(func.list)) {
    
    # process the simulated cluster using: process_sim_data
    df.proc <- 
      process_sim_data(model_data = mod.list[[i]], 
                       func.mat =  func.list[[j]], 
                       time_final = TRUE, 
                       species_abun = "pa")
    
    # run the different turnover approaches
    func.reps[[j]] <- turnover_tester(model_dat = df.proc, function_matrix = func.list[[j]], ses_reps = 1000) 
    
  }
  
  mod.turnover.test[[i]] <- bind_rows(func.reps, .id = "function_matrix")
  
}

# assign names to this object
names(mod.turnover.test) <- names(mod.list)

# bind this into one large data.frame
mod.turnover.df <- bind_rows(mod.turnover.test, .id = "parameter_combination")

# pull this into a longer data.frame
library(ggplot2)

df.plot <- 
  mod.turnover.df %>%
  pivot_longer(cols = starts_with("F_"),
               names_to = "function_id",
               values_to = "metric") %>%
  mutate(mod_id = paste(parameter_combination, function_matrix, sep = "_"), .before = 1)

df.plot$mod_id <- as.factor(df.plot$mod_id)
levels(df.plot$mod_id) <- paste("m.", 1:length(unique(df.plot$mod_id)), sep = "")
df.plot$mod_id <- as.character(df.plot$mod_id)

df.plot %>%
  ggplot(data = ., 
       mapping = aes(x = method, y = metric, colour = mod_id)) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic() +
  facet_wrap(~output_metric, scales = "free") +
  scale_colour_viridis_d() +
  theme(legend.position = "bottom",
        axis.text = element_text(colour = "black"))

### END
