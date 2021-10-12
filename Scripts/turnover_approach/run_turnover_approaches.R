
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Run the turnover approach test on 1000 neutral models for three function matrices

# load relevant libraries
library(here)

# link to scripts with the relevant functions
source(here("Scripts/process_model_data.R"))
source(here("Scripts/turnover_approach/turnover_approach_functions.R"))
source(here("Scripts/ecological_drift_model.R"))
source(here("Scripts/function_plotting_theme.R"))

# run an example of an ecological drift model
drift_exp <- 
  drift_model(lsp = c(2, 4, 6, 9),
              mono = "all",
              reps = 5,
              technical_reps = 2,
              rsp = 12,
              t_steps = 500,
              n0 = 500,
              prop_change = 0.025,
              n_repeats = 1)

# load packages
library(dplyr)
library(ggplot2)

# plot the example of patches in the ecological drift model
drift_exp %>%
  filter(patch %in% sample(x = unique(d_mod$patch), size = 4 )) %>%
  filter(abundance > 0) %>%
  ggplot(data = .,
         mapping = aes(x = time, y = abundance, colour = species)) +
  geom_line() +
  facet_wrap(~patch, scales = "free") +
  theme_meta() +
  theme(legend.position = "none")

# load the readr library
library(readr)

# list of simulated data cluster (e.g. neutral model with 1000 runs with same parameters)
# this list is generated using the 1_drift_model_number_functions.R script (number_of_functions folder)
mod.list <- here("data/drift_model_n_functions.csv")

# load Laura's function matrices
f.list <- list.files(path = here("data"))
fm_names <- f.list[grepl(pattern = "SpeciesIDs", f.list)]

# load these data.frames into a list
library(readr)

func.list <- vector("list", length = length(fm_names))
for (i in 1:length(func.list)) {
  
  x <- read_csv(file = paste(here("data"), fm_names[i], sep = "/") )
  colnames(x)[1] <- "species"
  
  func.list[[i]] <- x
  
}


# for each dataset in mod.list and for each function matrix in func.list, run the AIC and SES turnover approaches

# set the number of reps for the SES method
ses.r <- 100

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
    func.reps[[j]] <- turnover_tester(model_dat = df.proc, function_matrix = func.list[[j]], ses_reps = ses.r) 
    
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

View(df.plot)

# df.plot$mod_id <- as.factor(df.plot$mod_id)
# levels(df.plot$mod_id) <- paste("m.", 1:length(unique(df.plot$mod_id)), sep = "")
# df.plot$mod_id <- as.character(df.plot$mod_id)

df.plot %>%
  group_by(mod_id, model_run, method, output_metric) %>%
  summarise(mean_metric = mean(metric, na.rm = TRUE), .groups = "drop") %>%
  filter(output_metric == "prop_incorrect") %>%
  ggplot(data = .,
         mapping = aes(x = method, y = mean_metric)) +
  geom_jitter() +
  facet_wrap(~mod_id)

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
