
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Run the neutral model to test the number of functions question

# load relevant libraries
library(here)

# link to scripts with the relevant functions
source(here("Scripts/turnover_approach/ecological_drift_model.R"))


# how many model reps for each parameter combination?
n_reps <- 50

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


# process this model data

# load scripts with processing functions
source(here("Scripts/turnover_approach/process_model_data.R"))

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


# list of simulated data cluster (e.g. neutral model with 1000 runs with same parameters)
# this list is generated using the collate_model_data.R script
drift.mod.list

# for each dataset in mod.list and for each function matrix in func.list
# process the data

mod.out <- vector("list", length = length(drift.mod.list))
for (i in 1:length(drift.mod.list)) {
  
  func.reps <- vector("list", length = length(func.list))
  for (j in 1:length(func.list)) {
    
    # process the simulated cluster using: process_sim_data
    df.proc <- 
      process_sim_data(model_data = drift.mod.list[[i]], 
                       func.mat =  func.list[[j]], 
                       time_final = TRUE, 
                       species_abun = "raw")
    
    # add this to a list
    func.reps[[j]] <- df.proc
    
  }
  
  mod.out[[i]] <- bind_rows(func.reps, .id = "function_matrix")
  
}

# assign names to this object
names(mod.out) <- names(drift.mod.list)

# bind this into one large data.frame
mod.df <- bind_rows(mod.out, .id = "parameter_combination")
View(mod.df)
head(mod.df)

# add a variable for a unique identifier variable
mod.df <- 
  mod.df %>%
  mutate(mod_id = paste(parameter_combination, function_matrix, model_run, sep = ".")) %>%
  select(mod_id, parameter_combination:F_5)

length(unique(mod.df$mod_id))
head(mod.df)

# calculate the expected slope between local species pool diversity and abundance
# filter out first function matrix because the models for different function matrices are the same
library(purrr)
library(broom)

bef_slopes <- 
  mod.df %>%
  filter(function_matrix == first(function_matrix)) %>%
  group_by(parameter_combination, model_run) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(abundance ~ 1 + local_species_pool, data = .x) %>% 
                       tidy)) %>% 
  unnest(model) %>% 
  filter(term == "local_species_pool") %>%
  select(-data) %>%
  ungroup()


# make some preliminary figures
source(here("Scripts/function_plotting_theme.R"))

# plot first parameter combination
mod.df %>%
  filter(function_matrix == 1) %>%
  filter(parameter_combination == first(parameter_combination)) %>%
  filter(model_run %in% sample(unique(mod.df$model_run), 5)) %>%
  ggplot(data = .,
         mapping = aes(x = local_species_pool, y = abundance, colour = model_run)) +
  geom_jitter(width = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  theme_meta() +
  theme(legend.position = "none")

ggplot(data = filter(bef_slopes, parameter_combination == first(parameter_combination)) ,
       mapping = aes(x = estimate)) +
  geom_histogram() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("abundance ~ SR est.") +
  theme_meta() +
  theme(legend.position = "none")































