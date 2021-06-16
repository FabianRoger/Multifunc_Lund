
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

lm.cleaner <- function(data, response, explanatory) {
  
  x <- 
    lm(reformulate(explanatory, response), data = data) %>% 
    tidy %>% 
    filter(term == explanatory) %>% 
    select(!!paste("estimate_SR_", response, sep = "") := estimate )
  
  return(x)
  
  }

bef_slopes <- 
  mod.df %>%
  group_by(parameter_combination, model_run, function_matrix) %>% 
  nest() %>% 
  mutate(model_abun = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "abundance")),
         model_F1 = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_1")),
         model_F2 = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_2")),
         model_F3 = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_3")),
         model_F4 = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_4")),
         model_F5 = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_5")) ) %>%
  unnest(c("model_abun", paste("model_F", 1:5, sep = "" ))) %>% 
  select(-data) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("estimate"),
               names_to = "response_var",
               values_to = "estimate")
View(bef_slopes)


# make some preliminary figures
source(here("Scripts/function_plotting_theme.R"))
library(ggplot2)

# get the species richness ~ abundance slopes (do not depend on function matrix)
sr_abun <- 
  bef_slopes %>%
  filter(function_matrix == first(function_matrix)) %>%
  filter(response_var == "estimate_SR_abundance") %>%
  select(-function_matrix)

# calculate summary statistics
sr_abun_ss <- 
  sr_abun %>%
  summarise(n = n(), 
            mean = mean(estimate),
            se = sd(estimate)/sqrt(n())) %>%
  mutate(upp_ci = mean + qt(p = 0.975, df = n)*se,
         low_ci = mean - qt(p = 0.975, df = n)*se)

# plot a few relationships between species richness and abundance
# plot the abundance and SR relationship
names(mod.df)
mod.df %>%
  filter(mod_id %in% sample(unique(mod.df$mod_id), 5)) %>%
  ggplot(data = .,
         mapping = aes(x = local_species_pool, y = abundance, colour = mod_id)) +
  geom_jitter(width = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  theme_meta() +
  theme(legend.position = "none")

# plot the different models
ggplot(data = sr_abun,
         mapping = aes(x = estimate)) +
  geom_histogram(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(sr_abun_ss$low_ci, sr_abun_ss$upp_ci), colour = "red") +
  xlab("abundance ~ SR est.") +
  theme_meta() +
  theme(legend.position = "none")


# plot boxplots for each function
sr_funcs <- 
  bef_slopes %>%
  filter(response_var != "estimate_SR_abundance")

sr_funcs %>%
  group_by(response_var, function_matrix) %>%
  summarise(n = n())

ggplot(data = sr_funcs, 
       mapping = aes(x = function_matrix, y = estimate)) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~response_var) +
  ylab("function ~ SR est.") +
  xlab("function matrix") +
  theme_meta()

