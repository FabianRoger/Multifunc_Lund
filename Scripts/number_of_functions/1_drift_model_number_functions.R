
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Run the neutral model to test the number of functions question

# load relevant libraries
library(here)

# clear the current memory
rm(list = ls())

# link to scripts with the relevant functions
source(here("Scripts//ecological_drift_model.R"))
source(here("Scripts/MF_functions_collated.R"))

# how many model reps for each parameter combination?
n_reps <- 1000

# run the ecological drift model

# drift model parameters
p_change <- 0.025

drift.mod.list <- vector("list", length = n_reps)
for (i in 1:n_reps) {
  
  library(dplyr)
  
  drift.mod.list[[i]] <- 
  
    drift_model(lsp = c(2, 4, 6, 9),
                mono = "all",
                reps = 5,
                technical_reps = 2,
                rsp = 12,
                t_steps = 500,
                n0 = 500,
                prop_change = p_change,
                n_repeats = 1) %>%
    filter(time == max(time))
}

# bind the list into a data.frame
drift.mod.list <- bind_rows(drift.mod.list, .id = "model_run")
drift.mod.list$drift_parameter <- p_change

# reorganise the columns
drift.mod.list <- 
  drift.mod.list %>%
  select(drift_parameter, model_run, time, patch, local_species_pool, composition,
         species, abundance)

# write this model data.frame into a .csv file so the model does not have to be re-run
library(readr)

# make a folder to export the cleaned data
if(! dir.exists(here("data"))){
  dir.create(here("data"))
}

write_csv(x = drift.mod.list, file = here("data/drift_model_n_functions.csv"))

### END






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
  mod.out %>%
  group_by(drift_parameter, function_matrix, model_run, mod_id) %>% 
  mutate( across(.cols = starts_with("F_"), standardise) ) %>%
  mutate( abundance = standardise(abundance) ) %>%
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
            se = sd(estimate)/sqrt(n()),
            upp_ci = quantile(estimate, 0.975),
            low_ci = quantile(estimate, 0.025))

# plot a few relationships between species richness and abundance
# plot the abundance and SR relationship
names(mod.out)
fx.a <- 
  mod.out %>%
  filter(mod_id %in% sample(unique(mod.out$mod_id), 10)) %>%
  ggplot(data = .,
         mapping = aes(x = local_species_pool, y = abundance, colour = mod_id)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  theme_meta() +
  theme(legend.position = "none")

# plot the different models
fx.b <- 
  ggplot(data = sr_abun,
         mapping = aes(x = estimate)) +
  geom_density(alpha = 0.3, fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(sr_abun_ss$low_ci, sr_abun_ss$upp_ci), colour = "red") +
  xlab("abundance ~ SR est.") +
  theme_meta() +
  theme(legend.position = "none")
fx.b

# plot boxplots for each function
sr_funcs <- 
  bef_slopes %>%
  filter(response_var != "estimate_SR_abundance")

sr_funcs %>%
  group_by(response_var, function_matrix) %>%
  summarise(n = n())

f.ests <- unique(sr_funcs$response_var)

plots.fx3 <- vector("list", length = length(f.ests))
for (i in 1:length(f.ests)) {
  
  df.x <- 
    sr_funcs %>%
    filter(response_var == f.ests[i]) %>%
    ggplot(data = ., 
           mapping = aes(x = function_matrix, y = estimate, colour = function_matrix, fill = function_matrix)) +
    geom_jitter(width = 0.1, alpha = 0.05) +
    geom_boxplot(width = 0.1, outlier.shape = NA, colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_viridis_d() +
    scale_fill_viridis_d() +
    ylab("function ~ SR est.") +
    xlab("function matrix") +
    theme_meta() +
    theme(legend.position = "none")
  
  plots.fx3[[i]] <- df.x
  
}
names(plots.fx3) <- f.ests

sr_funcs %>%
  group_by(response_var, function_matrix) %>%
  summarise(mean_est = mean(estimate),
            median_est = median(estimate))

plots.fx3[[5]]


# combine these plots using patchwork
library(patchwork)

p.y <- fx.a/fx.b +
  plot_annotation(tag_levels = "a")

p.x <- 
  plots.fx3$estimate_SR_F_1 +  
  plots.fx3$estimate_SR_F_2 + 
  plots.fx3$estimate_SR_F_3 + 
  plots.fx3$estimate_SR_F_4 + 
  plots.fx3$estimate_SR_F_5 +
  plot_annotation(tag_levels = list(c("c", "d", "e", "f", "g") ) )

p.y - p.x + plot_layout(widths = c(1, 3.5))



### END

### Code for examples:

# to do this, we will need to run the models again

# output raw abundance data for plot from one random model
ran.runs <- sample(unique(drift.mod.list$drift_model_p_change_0.025$model_run), 1)

raw_abun <- 
  drift.mod.list$drift_model_p_change_0.025 %>%
  filter(model_run %in% ran.runs) %>%
  filter(abundance > 0)

View(raw_abun)

library(ggplot2)

cmp <- sample(unique(raw_abun$composition), 5)

raw_abun %>%
  filter(composition %in% cmp) %>%
  group_by(composition) %>%
  filter(patch == first(patch)) %>%
  ungroup() %>%
  filter(time %in% seq(1, 500, 10) ) %>%
  ggplot(data = ., 
         mapping = aes(x = time, y = abundance, colour = species)) +
  geom_line() +
  facet_wrap(~composition, scales = "free") +
  theme_classic() +
  theme(legend.position = "none")





