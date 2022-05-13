
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Plot a general overview of the behaviour of the ecological drift model

# load relevant libraries
library(here)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

# clear the current memory
rm(list = ls())

# get scripts to call functions from
source(here("Scripts/01_general_functions/function_plotting_theme.R"))
source(here("Scripts/01_general_functions/MF_functions_collated.R"))
source(here("Scripts/01_general_functions/linear_model_slope_function.R"))

# plot example of the ecological drift model
source(here("Scripts/01_general_functions/ecological_drift_model.R"))

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

# plot the example of patches in the ecological drift model
p1 <- 
  drift_exp %>%
  filter(patch %in% sample(x = unique(drift_exp$patch), size = 4 )) %>%
  filter(abundance > 0) %>%
  ggplot(data = .,
         mapping = aes(x = time, y = abundance, colour = species)) +
  geom_line() +
  ylab("Total abundance") +
  xlab("Time (generations)") +
  facet_wrap(~patch, scales = "free") +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  theme_meta() +
  theme(legend.position = "none")

ggsave(filename = here("Figures/fig_y1.png"), plot = p1, width = 20, height = 12,
       units = "cm")

# load the raw drift model data
drift_mod_dat <- read_csv(here("data/drift_model_n_functions_processed.csv"))
head(drift_mod_dat)

# calculate the slope between number of functions and multifuncitonal BEF slope 
# and other summary statistics for each multifunctionality metric
library(purrr)
library(broom)

raw_slopes <- 
  drift_mod_dat %>%
  group_by(function_matrix, mod_id) %>% 
  mutate(across(.cols = c("abundance", "F_1", "F_2", "F_3", "F_4", "F_5"), ~standardise_functions(x = ., method = "max") )) %>%
  nest() %>% 
  summarise(total_abun_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "abundance")),
         F1_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_1")),
         F2_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_2")),
         F3_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_3")),
         F4_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_4")),
         F5_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_5")) ) %>%
  unnest(ends_with("slope"))


# input data on total abundance for each function matrix is the same
# therefore, we can plot the abundance-richness correlation using just one function matrix
unique(raw_slopes$function_matrix)
plot_df <- 
  raw_slopes %>%
  filter(function_matrix == 1)
head(plot_df)
names(plot_df)

x1 <- 
  ggplot(data = plot_df,
       mapping = aes(x = local_species_pool_abundance)) +
  # geom_density(fill = "grey") +
  geom_histogram(colour = "transparent", alpha = 0.5) +
  geom_vline(xintercept = mean(plot_df$local_species_pool_abundance), colour = "red") +
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  ylab("Count") +
  xlab("Est. (richness ~ abundance)") +
  theme_meta()
x1

# plot the different individual function distributions for each function matrix
f1_5_plots <- 
  raw_slopes %>%
  dplyr::select(-local_species_pool_abundance)

names(f1_5_plots) <- c("mod_id", "function_matrix", "F1", "F2", "F3", "F4", "F5")

f1_5_plots <- 
  f1_5_plots %>%
  pivot_longer(cols = starts_with("F", ignore.case = FALSE),
               names_to = "variable",
               values_to = "slope")

f1_5_plots_sum <- 
  f1_5_plots %>%
  group_by(function_matrix, variable) %>%
  summarise(slope = mean(slope, na.rm = TRUE),
            n = n())

x2 <- 
  ggplot(data = f1_5_plots,
         mapping = aes(x = slope, fill = variable)) +
  geom_histogram(position = "identity", alpha = 0.2) + 
  geom_vline(data = f1_5_plots_sum, mapping = aes(xintercept = slope, colour = variable)) +
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  scale_fill_viridis_d(option = "C", end = 0.9) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  ylab(NULL) +
  xlab("Est. (richness ~ function)") +
  facet_wrap(~function_matrix) +
  theme_meta() +
  theme(legend.title = element_blank())
x2

# combine these plots using patchwork
library(patchwork)

x12 <- 
  x1 + x2 + plot_layout(ncol = 2, widths = c(1, 1.6)) +
  plot_annotation(tag_levels = "a")
x12

ggsave(filename = here("Figures/fig_x12.png"), plot = x12, width = 21, height = 10,
       units = "cm")

### END
