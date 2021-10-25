
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Analyse the drift model data for the number of functions

# Next steps:

# Figure out why these different metrics change with the number of functions considered

# load relevant libraries
library(here)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

# clear the current memory
rm(list = ls())

# get scripts to call functions from
source(here("Scripts/function_plotting_theme.R"))
source(here("Scripts/MF_functions_collated.R"))

### plot example of the ecological drift model
source(here("Scripts/ecological_drift_model.R"))

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
  scale_colour_viridis_d() +
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

# define function to efficiently output the slope
lm.cleaner <- function(data, response, explanatory, output_prefix = "x") {
  
  x <- 
    lm(reformulate(explanatory, response), data = data) %>% 
    tidy %>% 
    filter(term == explanatory) %>% 
    select(!!paste(output_prefix, response, sep = "") := estimate )
  
  return(x)
  
}

raw_slopes <- 
  drift_mod_dat %>%
  group_by(function_matrix, mod_id) %>% 
  mutate(across(.cols = c("abundance", "F_1", "F_2", "F_3", "F_4", "F_5"), standardise)) %>%
  nest() %>% 
  summarise(total_abun_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "abundance")),
         F1_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_1")),
         F2_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_2")),
         F3_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_3")),
         F4_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_4")),
         F5_slope = map(data, ~lm.cleaner(data = .x, explanatory = "local_species_pool", response = "F_5")) ) %>%
  unnest(ends_with("slope"))

# input data on total abundance for each function matrix is the same
plot_df <- 
  raw_slopes %>%
  filter(function_matrix == 1)

x1 <- 
  ggplot(data = plot_df,
       mapping = aes(x = xabundance)) +
  geom_histogram(colour = "transparent", alpha = 0.5) +
  geom_vline(xintercept = mean(plot_df$xabundance), colour = "red") +
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  ylab("Count") +
  xlab("Est. (richness ~ abundance)") +
  theme_meta()

# plot the different individual function distributions for each function matrix
f1_5_plots <- 
  raw_slopes %>%
  select(-xabundance)

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
  geom_histogram(position = "identity", alpha = 0.4) + 
  geom_vline(data = f1_5_plots_sum, mapping = aes(xintercept = slope, colour = variable)) +
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  ylab(NULL) +
  xlab("Est. (richness ~ function)") +
  facet_wrap(~function_matrix) +
  theme_meta() +
  theme(legend.title = element_blank())

# combine these plots using patchwork
library(patchwork)

x12 <- 
  x1 + x2 + plot_layout(ncol = 2, widths = c(1, 1.6)) +
  plot_annotation(tag_levels = "a")

ggsave(filename = here("Figures/fig_x12.png"), plot = x12, width = 21, height = 10,
       units = "cm")

### plot n-functions and multifunctionality

# read in the model data
sim.n.out <- read_csv(here("data/sim_n_functions.csv"))
View(head(sim.n.out))

# plot a few illustrative examples for each multifunctionality metric
mf.metric.list <- unique(sim.n.out$multifunctionality_metric)
mf.names <- c("sum MF", "ave. MF", "Pasari MF", "thresh 30 MF", "thresh 70 MF")

dat.fx1 <- 
  sim.n.out %>%
  filter(mod_id %in% sample(unique(sim.n.out$mod_id), 10 )) %>%
  mutate(function_matrix = as.character(function_matrix))

names(dat.fx1 )

plots.fx1 <- vector("list", length = length(mf.metric.list))
for(i in 1:length(mf.metric.list)) {
  
  df.x <- 
    dat.fx1 %>%
    filter(multifunctionality_metric == mf.metric.list[[i]]) %>%
    ggplot(data = .,
           mapping = aes(x = number_of_functions, y = diversity_mf_est, # change to diversity if running again
                         group = mod_id, colour = function_matrix)) +
    geom_jitter(width = 0.1, alpha = 0.2) +
    geom_smooth(method = "lm", se = FALSE, width = 0.5) +
    scale_colour_viridis_d() +
    ggtitle(mf.names[i]) +
    xlab("Number of functions") +
    ylab("MF BEF-slope") +
    theme_meta() +
    theme(legend.position = "none",
          plot.title = element_text(family = "sans", size = 10, margin=margin(-10,0,0,0)))
  
  plots.fx1[[i]] <- df.x
  
}
names(plots.fx1) <- mf.metric.list
plots.fx1$Pasari_MF


# plot the expectations from all simulations for each multifunctionality metric
head(sim.n.out)

# calculate the slope between number of functions and multifuncitonal BEF slope 
# and other summary statistics for each multifunctionality metric
nfunc_slopes <- 
  sim.n.out %>%
  group_by(drift_parameter, model_run, function_matrix, multifunctionality_metric) %>% 
  nest() %>% 
  mutate(nfunc.bef_slope = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "diversity_mf_est")),
         nfunc.cv = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "cv_MF")),
         nfunc.range = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "range_MF")),
         nfunc.min = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "min_MF")),
         nfunc.max = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "max_MF"))  ) %>%
  unnest(starts_with("nfunc"))  %>% 
  select(-data) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("x"),
               names_to = "response_var",
               values_to = "estimate")

View(nfunc_slopes)

nfunc_slopes$response_var %>% unique()

# get the range of multifunctional BEF slopes
nfunc_slopes %>%
  filter(response_var == "xdiversity_mf_est") %>%
  pull(estimate) %>%
  range()

plots.fx2 <- vector("list", length = length(mf.metric.list))
for(i in 1:length(mf.metric.list)) {
  
  df <- 
    nfunc_slopes %>%
    filter(response_var == "xdiversity_mf_est",
           multifunctionality_metric == mf.metric.list[i]) %>%
    mutate(function_matrix = as.character(function_matrix))
  
  df_sum <- 
    df %>%
    group_by(function_matrix) %>%
    summarise(mean_estimate = mean(estimate, na.rm = TRUE))
  
  df.x <- 
    ggplot(data = df,
           mapping = aes(x = estimate, fill = function_matrix)) +
    geom_histogram(alpha = 0.5, colour = "transparent") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(data = df_sum, mapping = aes(xintercept = mean_estimate),
               colour = "red") +
    scale_fill_viridis_d() +
    scale_x_continuous(limits = c(-0.11, 0.11)) +
    ggtitle("") +
    xlab("Estimate") +
    scale_colour_viridis_d() +
    facet_wrap(~function_matrix, scales = "free") +
    theme_meta() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, size = 10))
    
  
  plots.fx2[[i]] <- df.x
  
}
names(plots.fx2) <- mf.metric.list

# link these plots using patchwork
library(patchwork)

f.x1 <- 
  plots.fx1$ave._MF + plots.fx2$ave._MF + 
  plots.fx1$sum_MF + plots.fx2$sum_MF + 
  plots.fx1$Pasari_MF + plots.fx2$Pasari_MF + 
  plots.fx1$thresh.30_MF + plots.fx2$thresh.30_MF +
  plots.fx1$thresh.70_MF + plots.fx2$thresh.70_MF + 
  plot_layout(ncol = 2, widths = c(1, 2.5)) +
  plot_annotation(tag_levels = "a")

ggsave(filename = here("Figures/fig_x1.png"), plot = f.x1, width = 20, height = 27,
       units = "cm")

# what about among-function variation?
names(sim.n.out)
sim.n.out %>%
  select(parameter_combination, model_run, function_matrix, number_of_functions, sd_funcs) %>%
  distinct() %>%
  group_by(parameter_combination, model_run, function_matrix) %>% 
  nest() %>% 
  mutate(nfunc.among_sd = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "sd_funcs")) ) %>%
  unnest(nfunc.among_sd)  %>% 
  select(-data) %>%
  ungroup() %>%
  ggplot(data = .,
         mapping = aes(x = estimate_n_func_sd_funcs)) +
  geom_histogram(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~function_matrix, scales = "free_y") +
  theme_meta()





