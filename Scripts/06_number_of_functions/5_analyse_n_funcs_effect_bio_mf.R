
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Analyse the drift model data for the number of functions

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
source(here("Scripts/linear_model_slope_function.R"))


# plot n-functions and multifunctionality

# read in the model data
sim.n.out <- read_csv(here("data/sim_n_functions.csv"))
View(head(sim.n.out))
dim(sim.n.out)

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
library(purrr)
library(broom)

nfunc_slopes <- 
  sim.n.out %>%
  group_by(drift_parameter, model_run, function_matrix, multifunctionality_metric) %>% 
  nest() %>% 
  mutate(nfunc.bef_slope = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "diversity_mf_est"))  ) %>%
  unnest(starts_with("nfunc"))  %>% 
  select(-data) %>%
  ungroup()

View(nfunc_slopes)

plots.fx2 <- vector("list", length = length(mf.metric.list))
for(i in 1:length(mf.metric.list)) {
  
  df <- 
    nfunc_slopes %>%
    filter(multifunctionality_metric == mf.metric.list[i]) %>%
    mutate(function_matrix = as.character(function_matrix))
  
  df_sum <- 
    df %>%
    group_by(function_matrix) %>%
    summarise(mean_estimate = mean(number_of_functions_diversity_mf_est, na.rm = TRUE))
  
  df.x <- 
    ggplot(data = df,
           mapping = aes(x = number_of_functions_diversity_mf_est, fill = function_matrix)) +
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


# explore these patterns in more detail

# how does variation among functions change as more functions are considered?
View(sim.n.out)
names(sim.n.out)

sim.n.sub <- 
  sim.n.out %>%
  select(model_run, function_matrix, mod_id, func.comb.id, n.func.id, number_of_functions, sd_funcs, cv_funcs, range_funcs, sd_cor, cv_cor, range_cor, mean_cor) %>%
  distinct()
head(sim.n.sub)
names(sim.n.sub)

sim.n.slopes <- 
  sim.n.sub %>%
  group_by(mod_id, model_run, function_matrix) %>% 
  nest() %>% 
  mutate(nfunc.cor_sd = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "sd_cor")),
         nfunc.cor_mean = map(data, ~lm.cleaner(data = .x, explanatory = "number_of_functions", response = "mean_cor")) ) %>%
  unnest(starts_with("nfunc"))  %>% 
  select(-data) %>%
  ungroup()

ggplot(data = sim.n.slopes,
       mapping = aes(x = number_of_functions_sd_cor)) +
  geom_histogram(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~function_matrix, scales = "free_y") +
  theme_meta()

ggplot(data = sim.n.slopes,
       mapping = aes(x = number_of_functions_mean_cor)) +
  geom_histogram(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~function_matrix, scales = "free_y") +
  theme_meta()


