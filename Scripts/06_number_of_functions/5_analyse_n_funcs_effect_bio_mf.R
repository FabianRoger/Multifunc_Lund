
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
source(here("Scripts/01_general_functions/function_plotting_theme.R"))
source(here("Scripts/01_general_functions/MF_functions_collated.R"))
source(here("Scripts/01_general_functions/linear_model_slope_function.R"))


# plot n-functions and multifunctionality

# read in the model data
sim.n.out <- read_csv(here("data/sim_n_functions.csv"))
View(head(sim.n.out))
dim(sim.n.out)
unique(sim.n.out$function_matrix)

# plot a few illustrative examples for each multifunctionality metric
mf.metric.list <- unique(sim.n.out$multifunctionality_metric)
mf.names <- c("sum MF", "ave. MF", "sd MF", "Pasari MF", "ENF_MF", "thresh 30 MF", "thresh 70 MF")

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
    geom_smooth(method = "lm", se = FALSE, size = 0.5) +
    scale_colour_viridis_d(option = "C", end = 0.9) +
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
  dplyr::select(-data) %>%
  ungroup()
View(nfunc_slopes)

nfunc_slopes %>%
  group_by(multifunctionality_metric, function_matrix) %>%
  summarise(mean_slope = mean(number_of_functions_diversity_mf_est) )

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
    geom_histogram(alpha = 0.5, colour = "transparent", bins = 30) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(data = df_sum, mapping = aes(xintercept = mean_estimate),
               colour = "red") +
    scale_fill_viridis_d(option = "C", end = 0.9) +
    scale_x_continuous(limits = c(-0.11, 0.11)) +
    ggtitle("") +
    xlab("Estimate") +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    facet_wrap(~function_matrix, scales = "free") +
    theme_meta() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, size = 10))
  
  
  plots.fx2[[i]] <- df.x
  
}
names(plots.fx2) <- mf.metric.list

names(plots.fx2)
plots.fx2$sum_MF

# link these plots using patchwork
library(patchwork)

f.x1 <- 
  plots.fx1$ave_MF + plots.fx2$ave_MF + 
  plots.fx1$sum_MF + plots.fx2$sum_MF + 
  # plots.fx1$sd_MF + plots.fx2$sd_MF + 
  plots.fx1$Pasari_MF + plots.fx2$Pasari_MF + 
  # plots.fx1$ENF_MF + plots.fx2$ENF_MF + 
  plots.fx1$thresh_30_MF + plots.fx2$thresh_30_MF +
  plots.fx1$thresh_70_MF + plots.fx2$thresh_70_MF + 
  plot_layout(ncol = 2, widths = c(1, 2.5)) +
  plot_annotation(tag_levels = "a")
f.x1

ggsave(filename = here("Figures/fig_x1.png"), plot = f.x1, width = 18, height = 35,
       units = "cm")


# explore these patterns in more detail

# how does variation among functions change as more functions are considered?
View(sim.n.out)
names(sim.n.out)

sim.n.sub <- 
  sim.n.out %>%
  dplyr::select(model_run, function_matrix, mod_id, func.comb.id, n.func.id, number_of_functions, sd_funcs, cv_funcs, range_funcs, sd_cor, cv_cor, range_cor, mean_cor) %>%
  distinct()
head(sim.n.sub)
names(sim.n.sub)
View(sim.n.sub)

# analyse examples from the raw model data
source(here("Scripts/01_general_functions/process_model_data.R"))

# load the raw data
mod.out <- read_csv(here("data/drift_model_n_functions_processed.csv"))
head(mod.out)

# get one random model from function matrix 2
mod1 <- 
  mod.out %>%
  filter(function_matrix == 2)

mod1 <- mod1[mod1$mod_id == sample(mod1$mod_id, 1), ]
head(mod1)

# write the adf.data into a new data.frame
dat.in <- mod1

# use the function to get a list of all function combinations
list.func.names <- get.function.combinations(function.names = c("F_1", "F_2", "F_3", "F_4", "F_5"))

# set an output list
list.out <- vector("list", length(list.func.names))

# loop over each set of function names
for (i in 1:length(list.func.names)) {
  
  # get the vector of function names
  sample.func.names <- list.func.names[[i]]
  
  # calculate the multifunctionality metrics
  data.mf <- 
    dat.in %>%
    mutate(sum_MF = MF_sum(adf = dat.in, vars = sample.func.names),
           ave_MF = MF_av(adf = dat.in, vars = sample.func.names, stand_method = "z_score_abs"),
           sd_MF = MF_sd(adf = dat.in, vars = sample.func.names, stand_method = "z_score_abs"),
           Pasari_MF = MF_pasari(adf = dat.in, vars = sample.func.names),
           ENF_MF = hill_multifunc(adf = dat.in, vars = sample.func.names, scale = 1, HILL = TRUE),
           thresh_30_MF = single_threshold_mf(adf = dat.in, vars = sample.func.names, thresh = 0.3),
           thresh_70_MF  = single_threshold_mf(adf = dat.in, vars = sample.func.names, thresh = 0.7) 
    ) %>%
    mutate(n_func = length(sample.func.names),
           func_id = paste(sample.func.names, collapse = ""))
  
  
  data.mf$sd_func <- apply(dat.in[, sample.func.names], 1, sd)
    
  list.out[[i]] <- data.mf
  
  }

bind_rows(list.out) %>%
  select(n_func, func_id, local_species_pool, sd_func) %>%
  distinct() %>%
  mutate(n_func = as.character(n_func)) %>%
  ggplot(data = .,
         mapping = aes(x = local_species_pool, y = sd_func, colour = n_func)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_meta()

bind_rows(list.out) %>%
  pivot_longer(cols = ends_with("_MF"),
               names_to = "metric",
               values_to = "MF") %>%
  ggplot(data = .,
         mapping = aes(x = local_species_pool, y = MF, 
                       group = func_id, colour = as.character(n_func))) +
  geom_jitter(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~metric, scales = "free") +
  theme_meta() +
  theme(legend.position = "bottom")


# Does the variation decline with species richness more with more functions?

# Yes, with function matrix 2 i.e. different functional contributions of species
ggplot(data = sim.n.sub,
       mapping = aes(x = number_of_functions, y = sd_cor) ) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  facet_wrap(~function_matrix) +
  ylab("Spearman's rho(species richness ~ SD)") +
  xlab("Number of functions") +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  theme_meta()

ggsave(filename = here("Figures/fig_S2.png"), width = 16, height = 9,
       units = "cm")

# does the range in multifunctionality change with more functions?
ggplot(data = sim.n.out %>%
         mutate(function_matrix = as.character(function_matrix)),
       mapping = aes(x = number_of_functions, y = range_MF, colour = function_matrix) ) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  facet_wrap(~multifunctionality_metric, scales = "free") +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Number of functions") +
  ylab("Range of multifunctionality (max - min)") +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  guides(colour = guide_legend(title="Function matrix")) +
  geom_hline(yintercept = 0, colour = "black", linetype = "dashed") +
  theme_meta() +
  theme(legend.position = "top", 
        legend.key=element_blank(),
        legend.title = element_text(size = 11))

ggsave(filename = here("Figures/fig_S3.png"), width = 18, height = 20,
       units = "cm")

ggplot(data = sim.n.out %>%
         mutate(function_matrix = as.character(function_matrix)),
       mapping = aes(x = number_of_functions, y = sd_MF, colour = function_matrix) ) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  facet_wrap(~multifunctionality_metric, scales = "free") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, colour = "red") +
  theme_meta()

ggplot(data = sim.n.out %>%
         mutate(function_matrix = as.character(function_matrix)),
       mapping = aes(x = number_of_functions, y = min_MF, colour = function_matrix) ) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  facet_wrap(~multifunctionality_metric, scales = "free") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, colour = "red") +
  theme_meta()

sim.n.out %>%
  filter(range_MF == 0) %>%
  View()

sim.n.out %>%
  filter(multifunctionality_metric == "sd_MF") %>%
  View()

### END
