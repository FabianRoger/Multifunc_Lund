#'
#' @title Run the turnover approaches on the Jena and BIODEPTH data
#' 
#' @description Uses the AIC and SES - based turnover approaches to 
#' calculate the proportion of the species pool required to support different
#' ecosystem functions on the Jena data and two BIODEPTH datasets: Sweden
#' and Portugal
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)

# load relevant scripts
source("code/helper-plotting-theme.R")
source("code/helper-turnover-approach-functions.R")

# load the relevant datasets
jena_dat <- readRDS("data/jena_data.rds")
swe_dat <- readRDS("data/biodepth_swe_data.rds")
prt_dat <- readRDS("data/biodepth_prt_data.rds")

# get the species and function names for the different datasets

# Jena data

# get the species names
jena_cols <- names(jena_dat)
jena_sp <- jena_cols[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", jena_cols)]

# get the function names
jena_funcs <- c("biomass", "plantCN", "soilC", "soilorgC", "herbi", "micBMC", "phoact", "poll","rootBM")

# BIODEPTH: Sweden
swe_cols <- names(swe_dat)
swe_sp <- swe_cols[grepl(pattern = "[A-Z]{6}[1]", swe_cols)]

# get the function names
swe_funcs <- c("biomassY3", "root3", "cotton3", "wood3", "N.g.m2")
all(swe_funcs %in% swe_cols) 

# BIODEPTH: Portugal
prt_cols <- names(prt_dat)
prt_sp <- prt_cols[grepl(pattern = "[A-Z]{6}[1]", prt_cols)]

# get the function names
prt_funcs <- c("root3", "cotton3", "wood3", "N.Soil")
all(prt_funcs %in% prt_cols) 

# pull these datasets into a list
dat_list <- list(jena_dat, swe_dat, prt_dat)
dat_funcs <- list(jena_funcs, swe_funcs, prt_funcs)
dat_sp <- list(jena_sp, swe_sp, prt_sp)

# convert abundance data into presence absence
dat_list <- 
  mapply(function(D, S){
  D[,S] <- vegan::decostand(D[,S], method = "pa")
  return(D)
}, dat_list, dat_sp)

# set the number of randomisations to do for the null expectation
n_null <- 999

# run the AIC turnover approach on these datasets
aic_dat <- 
  
  mapply(function(D, Func, S) {
  
  # get the observed proportion of the species pool
  aic_obs <- 
    prop_species_pool(data = D, 
                      func_names = Func, 
                      sp_names = S, method = "AIC")
  
  # get the null proportion of the species pool
  aic_ran <- 
    prop_species_pool_random(data = D,
                             func_names = Func,
                             sp_names = S, method = "AIC", n = n_null)
  
  return( list(obs = aic_obs, null = aic_ran) )
  
}, dat_list, dat_funcs, dat_sp, SIMPLIFY = FALSE)

# set set the names of the datasets
names(aic_dat) <- c("jena", "swe", "prt")

# run the SES-based turnover approach
ses_dat <- 
  
  mapply(function(D, Func, S) {
    
    # get the observed proportion of the species pool
    ses_obs <- 
      prop_species_pool(data = D, 
                        func_names = Func, 
                        sp_names = S, method = "SES", n_ran = 999)
    
    # get the null proportion of the species pool
    ses_ran <- 
      prop_species_pool_random(data = D,
                               func_names = Func,
                               sp_names = S, method = "SES", n_ran = 999, 
                               n = n_null)
    
    return( list(obs = ses_obs, null = ses_ran) )
    
  }, dat_list, dat_funcs, dat_sp, SIMPLIFY = FALSE)

# set set the names of the datasets
names(ses_dat) <- c("jena", "swe", "prt")

# plot the results

# combine the AIC and SES lists
df_list <- c(aic_dat, ses_dat)

# set-up a list of titles
title <- c("AIC", " ", " ", "SES", " ", " ")
ylabs <- c("Prop. species pool", NA, NA, "Prop. species pool", NA, NA)
xlabs <- c(NA, NA, NA, "", "Number of functions", "")

# loop over the different datasets
plot_list <- vector("list", length = length(df_list))
for(i in 1:length(df_list)) {
  
  # get the observed data
  df_obs <- df_list[[i]][["obs"]]
  
  # get the number of functions
  n_funcs <- max(df_obs$num_funcs)
  
  # alter the factors of direction variable
  df_obs$eff_dir <- factor(df_obs$eff_dir)
  levels(df_obs$eff_dir) <- c("-ve effect", "+ve effect")
  
  # summarise the observed data into the mean line
  df_obs_sum <- 
    df_obs |>
    dplyr::group_by(num_funcs, eff_dir) |>
    dplyr::summarise(prop_sp_pool_m = mean(prop_sp_pool), .groups = "drop")
  
  # get the null distribution
  df_null <- df_list[[i]][["null"]]
  
  # alter the factors of direction variable
  df_null$eff_dir <- factor(df_null$eff_dir)
  levels(df_null$eff_dir) <- c("-ve effect", "+ve effect")
  
  # summarise the null into a 95% percentile
  df_null <- 
    df_null |>
    dplyr::group_by(num_funcs, eff_dir) |>
    summarise(PI_low = quantile(prop_sp_pool, 0.025),
              PI_high = quantile(prop_sp_pool, 0.975), .groups = "drop")
  
  p1 <- 
    ggplot() +
    geom_ribbon(data = df_null,
                mapping = aes(x = num_funcs, ymin = PI_low, ymax = PI_high, fill = eff_dir),
                alpha = 0.1, show.legend = FALSE, fill = "red") +
    geom_quasirandom(data = df_obs, 
                     mapping = aes(x = num_funcs, y = prop_sp_pool, colour = eff_dir),
                     width = 0.1, shape = 1, show.legend = FALSE,
                     alpha = 0.3, stroke = 0.25, colour = "red") +
    geom_line(data = df_obs_sum,
              mapping = aes(x = num_funcs, y = prop_sp_pool_m, colour = eff_dir),
              linewidth = 0.75, colour = "red") +
    scale_x_continuous(limits = c(0.9, n_funcs + 0.1), 
                       breaks = seq(1, n_funcs, length.out = if(ceiling(n_funcs/2) < 3){4} else{ceiling(n_funcs/2)} )) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ylab(if(is.na(ylabs[i])){NULL}else{ylabs[i]}) +
    xlab(if(is.na(xlabs[i])){NULL}else{xlabs[i]}) +
    facet_wrap(~eff_dir) +
    ggtitle(title[i]) + 
    theme_meta() +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          strip.text = element_text(size = 11),
          plot.title = element_text(vjust = -1))
  
  plot_list[[i]] <- p1
  
}

# combine into a single plot
f1 <- 
  cowplot::plot_grid(plotlist = plot_list, 
                     nrow = 2, ncol = 3,
                     rel_widths = c(1.4, 1.1, 1),
                     rel_heights = c(1, 1.075),
                     labels = letters[1:6],
                     label_size = 11, label_fontface = "plain")

ggplot2::ggsave(filename = "figures-paper-2/fig_1.svg", f1,
                units = "cm", width = 20, height = 14)

### END
