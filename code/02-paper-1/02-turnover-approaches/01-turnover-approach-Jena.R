#'
#' @title Run the turnover approaches on the Jena data
#' 
#' @description Uses the AIC and SES - based turnover approaches to 
#' calculate the proportion of the species pool required to support different
#' ecosystem functions on the Jena data
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

# get the species and function names for the different datasets

# Jena data

# get the species names
jena_cols <- names(jena_dat)
jena_sp <- jena_cols[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", jena_cols)]

# get the function names
jena_funcs <- c("biomass", "plantCN", "soilC", "soilorgC", "herbi", "micBMC", "phoact", "poll","rootBM")

# convert abundance data into presence absence
jena_dat[,jena_sp] <- vegan::decostand(jena_dat[,jena_sp], method = "pa")


# set the number of randomisations to do for the null expectation
n_null <- 10

# get the observed proportion of the species pool
aic_obs <- 
  prop_species_pool(data = jena_dat, 
                    func_names = jena_funcs, 
                    sp_names = jena_sp, method = "AIC")

# get the null proportion of the species pool
aic_ran <- 
  prop_species_pool_random(data = jena_dat, 
                           func_names = jena_funcs, 
                           sp_names = jena_sp, method = "AIC", n = n_null)

# combine into a list
aic_dat <- list( list(obs = aic_obs, null = aic_ran) )

# run the SES-based turnover approach

# get the observed proportion of the species pool
ses_obs <- 
  prop_species_pool(data = jena_dat, 
                    func_names = jena_funcs, 
                    sp_names = jena_sp, method = "SES", n_ran = 999)

# get the null proportion of the species pool
ses_ran <- 
  prop_species_pool_random(data = jena_dat, 
                           func_names = jena_funcs, 
                           sp_names = jena_sp, method = "SES", n_ran = 999, 
                           n = n_null)

# combine into a list
ses_dat <- list( list(obs = ses_obs, null = ses_ran) )

# plot the results

# combine the AIC and SES lists
df_list <- c(aic_dat, ses_dat)

# set-up a list of titles
title <- c("AIC", "SES")
ylabs <- c("Prop. species pool", NA)
xlabs <- c("Number of functions", "Number of functions")

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
p <- 
  cowplot::plot_grid(plotlist = plot_list, 
                     nrow = 1, ncol = 2,
                     rel_widths = c(1.075, 1),
                     rel_heights = c(1, 1),
                     labels = c("a", "b"),
                     label_size = 11, label_fontface = "plain")
plot(p)

ggplot2::ggsave(filename = "figures-paper-1/fig_4.svg", p,
                units = "cm", width = 18, height = 8)

### END
