#'
#' @title Plot the results of the turnover approach analysis Jena data
#' 
#' @description Load the output from the AIC-based and SES-based turnover
#' analyses.
#'

# load relevant scripts
source("code/helper-plotting-theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)
library(cowplot) 
library(ggbeeswarm)

# read the output files
aic_dat <- readRDS("code/02-paper-1/AIC_output.rds")
ses_dat <- readRDS("code/02-paper-1/SES_output.rds")

# combine the AIC and SES lists
df_list <- c(aic_dat, ses_dat)

# set-up a list of titles
title <- c("AIC - approach", "SES - approach")
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
              PI_high = quantile(prop_sp_pool, 0.975), 
              prop_sp_pool_m = mean(prop_sp_pool), .groups = "drop")
  
  p1 <- 
    ggplot() +
    geom_ribbon(data = df_null,
                mapping = aes(x = num_funcs, ymin = PI_low, ymax = PI_high, fill = eff_dir),
                alpha = 0.25, show.legend = FALSE, fill = "grey") +
    geom_line(data = df_null,
              mapping = aes(x = num_funcs, y = prop_sp_pool_m, colour = eff_dir),
              linewidth = 0.5, colour = "grey", linetype = "dashed", alpha = 0.8) +
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

# check a plot
plot_list[[1]]

# combine into a single plot
f1 <- 
  cowplot::plot_grid(plotlist = plot_list, 
                     nrow = 1, ncol = 2,
                     labels = letters[1:2],
                     label_size = 11, label_fontface = "plain")

ggplot2::ggsave(filename = "figures-paper-1/fig_4.svg", f1,
                units = "cm", width = 14, height = 7)




