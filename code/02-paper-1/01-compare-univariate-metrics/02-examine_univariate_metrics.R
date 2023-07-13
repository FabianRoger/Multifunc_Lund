#'
#' @title Compare univariate metrics
#' 
#' @description Compares univariate metrics for almost all possible combinations
#' of function distributions
#'

# load relevant scripts
source("code/helper-univariate-mf-functions.R")
source("code/helper-plotting-theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

# generate a data set with many different functional distributions
df <- 
  tidyr::crossing(F1 = seq(0, 1, 0.1),
                  F2 = seq(0, 1, 0.1),
                  F3 = seq(0, 1, 0.1),
                  F4 = seq(0, 1, 0.1))

# get a vector of function names
f_names <- names(df)

# calculate the mean among functions
df$MF_ave <- apply(df[,names(df) %in% f_names], 1, mean )

# calculate inverse simpson as a measure of evenness

# get max inverse simpson for these data
max_inv_simp <- vegan::diversity(x = c(1, 1, 1, 1), index = "invsimp")
print(max_inv_simp)

df$MF_inv_simp <- 
  apply(df[,names(df) %in% f_names], 1, function(x) { 
  
    y <- vegan::diversity(x = x, index = "invsimp") 
    z <- ifelse(is.infinite(y), 1, y)

    return(z)  
    
    } )

# use the distinct function to avoid duplicates
df <- dplyr::distinct(df)

# plot histograms of a few representative points
rep_points <- 
  dplyr::tibble(MF_ave = c(0.5, 0.5, 0.1, 0.75),
                MF_inv_simp = c(2, 4, 3, 3.2)
                )

# get points that match these criteria
rep_points_list <- vector("list", length = nrow(rep_points))
for(i in 1:nrow(rep_points)) {
  
  x <- 
    df |>
    dplyr::filter( dplyr::near(MF_ave, rep_points$MF_ave[i], 0.025 ),
                   dplyr::near(MF_inv_simp, rep_points$MF_inv_simp[i], 0.025 )) |>
    dplyr::sample_n(size = 1)
  
 rep_points_list[[i]] <- x 
  
}

# get the two extreme points
ex_p1 <- list( dplyr::filter(df, MF_ave == 1, MF_inv_simp == 4) )
ex_p2 <- list( dplyr::filter(df, MF_ave == 0, MF_inv_simp == 1) )
ex_points_df <- dplyr::bind_rows(ex_p1, ex_p2, .id = "id")

# bind this into a data.frame
rep_points_df <- dplyr::bind_rows(rep_points_list, .id = "id")

col_pal <- wesanderson::wes_palette(name = "Darjeeling1", n = 6, type = "continuous")

# plot the data
p1 <- 
  ggplot() +
  geom_jitter(data = df,
             mapping = aes(x = MF_ave, y = MF_inv_simp), 
             size = 1.5, alpha = 0.25, stroke = 0.1, shape = 1, width = 0.005) +
  geom_point(data = rep_points_df, 
             mapping = aes(x = MF_ave, y = MF_inv_simp, fill = id), 
             size = 4, shape = 23, colour = "white") +
  scale_fill_manual(values = col_pal[3:6]) +
  geom_point(data = ex_points_df, 
             mapping = aes(x = MF_ave, y = MF_inv_simp, colour = id), 
             size = 3.5, shape = 19, stroke = 1) +
  scale_colour_manual(values = col_pal[c(1, 2)]) +
  ylab("Inverse Simpson's index") +
  xlab("Average multifunctionality") +
  theme_meta() +
  theme(legend.position = "none")

# export to check the dimensions
ggsave(filename = "figures-paper-1/fig_S1_main.svg", 
       p1, units = "cm",
       width = 10, height = 10)

# plot the different column plots
col_plots_df <- dplyr::bind_rows(ex_points_df, rep_points_df)
col_plots_df$id <- 1:nrow(col_plots_df)

# pull into the long format
col_plots_df <- 
  col_plots_df |>
  dplyr::select(-MF_ave, -MF_inv_simp) |>
  tidyr::pivot_longer(cols = contains("F"),
                      names_to = "Function", 
                      values_to = "Value")

# loop over the six plots
for(i in 1:length( unique(col_plots_df$id) )) {
  
 p1 <- 
   ggplot(data = dplyr::filter(col_plots_df, id == i),
           mapping = aes(x = Function, y = Value)) +
    geom_col(width = 0.5, fill = col_pal[i], colour = "white") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1) ) +
    ylab("Function value (0-1)") +
    xlab(NULL) +
    theme_meta() +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 9))
  
  # export to check the dimensions
  ggsave(filename = paste0("figures-paper-1/", "fig_S1_sub", i, ".svg"), 
         p1, units = "cm",
         width = 4, height = 4.5)
  
}

# calculate the different multifunctionality metrics on the simulated data

# write a function with the relevant metrics
calculate_MF <- function(data, func.names) {
  
  data |>
    dplyr::mutate(`sum MF` = MF_sum(adf = data, vars = func.names, stand_method = "none"),
                  `ave. MF` = MF_av(adf = data, vars = func.names, stand_method = "none"),
                  `Pasari MF` = MF_pasari(adf = data, vars = func.names, stand_method = "none"),
                  `SAM MF` = MF_dooley(adf = data, vars = func.names,  stand_method = "none"),
                  `Simp. MF` = MF_inv_simpson(adf = data, vars = func.names, stand_method = "none"),
                  `Shannon MF` = MF_shannon(adf = data, vars = func.names, stand_method = "none"),
                  `ENF.Q0 MF` = multifunc::getMF_eff(data = data, vars = func.names, q = 0,
                                                     standardized = FALSE,
                                                     standardize_function = standardizeUnitScale,
                                                     D = NULL, tau = NULL),
                  `ENF.Q1 MF` = multifunc::getMF_eff(data = data, vars = func.names, q = 1,
                                                     standardized = FALSE,
                                                     standardize_function = standardizeUnitScale,
                                                     D = NULL, tau = NULL),
                  `thresh.30 MF` = MF_thresh(adf = data, vars = func.names, thresh = 0.3),
                  `thresh.70 MF` = MF_thresh(adf = data, vars = func.names, thresh = 0.7),
                  `PCA MF` = MF_pca(adf = data, vars = func.names) 
                  )
  
}

# use the function to calculate the MF metrics
df_mf <- calculate_MF(data = df, func.names = paste0("F", 1:4))

# convert to the long format
df_mf <-
  df_mf |>
  tidyr::pivot_longer(cols = contains(" MF"),
                      names_to = "Metric", 
                      values_to = "Value"
                      )

# add a column describing whether a metric is NA or not
df_mf <- 
  df_mf |>
  dplyr::mutate(NA_YN = ifelse(is.na(Value) | is.infinite(Value), "Y", "N"))

# how many NAs are there?
df_mf |>
  dplyr::filter(NA_YN == "Y")

# get a vector of the different names
metrics <- unique(df_mf$Metric)

# plot the different metrics
pal <- wesanderson::wes_palette("Darjeeling1", 50, type = "continuous")

# set-up the plots that need y axis labels
ylabs <- rep(0, length(metrics))
ylabs[c(1, 4, 7, 10)] <- 1

# set-up the plots that need x axis labels
xlabs <- rep(0, length(metrics))
xlabs[c(10, 11, 12)] <- 1

metric_plot_list <- vector("list", length = length(metrics))
for(i in 1:length(metrics)) {
  
  ylab <- if(ylabs[i] == 0) {
    NULL
  } else {
    "Inverse Simpson's index"
  }
  
  xlab <- if(xlabs[i] == 0) {
    NULL
  } else {
    "Average multifunctionality"
  }
  
  metric_plot_list[[i]] <- 
    ggplot(data = dplyr::filter(df_mf, Metric == metrics[i]), 
         mapping = aes(x = MF_ave, y = MF_inv_simp, 
                       colour = Value, shape = NA_YN)) +
    geom_jitter(mapping = aes(x = MF_ave, y = MF_inv_simp), 
                size = 1.2, alpha = 0.2, shape = 16, width = 0.01) +
    scale_shape_manual(values = c(1, 8)) +
    scale_colour_gradientn(colours = pal) +
    labs(colour = metrics[i]) +
    guides(color = guide_colourbar(frame.colour = "black", 
                                   ticks.colour = NA,
                                   title.vjust = 2.5,
                                   barwidth = 0.75,
                                   barheight = 3)) +
    ylab(ylab) +
    xlab(xlab) +
    theme_meta() +
    theme(legend.position = c(0.75, 0.3),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  
}


p1 <- 
  cowplot::plot_grid(plotlist = metric_plot_list,
                     align = "hv", axis = "tblr",
                     nrow = 4, ncol = 3,
                     label_size = 11,
                     labels = letters[1:12],
                     label_fontface = "plain"
                     )

# export to check the dimensions
ggsave(filename = "figures-paper-1/fig_S2.svg", 
       p1, units = "cm",
       width = 20, height = 25)

### END
