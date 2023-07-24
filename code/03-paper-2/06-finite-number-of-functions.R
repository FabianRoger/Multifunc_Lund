#'
#' @title Finite number of functions analaysis
#' 
#' @description Simulates an ecosystem under the assumption that it
#' has some finite number of functions present. It then tests how measuring
#' different numbers of functions leads to better estimates of the effect
#' of biodiversity (or any other variable) on multifunctionality.
#'

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load relevant scripts
source("code/helper-plotting-theme.R")
source("code/03-paper-2/02-helper-simulate-functions.R")

# set a seed for reproducibility
set.seed(54807)

# how many plots?
n <- 100

# how many functions?
n_func <- 20

# simulate a set of 20 functions
fsim <- sim_funcs(n_func = n_func, n = n, 
                  lambda = 10, 
                  mu_est = 0.1, sd_est = 0.1, 
                  error_sd = 0.5)

# bind the simulation data for efficient plotting
fsim_plot <- dplyr::bind_cols( dplyr::tibble(plot = 1:n, div = fsim[[1]]$div), 
                               dplyr::as_tibble(fsim[[2]]))

# rename the columns
names(fsim_plot) <- c("plot", "div", paste0("F", 1:n_func))

# pull into the long format
fsim_plot <- 
  fsim_plot |>
  tidyr::pivot_longer(cols = contains("F"),
                      names_to = "Function", 
                      values_to = "Value")

p1 <- 
  ggplot(data = fsim_plot,
       mapping = aes(x = div, y = Value, group = Function)) +
  geom_line(stat = "smooth", method = "lm", size = 0.5, alpha = 0.3) +
  ylab("Function value (0-1)") +
  xlab("Species richness") +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "black", size = 9.5),
        axis.text.y = element_text(colour = "black", size = 9.5),
        axis.title.x = element_text(size = 10.5),
        axis.title.y = element_text(size = 10.5))
plot(p1)

# plot the relationship between average EMF and diversity
fsim_mu <- dplyr::tibble(div = fsim[[1]]$div, ave_EMF = apply(fsim[[2]], 1, mean))

p2 <- 
  ggplot(data = fsim_mu,
       mapping = aes(x = div, y = ave_EMF)) +
  geom_jitter(width = 0.1, shape = 1, alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", size = 0.5, colour = "black") +
  ylab("True average EMF") +
  xlab("Species richness") +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "black", size = 9.5),
        axis.text.y = element_text(colour = "black", size = 9.5),
        axis.title.x = element_text(size = 10.5),
        axis.title.y = element_text(size = 10.5))
plot(p2)

p12 <- 
  cowplot::plot_grid(p1, p2, labels = c("a", "b"),
                     label_fontface = "plain", label_size = 11
                     )
plot(p12)

ggsave(filename = "figures-paper-2/fig_3.svg", p12,
       units = "cm", width = 15, height = 7.5)

# calculate the true average multifunctionality slope

# get the population slope
mu <- lm(apply(fsim[[2]], 1, mean) ~ fsim[[1]]$div)
mu <- summary(mu)

# extrac the true coefficients
true_B <- mu$coefficients[2,][["Estimate"]]
print(true_B)
true_SE <- mu$coefficients[2,][["Std. Error"]]
print(true_SE)

# get a vector of the number of functions to sample
sample_funcs <- rep(c(2, 6, 10, 14, 18), each = 10)

# get the individual slope of each function
samp_B <- vector(length = length(sample_funcs) )
samp_SE <- vector(length = length(sample_funcs) )
for(i in 1:length(sample_funcs)) {
  
  funcs <- sample(1:ncol(fsim[[2]]), sample_funcs[i])
  
  x <- lm(apply(fsim[[2]][,funcs], 1, mean) ~ fsim[[1]]$div)
  y <- summary(x)
  samp_B[i] <- y$coefficients[2,][["Estimate"]]
  samp_SE[i] <- y$coefficients[2,][["Std. Error"]]
  
}

# pull into a data.frame
samp_df <- dplyr::tibble(id = 1:length(samp_B),
                         n_funcs = sample_funcs,
                         slope = samp_B,
                         slope_SE = samp_SE
                         )

# add some random noise to the n_funcs variable for plotting
samp_df$n_funcs_dodge <- samp_df$n_funcs + rnorm(n = nrow(samp_df), 0, 0.2)

# pull the true effect into a data.frame
true_df <- dplyr::tibble(n_funcs = n_func,
                         slope = true_B,
                         slope_SE = true_SE)

p3 <- 
  ggplot() +
  geom_hline(yintercept = true_df$slope, linetype = "dashed", colour = "red") +
  geom_point(data = samp_df,
             mapping = aes(x = n_funcs_dodge, y = slope),
             shape = 16, size = 1.5, alpha = 0.3) +
  geom_errorbar(data = samp_df,
                mapping = aes(x = n_funcs_dodge, 
                              ymin = slope-slope_SE, ymax = slope+slope_SE),
                width = 0, alpha = 0.25) +
  geom_point(data = true_df,
             mapping = aes(x = n_funcs, y = slope),
             shape = 18, size = 3, alpha = 0.75, colour = "red") +
  geom_errorbar(data = true_df,
                mapping = aes(x = n_funcs, 
                              ymin = slope-slope_SE, ymax = slope+slope_SE),
                width = 0, alpha = 0.75, colour = "red") +
  scale_x_continuous(limits = c(0, 21), breaks = c(unique(sample_funcs), 20) ) +
  xlab("Number of functions measured") +
  ylab("Slope est. (+-SE)") +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "black", size = 9.5),
        axis.text.y = element_text(colour = "black", size = 9.5),
        axis.title.x = element_text(size = 10.5),
        axis.title.y = element_text(size = 10.5))
plot(p3)

ggsave(filename = "figures-paper-2/fig_4.svg", p3,
       units = "cm", width = 7.5, height = 7.5)

### END
