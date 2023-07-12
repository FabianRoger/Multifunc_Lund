#'
#' @title Test whether slope changes with number of functions
#' 
#' @description Here we simulate relationships between biodiversity and
#' ecosystem multifunctionality (EMF) and then we test whether the relationship
#' between biodiversity and EMF changes by considering different numbers
#' of functions.
#'

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load relevant functions
source("code/03-paper-2/02-helper-simulate-functions.R")
source("code/03-paper-2/03-helper-slope-estimate-slope-n-functions.R")
source("code/helper-plotting-theme.R")
source("code/helper-univariate-mf-functions.R")

# simulate a set of functions
fsim <- sim_funcs(n_func = 3, n = 100, 
                  lambda = 10, 
                  mu_est = 0.25, sd_est = 0.5, 
                  error_sd = 0.5)

# bind the simulation data for efficient plotting
fsim_plot <- dplyr::bind_cols( fsim[[1]], fsim[[2]])

# rename the columns
names(fsim_plot) <- c("plot", "div", paste0("F", 1:ncol(fsim[[2]])) )

# pull into the long format
fsim_plot <- 
  fsim_plot |>
  tidyr::pivot_longer(cols = contains("F"),
                      names_to = "Function", 
                      values_to = "Value")

p1 <- 
  ggplot(data = fsim_plot,
         mapping = aes(x = div, y = Value, group = Function)) +
  geom_smooth(method = "lm", size = 0.5, colour = "grey", se = FALSE) +
  ylab("Function value (0-1)") +
  xlab("Species richness") +
  theme_test() +
  theme(axis.text = element_text(colour = "black"))
plot(p1)

adf_dat <- dplyr::as_tibble(fsim[[2]])
names(adf_dat) <- paste0("F", 1:ncol(fsim[[2]]))

# estimate the slope between biodiveristy and functioning
# different numbers of functions
n_func_est(adf = adf_dat, vars = names(adf_dat), 
           div = fsim[[1]]$div, 
           metric = "thresh_30")


# add a few non-linear functions




