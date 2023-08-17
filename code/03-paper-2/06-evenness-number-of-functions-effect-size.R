#'
#' @title Why does the slope change with the number of functions for evenness
#' 
#' @description Simulate the case of considering more functions to illustrate
#' that the evenness among functions 
#'

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load relevant scripts
source("code/helper-plotting-theme.R")
source("code/03-paper-2/06-helper-gamfeldt-roger-functions.R")

# set a seed for reproducibility
set.seed(777)

# set the number of species and the number of functions
specnum <- 15
funcnum <- 9

# specify the distribution used to draw species' traits
distribution = "runif"

# get the species function matrix
FuncMat <- FunctionValue(specnum, funcnum, distribution, min = 0, max = 1)

# get function name and species name vectors
func_names <- as.character(unique(FuncMat$Functions))
spec_names <- as.character(unique(FuncMat$Species))

# set the maximum number of replicates of N functions and species combinations
maxrep <- 50

# set up the species combination matrix with the maximum replications
SpecMat <- SpeciesMatrix(specnum = specnum, maxrep = maxrep)

# set the method to calculate the function values in a plot based on species traits
method = "av"

# get the values of different functions for each species combination
AvFunc <- AverageFunction(SpecMat, FuncMat,
                          method = method,
                          CF = CF, 
                          compfunc = compfunc,
                          r = r)

# check the summary of the different variables
summary(AvFunc)

# empty dataframe to store results
RES_func <- tibble(div = numeric(),
                   MF_eff_Q1 = numeric(), 
                   MF_eff_Q2 = numeric(), 
                   EffN_Q1 = numeric(),
                   EffN_Q2 = numeric(), 
                   AV = numeric(),
                   nfunc = numeric(),
                   func_comb = numeric())

# loop over chosen subsets of all function of varying size
for (i in 2:funcnum) { 
  
  # all possible combination of i out of funcnum functions
  func_comb <- combn(func_names, i)
  
  # sample 50 random function combinations if more than 50 possible combinations
  if(ncol(func_comb) > maxrep) {
    func_comb <- func_comb[, sample(c(1:ncol(func_comb)), maxrep)]
  }
  
  # loop over all function combinations of size i
  for ( k  in seq_len(ncol(func_comb))) { 
    
    # effective number of functions of q = 1 EMF metric
    MF_eff_Q1 <- multifunc::getMF_eff(AvFunc, func_comb[ ,k], q = 1)
    
    # effective number of functions of q = 1 EMF metric
    MF_eff_Q2 <- multifunc::getMF_eff(AvFunc, func_comb[ ,k], q = 2)
    
    # get effective number of functions q = 1
    EffN_Q1 <- multifunc::eff_num_func(AvFunc, func_comb[ ,k], q = 1, standardized = FALSE)
    
    # get effective number of functions q = 2
    EffN_Q2 <- multifunc::eff_num_func(AvFunc, func_comb[ ,k], q = 2, standardized = FALSE)
    
    # get average MF
    AV <- multifunc::getStdAndMeanFunctions(AvFunc, func_comb[ ,k])
    
    temp <- tibble(div = AvFunc$Richness,
                   MF_eff_Q1 = MF_eff_Q1, 
                   MF_eff_Q2 = MF_eff_Q2, 
                   EffN_Q1 = EffN_Q1,
                   EffN_Q2 = EffN_Q2, 
                   AV = AV$meanFunction,
                   nfunc = i,
                   func_comb = k)
    
    RES_func <- rbind(RES_func, temp)
    
  }
  
}

# check the results
head(RES_func)

# standardise the effective number of functions
RES_func <- 
  RES_func |>
  dplyr::mutate(EffN_Q1_st = EffN_Q1 / nfunc,
                EffN_Q2_st = EffN_Q2 / nfunc)

# set-up the plotting options for the three plots
vars <- c("AV", "EffN_Q2", "MF_eff_Q2")
ylabs <- c("Average EMF", "Inv. Simpson EMF", "ENF-Q2 EMF")
xlabs <- c(NA, NA, "Species richness")
cols <- c("black", "white", "white")
sizes <- c(11, 2, 2)

# make an output list
plot_list <- vector("list", length = length(vars))

# loop over the different variables
for(i in 1:length(vars)) {
  
  p <- 
    ggplot(data = RES_func) +
    geom_jitter(mapping = aes_string(x = "div", y = vars[i]), 
                size = 0.5, colour = "grey", width = 0.1, shape = 1) +
    geom_smooth(mapping = aes_string(x = "div", y = vars[i], group = "func_comb"), 
                method = "lm", se = F, linewidth = 0.05, colour = "red", alpha = 0.1) +
    geom_smooth(mapping = aes_string(x = "div", y = vars[i]),
                method = "lm", se = F, size = 1, colour = "black") +
    xlab(if(is.na(xlabs[i])) { NULL } else { xlabs[i] } ) +
    ylab(ylabs[i]) +
    facet_wrap(~nfunc, nrow = 1, ncol = 8) +
    theme_meta() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = sizes[i], colour = cols[i]))
  
  plot_list[[i]] <- p
  
}

# combine the plots
p1 <- 
  cowplot::plot_grid(plotlist = plot_list,
                     nrow = 3, ncol = 1, align = "v",
                     labels = c("a", "b", "c"), label_fontface = "plain", 
                     label_size = 10)
plot(p1)

# calculate the slope of the relationship
RES_est <- 
  RES_func |>
  tidyr::pivot_longer(dplyr::one_of(c("AV", "EffN_Q2", "MF_eff_Q2")), names_to = "metric", values_to = "value") |> 
  dplyr::group_by(nfunc, func_comb, metric) |> 
  tidyr::nest() |>
  dplyr::mutate(
    model = purrr::map(data,  ~ lm(value ~ div, data = .)) 
  ) |> 
  dplyr::mutate(
    tidy_summary = purrr::map(model, broom::tidy)
  ) |> 
  tidyr::unnest(tidy_summary)

# get the div slope
RES_est <- dplyr::filter(RES_est, term == "div")

# set-up the plotting options for the three plots
vars <- c("AV", "EffN_Q2", "MF_eff_Q2")
ylabs <- c("Average EMF ~ Species richness (Est.)", 
           "Inv. Simpson EMF ~ Species richness (Est.)", 
           "ENF-Q2 EMF ~ Species richness (Est.)")
xlabs <- c(NA, "Number of functions", NA)

# make an output list
plot_list <- vector("list", length = length(vars))

# loop over the different variables
for(i in 1:length(vars)) {
  
  # subset the correct data
  data <- dplyr::filter(RES_est, metric == vars[i])
  
  # convert the number of functions to a character variable
  data$est_id <- with(data, paste(nfunc, func_comb, sep = "_"))
  
  p <- 
    ggplot(data = data) +
    geom_point(mapping = aes(x = nfunc, y = estimate, group = est_id), 
               shape = 1, alpha = 0.5, position = position_dodge(width = 0.5)) +
    geom_errorbar(mapping = aes(x = nfunc, 
                                ymin = estimate-std.error, ymax = estimate+std.error,
                                group = est_id),
                  width = 0, position = position_dodge(width = 0.5),
                  alpha = 0.2) +
    scale_x_continuous(breaks = c(2:9)) +
    geom_smooth(mapping = aes(x = nfunc, y = estimate),
                method = "lm", se = F, size = 0.5,colour = "red") +
    xlab(if(is.na(xlabs[i])) { "" } else { xlabs[i] } ) +
    ylab(ylabs[i]) +
    theme_meta()
  
  plot_list[[i]] <- p
  
}

# combine the plots
p1 <- 
  cowplot::plot_grid(plotlist = plot_list,
                     nrow = 1, ncol = 3, align = "v",
                     labels = c("a", "b", "c"), label_fontface = "plain", 
                     label_size = 10)
plot(p1)

### END
