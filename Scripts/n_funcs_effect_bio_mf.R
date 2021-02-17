
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Does the number of functions matter?

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)

# load functions from important scripts
source(here("Scripts/MF_functions_collated.R"))
source(here("Scripts/function_plotting_theme.R"))


# BEF slope and n-functions

# get the list of matrices of function combinations
function.combinations <- function(vector.func.names) {
  nested.list.matrices <- vector("list", length = (length(vector.func.names)-1) )
  for (i in 2:length(vector.func.names)){
    nested.list.matrices[[i-1]] <- combn(x = vector.func.names, m = i)
  }
  return(nested.list.matrices)
}

# write a function to flatten an individual part of a nested list
flatten.list.matrices <- function(nested.list.matrices){
  unnested.list <- split(nested.list.matrices, col(nested.list.matrices)) 
  names(unnested.list) <- NULL 
  return(unnested.list)
}

# combine function.combinations and flatten.list.matrices and loop over each n-func
get.function.combinations <- function(function.names){
  
  # get list of matrices with function combinations
  list.func.matrix <- function.combinations(vector.func.names = function.names)
  
  # flatten the first matrix in the list
  list.combination <- flatten.list.matrices(nested.list.matrices = list.func.matrix[[1]])
  
  # loop over this and bind into a list
  for (i in 2:length(list.func.matrix)){
    x <- flatten.list.matrices(nested.list = list.func.matrix[[i]])
    list.combination <- c(list.combination, x)
  }
  return(list.combination)
}


# using the defined functions:

# set up a function to calculate the realised diversity-multifunctionality slope
# for each combination of a set of ecosystem functions

# adf.data: data.frame where rows are plots and which contains ecosystem functions of interest
# mf.func.names: vector of function names to consider
# standardise_funcs = TRUE i.e. standardise functions (z-score) before calculating multifunctionality
# covariate_name = name of the variable to regress against the multifunctional BEF slope

get_BEF_mf_est <- function(adf.data, 
                           mf.func.names, 
                           standardise_funcs = TRUE,
                           covariate_name = "realised_diversity") {
  
  if (length(mf.func.names) < 2) {
    stop("must have more than two functions to adequately examine the effect of changing the number of functions")
  }
  
  # use the function to get a list of all function combinations
  list.func.names <- get.function.combinations(function.names = mf.func.names)
  
  # set an output list
  list.out <- vector("list", length(list.func.names))
  
  # loop over each set of function names
  for (i in 1:length(list.func.names)) {
    
    # get the vector of function names
    sample.func.names <- list.func.names[[i]]
    
    # standardise functions if TRUE
    if (standardise_funcs == TRUE) {
      adf.data <- 
        adf.data %>%
        mutate(across(.cols = all_of(sample.func.names) , ~standardise(.) ) )
      }
    
    # calculate the multifunctionality metrics
    data.mf <- multifunc_calculator(adf = adf.data, vars = sample.func.names,
                                    mf.functions = c("MF_sum", "MF_av", "MF_pasari", "single_threshold_mf", "single_threshold_mf"),
                                    mf.names = c("sum_MF", "ave._MF", "Pasari_MF", "thresh.30_MF", "thresh.70_MF"),
                                    add.args = list(NA, NA, NA, c(thresh = 0.3), c(thresh = 0.7))
               )
    
    # subset the multifunctionality metric names
    mf.names <- names(data.mf)[grepl(pattern = "MF", x = names(data.mf))]
    
    # for each multifunctionality metric, calculate the BEF-slope
    bef_mf_slope <- 
      sapply(mf.names, function(x){
        lm.x <- lm(reformulate(covariate_name, as.name(x) ), data = data.mf)
        lm.x$coefficients[2]
      })
    names(bef_mf_slope) <- NULL # remove names from the vector
    
    # bind this into a data.frame
    df.out <- 
      data.frame(func.comb.id = i,
                 number_of_functions = length(sample.func.names),
                 multifunctionality_metric = mf.names,
                 realised_diversity_mf_est = bef_mf_slope)
    
    list.out[[i]] <- df.out
  }
  
  bind_rows(list.out)
  
}


### load the Jena data

# load the cleaned Jena data
jena.raw <- read_csv(file = here("data/jena_data_cleaned.csv"))

# define variable groups
var.names <- names(jena.raw)

# (1) get species names
spp.p <- ( grepl("+\\.+", var.names) & nchar(var.names) == 7 )
spp.names <- var.names[spp.p]

# (2) get site identifiers
site.id <- c("year", "sowndiv", "plotcode", "realised_diversity")

# (3) get function names
jena.func.names <- var.names[!(var.names %in% ( c(spp.names, site.id) ))  ]

# subset the site.id and function names
jena.dat <- 
  jena.raw %>%
  select(all_of( c(site.id, jena.func.names) ))

# get realised diversity - number of functions from Jena data
jena.n.func <- get_BEF_mf_est(adf.data = jena.dat, 
                              mf.func.names = jena.func.names, 
                              standardise_funcs = TRUE,
                              covariate_name = "realised_diversity")


# plot the Jena data
p1 <- 
  ggplot(data = jena.n.func,
         mapping = aes(x = number_of_functions, 
                       y = realised_diversity_mf_est)) +
  geom_jitter(alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  ylab("realised diversity - MF (est)") +
  xlab("number of functions") +
  facet_wrap(~multifunctionality_metric, scales = "free") +
  theme_meta() +
  theme(strip.background =element_rect(fill = "white", colour = "black"))

p1

ggsave(filename = here("Figures/fig_3_jena.png"), plot = p1,
       width = 16, height = 12, units = "cm", dpi = 300)
  

### load the simulated data

# load the raw data
sim.raw <- read_csv(file = here("data/multifunctionality_data.csv"))

# define the variable groups
sim.vars <- names(sim.raw)

# subset the function names from that 
sim.func.names <- sim.vars[ (grepl("F_", sim.vars) & nchar(sim.vars) == 3) ]

# get the id variables
sim.id.vars <- c("patch", "time", "richness", "total_abundance", "local.sp.pool")

# subset the id variables and functions
sim.dat <- 
  sim.raw %>%
  select(all_of( c(sim.id.vars, sim.func.names) ))

sim.n.func <- 
  lapply(split(sim.dat, sim.raw$sim.id), function(data) {
  
  get_BEF_mf_est(adf.data = data, 
                 mf.func.names = sim.func.names, 
                 standardise_funcs = FALSE,
                 covariate_name = "richness")
  
})

# bind this list into a single data.frame
sim.n.func.dat <- bind_rows(sim.n.func, .id = "sim.id")
sim.n.func.dat$sim.id <- as.numeric(sim.n.func.dat$sim.id)

# load the simulated data parameters
params <- read_csv(file = here("data/parameters_sim.csv"))

# create a simulation category variable
# reps <- length(unique(params$rep.id))
# id <- unique(paste(params$a_mean, params$w.shape, params$w.scale, sep = "_"))
# id <- LETTERS[1:length(id)]

# add this to the params data
# params$sim.group <- rep(id, each = reps)

# join the parameter data to the n_function data
sim.n.func.dat <- full_join(sim.n.func.dat, params, by = "sim.id")


# plot the simulated data
p2 <- 
  ggplot(data = sim.n.func.dat,
       mapping = aes(x = number_of_functions, 
                     y = realised_diversity_mf_est,
                     group = sim.id,
                     colour = sim.group )) +
  geom_smooth(method = "lm", se = FALSE) +
  ylab("realised diversity - MF (est)") +
  xlab("number of functions") +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  facet_wrap(~multifunctionality_metric, scales = "free") +
  theme_meta() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        strip.background =element_rect(fill = "white", colour = "black")) 

ggsave(filename = here("Figures/fig_4_sim.png"), plot = p2,
       width = 16, height = 13.5, units = "cm", dpi = 300)






