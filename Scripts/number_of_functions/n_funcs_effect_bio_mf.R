
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
                 n.func.id = paste(sample.func.names, collapse = ""),
                 number_of_functions = length(sample.func.names),
                 multifunctionality_metric = mf.names,
                 realised_diversity_mf_est = bef_mf_slope)
    
    list.out[[i]] <- df.out
  }
  
  bind_rows(list.out)
  
}


### list with the ecological drift model

f.names <- names(mod.df)[grepl(pattern = "F_", names(mod.df))]

sim.n.func <- 
  lapply(split(mod.df, mod.df$mod_id), function(data) {
  
  get_BEF_mf_est(adf.data = data, 
                 mf.func.names = f.names, 
                 standardise_funcs = TRUE,
                 covariate_name = "local_species_pool")
  
})

sim.n.df <- bind_rows(sim.n.func, .id = "mod_id")
head(sim.n.df)

ggplot(data = sim.n.df,
       mapping = aes(x = number_of_functions, y = realised_diversity_mf_est, colour = mod_id)) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~multifunctionality_metric, scales = "free") +
  theme_classic() +
  theme(legend.position = "none")
  




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
  geom_jitter(width = 0.1, alpha = 0.3, shape = 16, size = 1) +
  geom_smooth(method = "lm", se = FALSE, 
              colour = viridis::viridis(n = 1, begin = 0.5, end = 0.5, option = "D")) +
  ylab("multifunctional BEF-slope") +
  xlab("number of functions") +
  facet_wrap(~multifunctionality_metric, scales = "free") +
  theme_meta() +
  theme(strip.background = element_rect(fill = "white", colour = "black"))


p1

ggsave(filename = here("Figures/fig_9_jena.png"), plot = p1,
       width = 16, height = 12, units = "cm", dpi = 450)



# questions to answer:

# why does the multifunctional BEF-slope increase with number of functions for the Pasari approach increase?
# why does the multifunctional BEF-slope increase and decrease with number of functions for the summing approach?

# make a copy of the sim.dat data
mf_slope_dat <- 
  lapply(split(sim.dat, sim.raw$sim.id), function(data) {
  
  # get combinations of function names
  list.func.names <- get.function.combinations(function.names = sim.func.names)
  
  func.out <- vector("list")
  
  # loop over each set of function names
  for (i in 1:length(list.func.names)) {
    
    # get the vector of function names
    sample.func.names <- list.func.names[[i]]
    
    site.mods <- 
      data %>%
      select("patch", "time", "richness", "total_abundance", "local.sp.pool")
    
    df.x <- 
      data %>%
      select(all_of(sample.func.names) )
    
    site.mods$mean_funcs <- apply(df.x, MARGIN = 1, mean)
    site.mods$sum_funcs <- apply(df.x, MARGIN = 1, sum)
    site.mods$sd_funcs <- apply(df.x, MARGIN = 1, sd)
    
    site.mods$n.func <- length(sample.func.names)
    site.mods$n.func.id <- paste(sample.func.names, collapse = "")
    
    func.out[[i]] <- site.mods
    
  }
  
  return(bind_rows(func.out))
  
  } )

# bind this into a data.frame
mf_slope_dat <- bind_rows(mf_slope_dat, .id = "sim.id")
View(head(mf_slope_dat))


# what changes with the sum_MF when we consider more functions?

# calculate average and range in sum of functions for each random draw
dat.sum <- 
  mf_slope_dat %>%
  group_by(sim.id, n.func, n.func.id) %>%
  summarise(average_sum = mean(sum_funcs, na.rm = TRUE),
            range_sum = max(sum_funcs, na.rm = TRUE) - min(sum_funcs, na.rm = TRUE),
            .groups = "drop")

# get slope between biodiversity and sum multifunctionality
sum.slopes <- 
  sim.n.func.dat %>%
  filter(multifunctionality_metric == "sum_MF") %>%
  mutate(sim.id = as.character(sim.id)) %>%
  select(sim.id, n.func.id, number_of_functions, realised_diversity_mf_est) %>%
  rename(n.func = number_of_functions)

mf.sum.exp <- 
  full_join(dat.sum, sum.slopes, by = c("sim.id", "n.func.id", "n.func"))

# overall sum_MF increases with number of functions
# this, alone, should only affect the intercept
ggplot(data = mf.sum.exp,
       mapping = aes(x = n.func, y = average_sum) ) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  facet_wrap(~sim.id, scales = "free") +
  theme_classic()

# but varying the number of functions increases with range of sum_MF
# however, the richness range stays the same and this will affect the slope
# this is not biological though
ggplot(data = mf.sum.exp,
       mapping = aes(x = n.func, y = range_sum) ) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  facet_wrap(~sim.id, scales = "free") +
  theme_classic()

# does the average sum or range sum correlate with the BEF slope?
ggplot(data = mf.sum.exp,
         mapping = aes(x = average_sum, y = realised_diversity_mf_est)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~sim.id, scales = "free") +
  theme_meta()

ggplot(data = mf.sum.exp,
       mapping = aes(x = range_sum, y = realised_diversity_mf_est)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~sim.id, scales = "free") +
  theme_meta()

# now we have to plot some of the raw data
id.input <- 27

mf.sum.exp %>%
  filter(sim.id == id.input) %>%
  ggplot(data = .,
         mapping = aes(x = n.func, y = realised_diversity_mf_est )) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_meta()

# func_ids <- 
  # mf_slope_dat %>%
  # filter(sim.id == id.input) %>%
  # group_by(n.func) %>%
  # slice_sample(n = 1) %>%
  # pull(n.func.id)

mf_slope_dat %>%
  filter(sim.id == id.input) %>%
  # filter(n.func.id %in% func_ids) %>%
  ggplot(data = .,
         mapping = aes(x = richness, y = sum_funcs, 
                       group = n.func.id, colour = n.func )) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, size = 0.75) +
  theme_meta()


# what about the Pasari approach?
dat.pas <- 
  mf_slope_dat %>%
  group_by(sim.id, n.func, n.func.id) %>%
  summarise(average_sd = mean(sd_funcs, na.rm = TRUE),
            average_mean = mean(mean_funcs, na.rm = TRUE),
            range_mean = diff(range(mean_funcs)),
            range_sd = diff(range(sd_funcs)),
            range_pasari = diff( range( (mean_funcs - sd_funcs)  )  ),
            .groups = "drop")

# get slope between biodiversity and sum multifunctionality
pas.slopes <- 
  sim.n.func.dat %>%
  filter(multifunctionality_metric == "Pasari_MF") %>%
  mutate(sim.id = as.character(sim.id)) %>%
  select(sim.id, n.func.id, number_of_functions, realised_diversity_mf_est) %>%
  rename(n.func = number_of_functions)

mf.pas.exp <- 
  full_join(dat.pas, pas.slopes, by = c("sim.id", "n.func.id", "n.func"))
  

# does the range in Pasari_MF change with n_funcs?
head(mf.pas.exp)

ggplot(data = mf.pas.exp,
       mapping = aes(x = n.func, y = range_pasari)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  facet_wrap(~sim.id) +
  theme_meta()

# why does this happen?
ggplot(data = mf.pas.exp,
       mapping = aes(x = n.func, y = average_sd)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  facet_wrap(~sim.id) +
  theme_meta()

ggplot(data = mf.pas.exp,
       mapping = aes(x = n.func, y = average_mean)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  facet_wrap(~sim.id) +
  theme_meta()

# now we have to plot some of the raw data
id.input <- 12

# func_ids <- 
  # mf_slope_dat %>%
  # filter(sim.id == id.input) %>%
  # group_by(n.func) %>%
  # slice_sample(n = 1) %>%
  # pull(n.func.id)

mf_slope_dat %>%
  filter(sim.id == id.input) %>%
  # filter(n.func.id %in% func_ids) %>%
  ggplot(data = .,
         mapping = aes(x = richness, y = (mean_funcs - sd_funcs),
                       group = n.func.id, colour = as.character(n.func) )) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_meta() +
  facet_wrap(~as.character(n.func))






