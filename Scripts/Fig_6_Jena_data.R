
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Implement the turnover approach using the multifunc package and an analogue to Meyer et al. (2017)'s null model

# arguments for the functions:

# vars: vector of function names
# species: vector of species names
# data: data.frame with species abundances/presence-absences and function data
# output: either "prop_species" for proportion of species affecting a function
# or "overlap" which reports the mean pairwise overlap in species positive and negative effects
# between all functions

# n: number of null replicates

turnover_aic <- function(func.names, species.names, adf.data, output = "prop_species") {
  
  if(! "multifunc" %in% installed.packages()[,1]) stop(
    "this function requires the multifunc package to be installed"
  )
  
  # load the multifunc package
  library(multifunc)
  
  # use getRedundancy to get the effect of different species on each function
  redund.dat <- getRedundancy(vars = func.names, species = species.names, data = adf.data)
  
  # specify what to output
  if (output == "prop_species") {
    
    # calculate the number of species that positively affect different numbers of functions
    posCurve <- divNeeded(redund.dat, type = "positive")
    posCurve$div <- posCurve$div/ncol(redund.dat)
    row.names(posCurve) <- NULL
    posCurve$effect_direction <- "positive"
    
    # calculate the number of species that negative affect different numbers of functions
    negCurve<-divNeeded(redund.dat, type = "negative")
    negCurve$div<-negCurve$div/ncol(redund.dat)
    row.names(negCurve) <- NULL
    negCurve$effect_direction <- "negative"
    
    # bind these positive and negative effects together
    species_effects <- rbind(posCurve, negCurve)
    species_effects <- cbind(row.id = 1:nrow(species_effects), species_effects)
    
    return(species_effects)
  }
  
  else if(output == "overlap") {
    
    # calculate average sorensen overlap among all pairs of functions
    
    # positive function effects
    pos.overlap <- getOverlapSummary(redund.dat, 
                                     m=2, type = "positive", 
                                     index = "sorensen", denom = "set")[1]
    
    # negative function effects
    neg.overlap <- getOverlapSummary(redund.dat, m=2, 
                                     type = "negative", 
                                     index = "sorensen", denom = "set")[1]
    
    # bind this into a data.frame
    overlap_effects <- 
      data.frame(direction = c("positive", "negative"),
                 mean_sorensen_overlap = c(pos.overlap, neg.overlap))
    
    return(overlap_effects)
  } 
  
  else { print("error, specify an output") }
  
}

# function to mix-up relationship between function and species
row.randomiser <- function(func.names, species.names, adf.data) {
  
  # subset a matrix of functions
  func.mat <- adf.data[, func.names]
  func.mat.nrow <- nrow(func.mat) # calculat the number of rows in that matrix
  
  # use sample to get random row ids from the function matrix
  random.row.ids <- sample(x = 1:func.mat.nrow, size = func.mat.nrow , replace = FALSE)
  
  # randomise the function matrix
  func.mat.random <- func.mat[random.row.ids, ]
  
  # get a matrix of species data
  spec.mat <- adf.data[, species.names]
  
  # bind this randomised data together
  adf.data.random <- cbind(spec.mat, func.mat.random)
  
  return(as.data.frame(adf.data.random))
  
}


# use the defined functions to calculate turnover for randomised data and observed data
turnover_aic_null <- function(func.names, 
                              species.names, 
                              adf.data, 
                              output = "prop_species",
                              n = 100) {
  
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires the dplyr package to be installed"
  )
  
  # randomise the data and calculate turnover indices n times
  null.out <- vector("list", length = n)
  for (i in 1:n){
    
    # randomise the function positions relative to the species using row.randomiser
    data.random <- row.randomiser(func.names = func.names, 
                                  species.names = species.names, 
                                  adf.data = adf.data)
    
    # apply the approach to the randomised data
    turnover.out <- turnover_aic(func.names = func.names, 
                                 species.names = species.names, 
                                 adf.data = data.random,
                                 output = output)
    
    # write the n_func proportion of the species pool to a list
    null.out[[i]] <- turnover.out
  }
  
  # bind the loop output into a data.frame
  null.out <- dplyr::bind_rows(null.out, .id = "null.rep")
  
  # get observed value
  obs.out <- turnover_aic(func.names = func.names, 
                          species.names = species.names, 
                          adf.data = data.random,
                          output = output)
  
  # bind this into a list
  return(list(null.out, obs.out))
  
}


# apply the aic_turnover_null function to the Jena data

# load relevant libraries
library(dplyr)
library(readr)
library(here)

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

# convert to presence-absence instead of relative abundance
jena.dat <- 
  jena.raw %>%
  mutate(across(.cols = all_of(spp.names), ~if_else(.x > 0, 1, 0)))

# test the turnover_aic_null function and plot the results
# this function can take a long time to run (e.g. 20-30 minutes)
# jena.turn <- turnover_aic_null(func.names = jena.func.names, 
# species.names = spp.names, 
# adf.data = jena.dat,
# output = "prop_species",
# n = 1000)

# output these files to a .csv
# write_csv(x = jena.turn[[1]], file = here("data/jena_null_model_data.csv"))
# write_csv(x = jena.turn[[2]], file = here("data/jena_null_observed_data.csv"))

# read in the null model data
jena.null <- read_csv(file = here("data/jena_null_model_data.csv"))
jena.obs <- read_csv(file = here("data/jena_null_observed_data.csv"))

# check the max of the null data to see if it is below 1
max(jena.null$div)

# calculate the 97.5% and 2.5% quantiles for each number of functions
# plot this as well
jena.null.q <- 
  jena.null %>% 
  group_by(effect_direction, null.rep, nfunc) %>%
  summarise(mean_div = mean(div, na.rm = TRUE), .groups = "drop") %>%
  group_by(effect_direction, nfunc) %>%
  summarise(quant_97.5 = quantile(mean_div, probs = c(0.975)),
            quant_2.5 = quantile(mean_div, probs = c(0.025)),
            .groups = "drop")

# calculate mean and sd of the observed data
jena.obs.s <- 
  jena.obs %>%
  group_by(effect_direction, nfunc) %>%
  summarise(mean_div = mean(div, na.rm = TRUE),
            sd_div = sd(div, na.rm = TRUE))

# plot the null expectations versus the observed data
library(ggplot2)
g1 <- 
  ggplot() +
  geom_ribbon(data = jena.null.q,
              mapping = aes(x = nfunc, ymax = quant_97.5, ymin = quant_2.5),
              alpha = 0.25) +
  geom_jitter(data = jena.obs,
              mapping = aes(x = nfunc, y = div), 
              width = 0.1, colour = "red", shape = 16, alpha = 0.1, size = 2) +
  geom_line(data = jena.obs.s,
            mapping = aes(x = nfunc, y = mean_div), size = 1, colour = "red") +
  facet_wrap(~effect_direction, scales = "free") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  ylab("number of species") +
  xlab("number of functions") +
  theme_meta()

g1

ggsave(filename = here("Figures/aic_turnover_null.png"), plot = g1,
       width = 11, height = 7.5, units = "cm", dpi = 450)

# randomly assigning functions to plots generates the same pattern as the empirical data
# why is this the case?

# maybe it is because the turnover approach is anti-conservative...

### END
