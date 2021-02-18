
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Implement the turnover approach using the multifunc package and an analogue to Meyer et al. (2017)'s null model

# install and load the multifunc package
# library(devtools)
# install_github("jebyrnes/multifunc")
library(multifunc)

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

# test the turnover_aic_null function using BIODEPTH data
data(all_biodepth)
allVars <- qw(biomassY3, root3, N.g.m2,  light3, N.Soil, wood3, cotton3)

# subset out the german site
germany <- subset(all_biodepth, all_biodepth$location=="Germany")
head(germany)

# get vectors for functions and for species
vars <- whichVars(germany, allVars)
species <- relevantSp(germany,26:ncol(germany))

# test the turnover_aic_null function and plot the results
x <- turnover_aic_null(func.names = vars, 
                       species.names = species, 
                       adf.data = germany,
                       output = "prop_species",
                       n = 100)
                  
# plot the null expectations versus the observed data
library(ggplot2)
ggplot() +
  stat_smooth(data = x[[1]],
              mapping = aes(x = nfunc, y = div, group = null.rep),
              geom='line', alpha=0.025, size = 1, se=FALSE, method = "lm") +
  geom_smooth(data = x[[2]],
              mapping = aes(x = nfunc, y = div), colour = "red",
              method = "lm", se = FALSE) +
  facet_wrap(~effect_direction, scales = "free") +
  theme_classic()

# randomly assigning functions to plots generates the same pattern as the empirical data
# why is this the case?

# maybe it is because the turnover approach is anti-conservative


# test the turnover approach using the simulated data

# can it detect positive and negative effects on functions?
# if so, at what rate do errors arise?


# next step: wrap the simulations into a list to run for each simulation

# wrap the turnover and AIC approaches into pipelines
# calculate error associated with both of them

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(here)

# load the simulated data
par.dat <- read_csv(here("data/parameters_sim.csv"))
abun.dat <- read_csv(here("data/raw_abundance_data.csv"))
mf.dat <- read_csv(here("data/multifunctionality_data.csv"))
func.dat <- read_csv(here("data/function_values_per_species.csv"))


# this must run for each simulation

# choose a simulation id
s <- 11

# check the parameters
par.dat[s,]


# prepare the abundance data

# subset the chosen simulation id
abun <- 
  abun.dat %>%
  filter(sim.id == s) %>%
  pivot_wider(names_from = "species", values_from = "abundance")

# get a vector of species names that are present for the chosen simulation
spp.present <- sapply(abun[, grepl("sp_", names(abun)) ], function(x) sum(ifelse(x > 0, 1, 0)))
spp.present <- names(spp.present[spp.present > 0])

# subset the species that are present in the simulation
abun <- 
  abun %>%
  select(sim.id, patch, time, all_of(spp.present))


# prepare the function data

# subset the multifunctionality data
mf <- 
  mf.dat %>%
  filter(sim.id == s) %>%
  select(-richness, -total_abundance, -local.sp.pool, -ends_with("MF"))

# join the raw abundance data to the function data
# mf.dat contains function data and species abundances in one data.frame
mf.a <- full_join(mf, abun, by = c("sim.id", "patch", "time"))

# get species effect on each function using presence-absence
# first convert the species abundance data to presence-absence data
mf.pa <- 
  mf.a %>%
  mutate(across(.cols = all_of(spp.present), ~if_else(. > 0, 1, 0) ) )

# get a list of function names
f.names <- 
  mf.a %>%
  select(starts_with("F_")) %>%
  names(.)


# function trait data

# subset the functional trait data
func.dat <- 
  func.dat %>%
  filter(sim.id == s)

# subset the species that are present
func.dat <- func.dat[func.dat$species %in% spp.present, ]

# test if order is preserved: if FALSE then order is preserved
# this works because R tests element per element
# e.g. c(1, 2, 3) == c(1, 2, 3)
any(func.dat$species != spp.present)


# write two functions:

# (1) calculate species effects using the AIC-based turnover approach

# (2) calculate species effects using Gotelli et al. (2011)'s method

a.out <- vector("list", length = length(f.names))
pa.out <- vector("list", length = length(f.names))
for (i in 1:length(f.names)){
  
  # get species effect on each function using abundances
  # this outputs a vector of species effects (-1, 0 or 1) on the function i
  r.abun <- getRedundancy(vars = f.names[i], species = spp.present, data = mf.a)
  a.out[[i]] <- sapply(r.abun, function(x)(x))
  
  # get species effect on each function using presence-absence
  # this outputs a vector of species effects (-1, 0 or 1) on the function i
  r.pa <- getRedundancy(vars = f.names[i], species = spp.present, data =  mf.pa)
  pa.out[[i]] <- sapply(r.pa, function(x)(x))
  
  }

do.call(cbind, a.out)
do.call(cbind, pa.out)


# proportion of incorrect directions?
# true.function.vals: real function values (i.e. linear coefficient)
# inferred.effects: inferred species effect (-1, 0 or 1)
# ignores species that were non-significant

# get the function values for each species for the given function i
f.vals <- func.dat[[f.names[i]]]
names(f.vals) <- func.dat$species

# outputs proportion of IN-correctly inferred directions
compare.directions.vectors <- function(true.function.vals, inferred.effects){
  x <- sum(true.function.vals[inferred.effects < 0] > 0)
  y <- sum(true.function.vals[inferred.effects > 0] < 0)
  (x + y)/sum(inferred.effects != 0)
}

# compare presence absence versus abundance accuracy
# error.pa is the proportion of incorrectly inferred directions
error.pa <- compare.directions.vectors(true.function.vals = f.vals, inferred.effects = r.pa)

# ifelse(f.vals < 0, -1, 1)[r.pa != 0]
# r.pa[r.pa != 0]

error.abun <- compare.directions.vectors(true.function.vals = f.vals, inferred.effects = r.abun)

# ifelse(f.vals < 0, -1, 1)[r.abun != 0]
# r.abun[r.abun != 0]


### implement Gotelli et al. (2011)'s turnover metric

# use row.randomiser to reshuffle function rows
# then calculate D-values for the shuffled data

# set number of replicates
n <- 1000

# data: species presence-absence data

SES_score <- function(data, function_names, species_names, n_ran) {
  
  ran.d.vals <- vector("list", length = n_ran)
  for (i in 1:n_ran) {
    
    # randomise function values across plots
    data.random <- row.randomiser(func.names = function_names,
                                  species.names = species_names, 
                                  adf.data = data)
    
    # calculate difference values for each species for each function
    ran.d.vals[[i]] <-  
      data.random %>%
      pivot_longer(cols = species_names,
                   names_to = "species",
                   values_to = "pa") %>%
      group_by(species, pa) %>%
      summarise(across(.cols = all_of(function_names), mean ), .groups = "keep") %>%
      summarise(across(.cols = all_of(function_names), diff ), .groups = "drop")
    
  }
  
  # bind the list into a data.frame  
  ran.d.vals <- bind_rows(ran.d.vals, .id = "run")
  
  # calculate mean and sd of d-values of randomised data
  ran.d.vals <- 
    ran.d.vals %>%
    pivot_longer(cols = function_names,
                 names_to = "function_name",
                 values_to = "d_value") %>%
    arrange(species, function_name) %>%
    group_by(species, function_name) %>%
    summarise(mean_d_value = mean(d_value, na.rm = TRUE),
              sd_d_value = sd(d_value, na.rm = TRUE), .groups = "drop")
  
  # calculate observed d values
  obs.d.vals <- 
    data %>%
    pivot_longer(cols = species_names,
                 names_to = "species",
                 values_to = "pa") %>%
    group_by(species, pa) %>%
    summarise(across(.cols = all_of(function_names), mean ), .groups = "keep") %>%
    summarise(across(.cols = all_of(function_names), diff ), .groups = "drop") %>%
    pivot_longer(cols = starts_with("F_"),
                 names_to = "function_name",
                 values_to = "d_value_obs") %>%
    arrange(species, function_name)
  
  # join the observed d.vals to the randomised d.vals
  d.val.dat <- full_join(ran.d.vals, obs.d.vals, by = c("species", "function_name"))
  
  # calculate species importance scores
  sp.SES <- 
    d.val.dat %>%
    mutate(SES = (d_value_obs - mean_d_value)/(sd_d_value)) %>%
    mutate(SES_effect = if_else(SES > 2, 1, 
                                if_else(SES < -2, -1, 0))) %>%
    select(species, function_name, SES_effect) %>%
    pivot_wider(id_cols = "species",
                names_from = "function_name",
                values_from = "SES_effect")
  
  return(sp.SES)
  
}

SES_score(data = mf.pa, function_names = f.names, species_names = spp.present, n_ran = 10)






















