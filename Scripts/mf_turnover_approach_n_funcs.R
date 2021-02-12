
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
  random.row.ids <- sample(x = 1:func.mat.nrow , size = func.mat.nrow , replace = FALSE)
  
  # randomise the function matrix
  func.mat.random <- func.mat[random.row.ids, ]
  
  # get a matrix of species data
  spec.mat <- adf.data[, species.names]
  
  # bind this randomised data together
  adf.data.random <- cbind(spec.mat, func.mat.random)
  
  return(adf.data.random)
  
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



# test the turnover approach using the simulated data

# can it detect positive and negative effects on functions?
# if so, at what rate are there errors?

# next step: wrap the simulations into a list to run for each simulation

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(here)

# choose a simulation id
s <- 11

# read some of the simulated data
par.dat <- read_csv(here("data/parameters_sim.csv"))
abun.dat <- read_csv(here("data/raw_abundance_data.csv"))
mf.dat <- read_csv(here("data/multifunctionality_data.csv"))
func.dat <- read_csv(here("data/function_values_per_species.csv"))

# check the parameters
par.dat[s,]

# subset the chosen simulation id
abun.dat <- 
  abun.dat %>%
  filter(sim.id == s) %>%
  pivot_wider(names_from = "species", values_from = "abundance")

spp.present <- sapply(abun.dat[, grepl("sp_", names(abun.dat)) ], function(x) sum(ifelse(x > 0, 1, 0)))
spp.present <- names(spp.present[spp.present > 0])

abun.dat <- 
  abun.dat %>%
  select(sim.id, patch, time, all_of(spp.present))

# subset the multifunctionality data
mf.dat <- 
  mf.dat %>%
  filter(sim.id == s) %>%
  select(-richness, -total_abundance, -local.sp.pool, -ends_with("MF"))

# join the abundance data to the function data
mf.dat <- full_join(mf.dat, abun.dat, by = c("sim.id", "patch", "time"))


# subset the function data
func.dat <- 
  func.dat %>%
  filter(sim.id == s)

# subset the species that are present
func.dat <- func.dat[func.dat$species %in% spp.present, ]

# test if order is preserved: if FALSE then order is preserved
any(func.dat$species != spp.present)

# list of function names
f.names <- 
  mf.dat %>%
  select(starts_with("F_")) %>%
  names(.)


list.out <- vector("list", length = length(f.names))
for (i in 1:length(f.names)){
  
  # get correct function
  f.vals <- func.dat[[f.names[i]]]
  names(f.vals) <- func.dat$species
  
  # get species effect on each function using abundances
  r.abun <- getRedundancy(vars = f.names[i], species = spp.present, data = mf.dat)
  r.abun <- sapply(r.abun, function(x)(x))
  
  # get species effect on each function using presence absence
  mf.dat2 <- 
    mf.dat %>%
    mutate(across(.cols = all_of(spp.present), ~if_else(. > 0, 1, 0) ) )
  
  r.pa <- getRedundancy(vars = f.names[i], species = spp.present, data =  mf.dat2)
  r.pa <- sapply(r.pa, function(x)(x))
  
  # proportion of incorrect directions?
  # vec1: real function values
  # vec2: inferred species effect (-1, 0 or 1)
  # ignores species that were non-significant
  
  # outputs proportion of incorrectly inferred directions
  compare.directions.vectors <- function(vec1, vec2){
    x <- sum(vec1[vec2 < 0] > 0)
    y <- sum(vec1[vec2 > 0] < 0)
    (x + y)/sum(vec2 != 0)
  }
  
  # compare presence absence versus abundance accuracy
  error.pa <- compare.directions.vectors(vec1 = f.vals, vec2 = r.pa)
  
  # ifelse(f.vals < 0, -1, 1)[r.pa != 0]
  # r.pa[r.pa != 0]
  
  error.abun <- compare.directions.vectors(vec1 = f.vals, vec2 = r.abun)
  
  # ifelse(f.vals < 0, -1, 1)[r.abun != 0]
  # r.abun[r.abun != 0]
  
  # calculate proportion of incorrect directions
  cor.pa <- cor(f.vals, r.pa, method = "spearman")
  cor.abun <- cor(f.vals, r.abun, method = "spearman")
  
  df.out <- 
    data.frame(func.id = f.names[i],
               error.pa = error.pa,
               error.abun = error.abun,
               cor.pa = cor.pa,
               cor.abun = cor.abun)
  
  list.out[[i]] <- df.out
}

bind_rows(list.out)




























