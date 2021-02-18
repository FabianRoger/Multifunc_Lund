
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Implement the turnover approach using the multifunc package and an analogue to Meyer et al. (2017)'s null model

# install and load the multifunc package
# library(devtools)
# install_github("jebyrnes/multifunc")

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
library(multifunc)
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
g1 <- 
  ggplot() +
  stat_smooth(data = x[[1]],
              mapping = aes(x = nfunc, y = div, group = null.rep),
              geom='line', alpha=0.025, size = 1, se=FALSE, method = "lm") +
  geom_smooth(data = x[[2]],
              mapping = aes(x = nfunc, y = div), colour = "red",
              method = "lm", se = FALSE) +
  facet_wrap(~effect_direction, scales = "free") +
  ylab("number of species") +
  xlab("number of functions") +
  theme_meta()

ggsave(filename = here("Figures/aic_turnover_null.png"), plot = g1,
       width = 11, height = 7.5, units = "cm", dpi = 450)

# randomly assigning functions to plots generates the same pattern as the empirical data
# why is this the case?

# maybe it is because the turnover approach is anti-conservative...



### test the turnover approach using the simulated data

# can it detect positive and negative effects on functions?
# if so, at what rate do errors arise?


# functions to implement aic-based species effects and Gotelli et al. (2011) species effects

# 1. calculate species effects using the AIC-based turnover approach

# data: data.frame with functions and either species abundances or presence-absences as columns and samples as rows
# function_names: vector of function names
# species_names: vector of species names

AIC_sp <- function(data, function_names, species_names) {
  
  if(! "multifunc" %in% installed.packages()[,1]) stop(
    "this function requires the multifunc package to be installed"
  )
  
  # load the multifunc package
  library(multifunc)
  
  sp.effect.aic <- vector("list", length = length(function_names))
  for (i in 1:length(function_names)){
    
    # get species effect on each function using abundances
    # this outputs a vector of species effects (-1, 0 or 1) on the function i
    redun.out <- getRedundancy(vars = function_names[i], species = species_names, data = data)
    sp.effect.aic[[i]] <- sapply(redun.out, function(x)(x))
    
  }
  
  # prepare the output
  sp.effect.aic <- as.data.frame(do.call(cbind, sp.effect.aic))
  names(sp.effect.aic) <- function_names
  sp.effect.aic <- cbind(species = rownames(sp.effect.aic), data.frame(sp.effect.aic, row.names = NULL))
  
  # return the data.frame with species effects
  return(sp.effect.aic)
  
}

# AIC_sp(data = mf.pa, function_names = f.names, species_names = spp.present)


# 2. calculate species effects using Gotelli et al. (2011)'s method

# data: data.frame with functions and species presence-absences as columns and samples as rows
# function_names: vector of function names
# species_names: vector of species names
# n_ran: number of randomisations

SES_score <- function(data, function_names, species_names, n_ran) {
  
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires the dplyr package to be installed"
  )
  if(! "tidyr" %in% installed.packages()[,1]) stop(
    "this function requires the tidyr package to be installed"
  )
  
  # load dplyr and tidyr packages
  library(dplyr)
  library(tidyr)
  
  # randomise the data and calculate d-values for each randomisation
  ran.d.vals <- vector("list", length = n_ran)
  for (i in 1:n_ran) {
    
    # randomise function values across plots
    data.random <- row.randomiser(func.names = function_names,
                                  species.names = species_names, 
                                  adf.data = data)
    
    # calculate difference values for each species for each function
    ran.d.vals[[i]] <-  
      data.random %>%
      pivot_longer(cols = all_of(species_names),
                   names_to = "species",
                   values_to = "pa") %>%
      group_by(species, pa) %>%
      summarise(across(.cols = all_of(function_names), mean ), .groups = "drop") %>%
      group_by(species) %>%
      summarise(across(.cols = all_of(function_names), diff ), .groups = "drop")
    
  }
  
  # bind the list into a data.frame  
  ran.d.vals <- bind_rows(ran.d.vals, .id = "run")
  
  # calculate mean and sd of d-values of randomised data
  ran.d.vals <- 
    ran.d.vals %>%
    pivot_longer(cols = all_of(function_names),
                 names_to = "function_name",
                 values_to = "d_value") %>%
    arrange(species, function_name) %>%
    group_by(species, function_name) %>%
    summarise(mean_d_value = mean(d_value, na.rm = TRUE),
              sd_d_value = sd(d_value, na.rm = TRUE), .groups = "drop")
  
  # calculate observed d values
  obs.d.vals <- 
    data %>%
    pivot_longer(cols = all_of(species_names),
                 names_to = "species",
                 values_to = "pa") %>%
    group_by(species, pa) %>%
    summarise(across(.cols = all_of(function_names), mean ), .groups = "drop") %>%
    group_by(species) %>%
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

# SES_score(data = mf.pa, function_names = f.names, species_names = spp.present, n_ran = 10)


# write functions to test reliability

# Q1: for detected species effects, how many species effects are in the incorrect direction?

# true.function.vals: real function values (i.e. linear coefficient)
# inferred.effects: inferred species effect (-1, 0 or 1)

# note that this function ignores non-significant species effects (i.e. effect = 0)

# outputs proportion of IN-correctly inferred directions
compare.directions <- function(true.function.vals, inferred.effects){
  n.incorrect.negative <- sum(true.function.vals[inferred.effects < 0] > 0)
  n.incorrect.positive <- sum(true.function.vals[inferred.effects > 0] < 0)
  n.incorrect.total <- (n.incorrect.negative + n.incorrect.positive)
  
  return( n.incorrect.total/sum(inferred.effects != 0) )
}


# Q2: do the approaches detect the correct direction for species 
# above and below certain percentiles of function values?

# true.function.vals: real function values (i.e. linear coefficient)
# inferred.effects: inferred species effect (-1, 0 or 1)
# low.q: lower quantile
# upp.q: upper quantile

# note that this function only includes values below the lower quantile that are negative
# and values above the upper quantile that are positive

# proportion_undetected_negatives: proportion of true negatives that were falsely specified
# proportion_undetected_positives: proportion of true positives that were falsely specified

# true negatives/positives are based on the quantile values
# we assume that function values below the low.q and upp.q (that are negative and positive respectively)
# indicate true negative effects

# outputs proportion of IN-correctly inferred directions
compare.quantiles <- function(true.function.vals, inferred.effects, low.q = 0.10, upp.q = 0.90){
  
  # get upper and lower quantiles of true.function.vals
  quantiles <- quantile(true.function.vals, probs = c(low.q, upp.q))
  
  # test if lowest quantiles show negative effects
  true.negative <- (true.function.vals < quantiles[1]) & (true.function.vals < 0) 
  undetected.negatives <- sum(inferred.effects[true.negative] < 0)/sum(true.negative)
  
  # test if high quantiles show positive effects
  true.positive <- (true.function.vals > quantiles[2]) & (true.function.vals > 0) 
  undetected.positives <- sum(inferred.effects[true.positive] > 0)/sum(true.positive)
  
  return(c("proportion_undetected_negatives" = undetected.negatives, 
           "proportion_undetected_positives" = undetected.positives))
  
}


# pull comparisons into a function that outputs a cleaned data.frame

# aic_dat: output from AIC_sp
# ses_dat: output from SES_score

# true_function_data: simulated function coefficients
# function_names: vector of function names

# note that these different objects must correspond to the same input data
compare_turnover_approaches <- function(aic_dat, ses_dat, true_function_data, function_names) {
  
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires the dplyr package to be installed"
  )
  if(! "tidyr" %in% installed.packages()[,1]) stop(
    "this function requires the tidyr package to be installed"
  )
  
  # load dplyr and tidyr packages
  library(dplyr)
  library(tidyr)
  
  # Q1:
  # aic approach
  q1.aic <- mapply(compare.directions, true_function_data[, function_names], aic_dat[, function_names])
  
  # ses approach
  q1.ses <- mapply(compare.directions, true_function_data[, function_names], ses_dat[, function_names])
  
  # Q2:
  # aic approach
  q2.aic <- mapply(compare.quantiles, true_function_data[, function_names], aic_dat[, function_names])
  
  # ses approach
  q2.ses <- mapply(compare.quantiles, true_function_data[, function_names], ses_dat[, function_names])
  
  # Q3:
  # aic approach
  q3.aic <- mapply(cor, true_function_data[, function_names], aic_dat[, function_names], method = "spearman")
  
  # ses approach
  q3.ses <- mapply(cor, true_function_data[, function_names], ses_dat[, function_names], method = "spearman")
  
  # create a vector of aic vs. ses
  approach <- c("aic", "ses", "aic", "aic", "ses", "ses", "aic", "ses")
  
  # create a vector of labels
  e.labs <- 
    c(rep("prop_incorrect_direction", 2),
      rep(c("proportion_undetected_negatives", "proportion_undetected_positives"), 2),
      rep("spearman_correlation", 2))
  
  # pull into a list
  q.list <- 
    list(q1.aic, q1.ses, 
         q2.aic[1,], q2.aic[2 ,], q2.ses[1, ], q2.ses[2,], 
         q3.aic, q3.ses)
  
  sp.effect.performance <- cbind(approach = approach, effect = e.labs, bind_rows(q.list, .id = "num") )
  sp.effect.performance <- as.data.frame(sp.effect.performance)
  
  # make this data longer
  sp.effect.performance <- 
    sp.effect.performance %>%
    pivot_longer(cols = all_of(f.names),
                 names_to = "function_name",
                 values_to = "value")
  
  return(sp.effect.performance)
}


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

# get a vector of simulation ids
sim.ids <- unique(par.dat$sim.id)

# loop over each simulation
compare.sims <- vector("list", length = length(sim.ids))
for (i in 1:length(sim.ids)) {
  
  # prepare the abundance data
  
  # subset the chosen simulation id
  abun <- 
    abun.dat %>%
    filter(sim.id == sim.ids[i]) %>%
    pivot_wider(names_from = "species", values_from = "abundance")
  
  # get a vector of species names that are present for the chosen simulation
  sp.present <- sapply(abun[, grepl("sp_", names(abun)) ], function(x) sum(ifelse(x > 0, 1, 0)))
  sp.present <- names(sp.present[sp.present > 0])
  
  # subset the species that are present in the simulation
  abun <- 
    abun %>%
    select(sim.id, patch, time, all_of(sp.present))
  
  # prepare the function data
  
  # subset the multifunctionality data
  mf <- 
    mf.dat %>%
    filter(sim.id == sim.ids[i]) %>%
    select(-richness, -total_abundance, -local.sp.pool, -ends_with("MF"))
  
  # join the raw abundance data to the function data
  # mf.dat contains function data and species abundances in one data.frame
  mf.a <- full_join(mf, abun, by = c("sim.id", "patch", "time"))
  
  # get species effect on each function using presence-absence
  # first convert the species abundance data to presence-absence data
  mf.pa <- 
    mf.a %>%
    mutate(across(.cols = all_of(sp.present), ~if_else(. > 0, 1, 0) ) )
  
  # get a list of function names
  f.names <- 
    mf.pa %>%
    select(starts_with("F_")) %>%
    names(.)
  
  # function trait data
  
  # subset the functional trait data
  func <- 
    func.dat %>%
    filter(sim.id == sim.ids[i])
  
  # subset the species that are present
  func <- func[func$species %in% sp.present, ]
  
  # test if order is preserved: if FALSE then order is preserved
  # this works because R tests element per element
  # e.g. c(1, 2, 3) == c(1, 2, 3)
  if( any(func$species != sp.present) == TRUE){
    print("WARNING! species order is not preserved")
  }
  
  # implement the aic-based approach to get species effects on each function
  aic.x <- AIC_sp(data = mf.pa, function_names = f.names, species_names = sp.present)
  
  # implement the ses-based approach to get species effects on each function
  ses.x <- SES_score(data = mf.pa, function_names = f.names, species_names = sp.present, n_ran = 1000)
  
  # compare these two approaches for each of the 30 simulated data.sets
  compare.sims[[i]] <- 
    compare_turnover_approaches(aic_dat = aic.x, 
                                ses_dat = ses.x, 
                                true_function_data = func, 
                                function_names = f.names)
  
}

# bind into a data.frame
compare.df <- bind_rows(compare.sims, .id = "simulation.id")

# plot the error rates

# load the plotting theme:
source(here("Scripts/function_plotting_theme.R"))

p1 <- 
  compare.df %>%
  group_by(simulation.id, approach, effect) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE), .groups = "drop") %>%
  ggplot(data = .,
         mapping = aes(x = approach, y = mean_value, colour = simulation.id)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(mapping = aes(x = approach, 
                              ymin = mean_value-sd_value, 
                              ymax = mean_value+sd_value, 
                              colour = simulation.id),
                width = 0.05,
                position=position_dodge(width=0.5)) +
  facet_wrap(~effect, scales = "free") +
  ylab("") +
  scale_colour_viridis_d(end = 0.9, option = "C") +
  theme_meta() +
  theme(legend.position = "right",
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))

ggsave(filename = here("Figures/turnover_reliability.png"), plot = p1,
       width = 18.5, height = 16, units = "cm", dpi = 450)











