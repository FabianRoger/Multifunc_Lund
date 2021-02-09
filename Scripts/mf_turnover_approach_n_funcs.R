
library(devtools)
install_github("jebyrnes/multifunc")
library(multifunc)

# test this package using the BIODEPTH data
data(all_biodepth)

# function names
allVars<-qw(biomassY3, root3, N.g.m2,  light3, N.Soil, wood3, cotton3)
allVars

# index of the functions
varIdx<-which(names(all_biodepth) %in% allVars)
varIdx

# subset the germany data
germany<-subset(all_biodepth, all_biodepth$location=="Germany")

# make variable lists
vars<-whichVars(germany, allVars)
vars

species<-relevantSp(germany,26:ncol(germany))
species

spIDX <- which(names(germany) %in% species) #in case we need these


### overlap approach

# reflect the soil nutrient function to make sure it increases when nutrient depletion is high
germany$N.Soil<- -1*germany$N.Soil + max(germany$N.Soil, na.rm=T)

# use the sAICfun function to get coefficients for one funciton
# species: vector of species names
# germany: dataset with species presence-absence data and the function values

spList <- sAICfun("biomassY3", species, germany)
spList

# we apply this function for all functions

# vars: vector of function names
# species: vector of species names
# germany: data.frame with species presence, absence and function values
redund <- getRedundancy(vars, species, germany)

# this outputs an effect matrix i.e. -1, 0 or 1 corresponding to neg, neu or positive effects for each function
redund

# plot the num. functions by fraction of the species pool needed
posCurve <- divNeeded(redund, type="positive")
posCurve$div <- posCurve$div/ncol(redund)

negCurve<-divNeeded(redund, type="negative")
negCurve$div<-negCurve$div/ncol(redund)
negCurve

getOverlapSummary(redund, 
                  m=2, type = "positive", 
                  index = "sorensen", denom = "set")[1]



# write a function to automate this

# three inputs:
# vars: vector of function names
# species: vector of species names
# data: data.frame with species abundances/presence-absences and function data
# output: either "prop_species" for proportion of species affecting a function
# or "overlap" which reports the mean pairwise overlap in species positive and negative effects
# between all functions

turnover_aic <- function(func.names, species.names, adf.data, output = "prop_species") {
  
  if(! "multifunc" %in% installed.packages()[,1]) stop(
    "this function requires the multifunc package to be installed"
  )
  
  # use getRedundancy to get the effect of different species on each function
  redund.dat <- getRedundancy(vars = func.names, species = species.names, data = adf.data)
  
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
  
  # calculate average sorensen overlap among all pairs of functions
  pos.overlap <- getOverlapSummary(redund.dat, 
                                   m=2, type = "positive", 
                                   index = "sorensen", denom = "set")[1]
  
  neg.overlap <- getOverlapSummary(redund.dat, m=2, 
                                   type = "negative", 
                                   index = "sorensen", denom = "set")[1]
  
  overlap_effects <- 
    data.frame(direction = c("positive", "negative"),
               mean_sorensen_overlap = c(pos.overlap, neg.overlap))
  
  # specify what to output
  if (output == "prop_species") {
    return(species_effects)
  } 
  
  else if(output == "overlap") {
    return(overlap_effects)
  } 
  
  else {print("error, specify an output")}
  
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

# wrap this loop into a function
# run this on different data.sets

# potentially also add the observed values...

n = 5

null.out <- vector("list", length = n)
for (i in 1:n){
  
  # randomise the function positions relative to the species using row.randomiser
  data.random <- row.randomiser(func.names = vars, 
                                species.names = species, 
                                adf.data = germany)
  
  # apply the approach to the randomised data
  turnover.out <- turnover_aic(func.names = vars, 
                               species.names = species, 
                               adf.data = data.random,
                               output = "prop_species")
  
  # write the n_func proportion of the species pool to a list
  null.out[[i]] <- turnover.out
  
}

# write this output into a dataframe
null.func.dat <- dplyr::bind_rows(null.out, .id = "null.rep")
View(null.func.dat)

# get observed value
x <- turnover_aic(func.names = vars, species.names = species, adf.data = germany,
                  output = "prop_species")
x

ggplot() +
  geom_smooth(data = null.func.dat %>% filter(effect_direction == "positive"),
              mapping = aes(x = nfunc, y = div, group = null.rep),
              method = "lm", se = FALSE) +
  geom_smooth(data = x %>% filter(effect_direction == "positive"),
              mapping = aes(x = nfunc, y = div), colour = "red",
              method = "lm")



























