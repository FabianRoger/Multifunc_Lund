
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Add multifunctionality to the abundance data from the stachova and leps (2010) model

# select scripts to call functions from
library(here)
source(here("Scripts/multispecies_lotka_volterra_model/stachova_leps_2010_LK_model_function.R"))
source(here("Scripts/MF_functions_collated.R"))

# combinations of simulations

# start with four simulation types:
# (1) stabilising competition + high specialisation
# (2) stabilising competition + low specialisation
# (3) strong competition + high specialisation
# (4) strong competition + low specialisation
# (5) equal competition + high specialisation
# (6) equal competition + low specialisation

# run each of these types 5 times
# stick with 500 time-steps to start with

# what to output?
# (1) multifunctionality dataframe
# (2) raw species abundances
# (3) competition coefficients
# (4) k values
# (5) r values
# (6) function matrix

# libraries for data manipulation
library(dplyr)
library(tidyr)
library(readr)

# set-up a data.frame of parameters for each run

# fixed parameters

# lotka-volterra model
lsp = c(4, 6, 8, 10, 12)
reps = 10
rsp = 20
t_steps = 10
n0 = 20
a_min = 0
a_max = 0.75
sim.comp = "sym"

# set the number of functions
func.n = 9

# set the probability that a species negatively affects a function
prob.neg = 0.1

# varying parameters

# number of replicates per simulation
sim.reps <- 5

# set the average level of interspecific competition
a_mean <- c(0.15, 0.7, 0.05)

# set the standard deviation of intraspecific competition
a_sd <- c(0.05, 0.05, 0)

# set the intraspecific competition value
a_spp <- c(1, 1, 0.05)

# set the min and max k-values
k_min = c(20, 20, 50)
k_max = c(150, 150, 50)

# set the min and max r-values
r_min = c(0.1, 0.1, 0.25)
r_max = c(0.5, 0.5, 0.25)

# set the level of species specialisation using parameters of the Weibull distribution
w.shape <- c(0.5, 3)
w.scale <- c(0.25, 0.5)

# set-up a data.frame of parameter combinations for each simulation
params <- data.frame(a_mean = a_mean, a_sd = a_sd, a_spp = a_spp,
           k_min = k_min, k_max = k_max,
           r_min = r_min, r_max = r_max)

params <- merge(params, data.frame(w.shape = w.shape))
params <- merge(data.frame(rep.id = 1:sim.reps), params )

params$w.scale <- rep(w.scale, each = length(a_mean)*sim.reps)

# set an id column
params <- cbind(sim.id = 1:nrow(params), params)

# create a simulation category variable
reps <- length(unique(params$rep.id))
id <- unique(paste(params$a_mean, params$w.shape, params$w.scale, sep = "_"))
id <- LETTERS[1:length(id)]

# add this to the params data
params$sim.group <- rep(id, each = reps)

# write the parameter combinations to a .csv file
# write_csv(x = params, here("data/parameters_sim.csv"))

# run each simulation

# create output lists
mf.dataframe <- vector("list", length = nrow(params))
raw.abundances <- vector("list", length = nrow(params))
alpha.dat <- vector("list", length = nrow(params))
k.vals <- vector("list", length = nrow(params))
r.vals <- vector("list", length = nrow(params))
function.vals <- vector("list", length = nrow(params))

for (i in 1:nrow(params)) {
  
  # simulate a set of communities
  df.x <- 
    s_l_function(lsp = lsp,
                 reps = reps,
                 rsp = rsp,
                 t_steps = t_steps,
                 n0 = n0,
                 a_mean = params$a_mean[i], a_sd = params$a_sd[i], 
                 a_min = a_min, a_max = a_max, 
                 a_spp = params$a_spp[i], sim.comp = sim.comp,
                 k_min = params$k_min[i], k_max = params$k_max[i],
                 r_min = params$r_min[i], r_max = params$r_max[i]
    )
  
  # get the raw species abundance data
  df.spp <- df.x$data.raw
  
  # make a species list
  spp.list <- unique(df.spp$species)
  
  # make a function name list
  func.names <- paste("F_", 1:func.n, sep = "")
  
  # generate a function matrix for the relationship between each species abundance and each function
  func.mat <- 
    lapply(spp.list, function(x){
      x <- rweibull(n = func.n, shape = params$w.shape[i], scale = params$w.scale[i])
      y <- x - quantile(x, prob.neg)
      z <- (x*y)
      round(z, digits = 4)
    })
  
  func.mat <- data.frame(do.call(rbind, func.mat))
  func.mat <- cbind(spp.list, func.mat)
  names(func.mat) <- c("species", func.names)
  
  # join the function data to the species abundance data
  multi.func <- full_join(df.spp, func.mat, by = "species")
  
  # multiply abundance by these function values
  multi.func <- 
    multi.func %>%
    mutate(across(.cols = func.names, ~(.*abundance) ))
  
  # function to standardise and then translate the functions
  func.std.trans <- function(x){
    y <- (x - mean(x))/sd(x)
    return(y + abs(min(y)))
  }
  
  # calculate each function as the sum of all species-specific function values
  multi.func <- 
    multi.func %>%
    group_by(patch) %>%
    summarise(across(.cols = func.names, sum), .groups = "drop") %>%
    mutate(across(.cols = func.names, func.std.trans ))
  
  # subset the functions from the multi.func data.frame
  adf.func <- select(multi.func, func.names)
  
  # calculate the multifunctionality metrics and add them to the multi.func data.frame
  multi.func <- 
    multi.func %>%
    mutate(`scal. MF` = MF_jing(adf = adf.func, vars = func.names),
           `sum MF` = MF_sum(adf = adf.func, vars = func.names),
           `ave. MF` = MF_av(adf = adf.func, vars = func.names),
           `MESLI MF` = MF_mesli(adf = adf.func, vars = func.names),
           `Pasari MF` = MF_pasari(adf = adf.func, vars = func.names),
           `SAM MF` = MF_dooley(adf = adf.func, vars = func.names),
           `ENF MF` = as.numeric(hill_multifunc(adf = adf.func, vars = func.names, scale = 1, HILL = TRUE)) ,
           `Simp. MF` = MF_simpsons_div(adf = adf.func, vars = func.names),
           `Manning.30 MF` = manning_multifunc(adf = adf.func, vars = func.names, thresh = 0.3),
           `Manning.50 MF` = manning_multifunc(adf = adf.func, vars = func.names, thresh = 0.5),
           `Manning.70 MF` = manning_multifunc(adf = adf.func, vars = func.names, thresh = 0.7),
           `thresh.30 MF` = single_threshold_mf(adf = adf.func, vars = func.names, thresh = 0.3),
           `thresh.50 MF` = single_threshold_mf(adf = adf.func, vars = func.names, thresh = 0.5),
           `thresh.70 MF` = single_threshold_mf(adf = adf.func, vars = func.names, thresh = 0.7),
           `Slade.10.90 MF` = MF_slade(adf = adf.func, vars = func.names, A_quant = 0.10, B_quant = 0.90),
           `Slade.40.60 MF` = MF_slade(adf = adf.func, vars = func.names, A_quant = 0.40, B_quant = 0.60),
           `PCA MF` = pca_multifunc(adf = adf.func, vars = func.names, standardise = FALSE) )
  
  # join this multifunctionality data to the summary data
  mf.sim <- full_join(df.x$data.summary, multi.func, by = c("patch"))
  
  # output the following:
  
  # (1) multifunctionality data.frame
  mf.dataframe[[i]] <- mf.sim
  
  # (2) raw species abundances
  raw.abundances[[i]] <- df.x$data.raw
  
  # (3) competition coefficients
  alpha.dat[[i]] <- as.data.frame(df.x$spp.info$competition.coefficients)
  
  # (4) carrying capacities
  k.vals[[i]] <- data.frame(K = df.x$spp.info$k.vals)
  
  # (5) intrinsic growth rates
  r.vals[[i]] <- data.frame(R = df.x$spp.info$r.vals)
  
  # (6) function matrix
  function.vals[[i]] <- func.mat
}


# put these outputs into a list
sim.outputs <- 
  list(mf.dataframe,
       raw.abundances,
       alpha.dat, 
       k.vals,
       r.vals,
       function.vals)


# set a vector of data.names
data.names <- c("multifunctionality_data",
                "raw_abundance_data",
                "competition_coefficients",
                "carrying_capacities",
                "intrinsic_growth_rates",
                "function_values_per_species")


# for each of the outputs, write into a data.frame
data.loc <- here("data")
for (i in 1:length(sim.outputs)) {
  
  x.dat <- bind_rows(sim.outputs[[i]], .id = "sim.id")
  write_csv(x = x.dat, paste(data.loc, "/", data.names[i], ".csv", sep = "" ))
  
}

