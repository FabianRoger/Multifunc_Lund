
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
# stick with 50 time-steps to start with

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

# choose which multifunctionality metrics to calculate

# set the function names to call
mf.metric.function <- c("MF_jing", "MF_sum", "MF_av", 
                        "MF_pasari", "MF_mesli", "MF_dooley", "hill_multifunc", "MF_simpsons_div",
                        "manning_multifunc", "manning_multifunc", "manning_multifunc",
                        "single_threshold_mf", "single_threshold_mf", "single_threshold_mf",
                        "MF_slade", "MF_slade",
                        "pca_multifunc" )

# set the metric names to be outputted
mf.metric.names = c("scal._MF", "sum_MF", "ave._MF", 
                    "Pasari_MF", "MESLI_MF", "SAM_MF", "ENF_MF", "Simp_MF",
                    "Manning.30_MF", "Manning.50_MF", "Manning.70_MF",
                    "thresh.30_MF", "thresh.50_MF", "thresh.70_MF",
                    "Slade.10.90 MF", "Slade.40.60 MF",
                    "PCA_MF")

# set the additional arguments list for each function
additional.mf.args = list(NA, NA, NA, 
                          NA, NA, NA, NA, NA,
                          c(thresh = 0.3), c(thresh = 0.5), c(thresh = 0.7),
                          c(thresh = 0.3), c(thresh = 0.5), c(thresh = 0.7),
                          c(A_quant = 0.10, B_quant = 0.90), c(A_quant = 0.40, B_quant = 0.60),
                          NA )



# set-up a data.frame of parameters for each run

# fixed parameters

# lotka-volterra model
lsp = c(1, 2, 4, 6, 8, 10, 12)
reps = 20
rsp = 20
t_steps = 50
n0 = 100
a_min = 0
a_max = 1
sim.comp = "sym"
a_scale = 1

# set the number of functions
func.n = 9

# set the probability that a species negatively affects a function
prob.neg = 0.1

# varying parameters

# number of replicates per simulation
sim.reps <- 5

# set the average level of interspecific competition
a_mean <- c(0.25, 0.5, 0.75)

# set the standard deviation of intraspecific competition
a_sd <- c(0.1, 0.1, 0)

# set the intraspecific competition value
a_spp <- c(1, 1, 0.75)

# set the min and max k-values
k_min = c(100, 100, 250)
k_max = c(300, 300, 250)

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
reps.sim <- length(unique(params$rep.id))
id <- unique(paste(params$a_mean, params$w.shape, params$w.scale, sep = "_"))
id <- LETTERS[1:length(id)]

# add this to the params data
params$sim.group <- rep(id, each = reps.sim)

# write the parameter combinations to a .csv file
write_csv(x = params, here("data/parameters_sim.csv"))

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
      x[x > 10] <- rnorm(n = sum(x > 10), mean = 10, sd = 2) # constrain values to around 10
      y <- x - quantile(x, prob.neg)
      z <- round(y, digits = 4)
    })
  
  func.mat <- data.frame(do.call(rbind, func.mat))
  func.mat <- cbind(spp.list, func.mat)
  names(func.mat) <- c("species", func.names)
  
  # join the function data to the species abundance data
  multi.func <- full_join(df.spp, func.mat, by = "species")
  
  # multiply abundance by these function values
  # add small normally distributed error (i.e. sd = 0.1)
  multi.func <- 
    multi.func %>%
    mutate(across(.cols = all_of(func.names), ~(.*abundance) )) %>%
    mutate(across(.cols = all_of(func.names), ~(. + rnorm(n = length(.), mean = 0, sd = 0.1)) ) )
  
  # calculate each function as the sum of all species-specific function values
  # standardise the functions
  multi.func <- 
    multi.func %>%
    group_by(patch) %>%
    summarise(across(.cols = all_of(func.names), sum), .groups = "drop") %>%
    mutate(across(.cols = all_of(func.names), standardise ))
  
  # subset the functions from the multi.func data.frame
  adf.func <- select(multi.func, func.names)
  
  # calculate the multifunctionality metrics and add them to the multi.func data.frame
  
  # calculate the multifunctionality metrics
  multi.func <- 
    multifunc_calculator(adf = multi.func, vars = func.names,
                         mf.functions = mf.metric.function,
                         mf.names = mf.metric.names,
                         add.args = additional.mf.args)
  
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

# check a few random outputs
check.id <- 26

params[check.id,]
mf.dataframe[[check.id]]$richness
mf.dataframe[[check.id]]$local.sp.pool
length(unique(mf.dataframe[[check.id]]$patch))


# output these data as .csv files into the 'data' folder

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

### END

