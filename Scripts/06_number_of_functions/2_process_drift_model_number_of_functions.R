
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Process neutral model data to test the number of functions

# load relevant libraries
library(here)
library(readr)
library(MASS)

# clear the current memory
rm(list = ls())

# process this model data

# read the model data from the data folder
# this list is generated using the drift_model_number_functions.R script
drift.mod.list <- read_csv(here("data/drift_model_n_functions.csv"))
head(drift.mod.list)
dim(drift.mod.list)

# make the analysis repeatable
set.seed(1543875)

# write a function to generate function matrices for the species
func_mat_generator <- function(mu = 5, sigma = 1, 
                               n_funcs = 5, sp_n = 12,
                               type = "equal") {
  
  # generate a list of species names
  sp <- paste("sp_", 1:sp_n, sep = "")
  
  # function matrix 1: all species equal
  
  if (type == "equal") {
    
    x <- rnorm(n = n_funcs, mean = mu, sd = sigma)
    y <- lapply(sp, function(i) x)
    z <- do.call("rbind" , y)
    print(mean((cor(t(z)))))
    
  } else if (type == "random") {
    
    # function matrix 2: random species values
    y <- lapply(sp, function(i) rnorm(n = n_funcs, mean = mu, sd = sigma))
    z <- do.call("rbind" , y)
    print(mean((cor(t(z)))))
    
  } else {
    
    stop("choose appropriate type of function matrix")
    
  }
  
  # package into a decent-looking output
  func_mat <- data.frame(z)
  f.names <- paste("F_", 1:n_funcs, sep = "")
  names(func_mat) <- f.names
  func_mat$species <- sp
  func_mat <- func_mat[,c("species", f.names)]
  
  # quantity to output
  return(func_mat)
  
}

# load scripts with processing functions
source(here("Scripts/01_general_functions/process_model_data.R"))

# generate an equal matrix (i.e. all species identical)
f.equal <- func_mat_generator(mu = 5, sigma = 1, n_funcs = 5, sp_n = 12, type = "equal")

# generate a random matrix (i.e. species differ in their functional contributions)
f.random <- func_mat_generator(mu = 5, sigma = 1, n_funcs = 5, sp_n = 12, type = "random")

# compile these matrices into a list
func.list <- list(f.equal, f.random)

# for each dataset in drift.mod.list and for each function matrix in func.list
# process the data

func.reps <- vector("list", length = length(func.list))
for (j in 1:length(func.list)) {
  
  # process the simulated cluster using: process_sim_data
  df.proc <- 
    process_sim_data(model_data = drift.mod.list, 
                     func.mat =  func.list[[j]], 
                     time_final = TRUE, 
                     species_abun = "raw")
  
  # add this to a list
  func.reps[[j]] <- df.proc
  
}

# bind this list into a data.frame  
mod.out <- bind_rows(func.reps, .id = "function_matrix")
View(mod.out)

# add a unique model-ID variable
mod.out <- 
  mod.out %>%
  mutate(mod_id = paste(drift_parameter, function_matrix, model_run, sep = ".")) %>%
  select(mod_id, drift_parameter, function_matrix, model_run:F_5)

# check the model parameters
length(unique(mod.out$mod_id))
head(mod.out)

# write this model data.frame into a .csv file so the model does not have to be re-run

# make a folder to export the cleaned data
if(! dir.exists(here("data"))){
  dir.create(here("data"))
}

write_csv(x = mod.out, file = here("data/drift_model_n_functions_processed.csv"))

### END
