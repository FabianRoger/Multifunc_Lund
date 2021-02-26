
# Prepare data for multivariate-diversity interactions models

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# I always use relative paths using the package 'here'
# if you are not using an R-project
# the best thing to do for this code to run is is to set your working directory now
setwd(dir = "")

# then load the here package
library(here)

# then we make a folder called data in the current working directory
if(! dir.exists(here("data"))){
  dir.create(here("data"))
}

# add the two full datasets that I sent into this data folder
# 1. multifunctionality_data.csv
# 2. raw_abundance_data.csv

# if you do not want to work with relative paths and here, 
# just be sure to specify the correct paths when you read and write the data


### load the raw data (if not using here() specify path manually)

# import the .csv file called: 'multifunctionality_data.csv'
mf.dat <- read_csv(here("data/multifunctionality_data.csv"))

# remove univariate multifunctionality metrics
mf.dat <- 
  mf.dat %>%
  select(-ends_with("MF"))

# import the .csv file called 'raw_abundance_data.csv'
a.dat <- read_csv(here("data/raw_abundance_data.csv"))


### view the data for one random simulation

# you can run this code block (line 42-63) multiple times to see different simulations
ran.sim <- 
  mf.dat %>%
  filter(sim.id == sample(unique(mf.dat$sim.id), 1))

# plot species richness versus total abundance (i.e. number of individuals)
ggplot(data = ran.sim,
       mapping = aes(x = richness, y = total_abundance)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  theme_classic()

# plot species richness versus each of nine different ecosystem functions
ran.sim %>%
  pivot_longer(cols = starts_with("F_"),
               names_to = "function.id",
               values_to = "function.value") %>%
  ggplot(data = .,
         mapping = aes(x = richness, y = function.value, colour = function.id)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "bottom")


### put data into format specified by Laura and Caroline

# split a.dat into a list
a.list <- split(a.dat, a.dat$sim.id)

# convert to species x site matrix for each simulation id
a.list <- 
  lapply(a.list, function(data){
    
    x <- 
      data %>%
      filter(abundance > 0) %>%
      arrange(species) %>%
      pivot_wider(id_cols = c("sim.id", "patch", "time"),
                  names_from = "species",
                  values_from = "abundance",
                  values_fill = 0) %>%
      arrange(patch)
    
    return(x)
    
  })

# test manipulation for a few random data points until confident

# select a random dataset
test.s <- sample(unique(a.dat$sim.id), 1)

# select a random patch
test.p <- sample(unique(a.dat$patch), 1)

y.test <- 
  a.dat %>%
  filter(sim.id == test.s) %>%
  filter(abundance > 0) %>%
  arrange(patch, species) %>%
  filter(patch == test.p) %>%
  pull(abundance)

z.test <- 
  a.list[[test.s]] %>%
  select(starts_with("sp_")) %>%
  slice(test.p) %>%
  unlist(.)

# any y.test not equal to z.test? if FALSE, then data are correct
any(y.test != z.test[z.test > 0])


### join the species abundance data to the function data

# convert the multifunctionality data into a list
mf.list <- split(mf.dat, mf.dat$sim.id)

# join these datasets together and output a list
mdi_list <- 
  mapply(function(x, y) { full_join(x, y, c("sim.id", "patch", "time")) }, 
         mf.list, a.list, SIMPLIFY = FALSE)

# view a single simulated dataset from this output list
View(mdi_list[[1]])


# save this list as an RDS file (I sent this file directly to both of you)

# saveRDS(mdi_list, here("data/multi_div_model_list.rds"))
# l_load <- readRDS(here("data/multi_div_model_list.rds")) 


# code to output each simulated dataset (1:30) as a separate .csv into the current working directory
# if not using here(), specify path manually in this loop
for(i in 1:length(mdi_list)) {
  
  file_name <- paste("simulation_", i, sep = "")
  write_csv(x = mdi_list[[i]], here("file_name"))
  
}

### END
