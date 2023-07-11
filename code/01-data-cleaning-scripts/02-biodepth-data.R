#'
#' @title Download and clean the BIODEPTH data
#' 
#' @description Load the BIODEPTH data taken from Byrnes et al. (2014, Methods 
#' in Ecology and Evolution), clean it and output a cleaned version for 
#' further analyses
#'

# load relevant libraries
library(dplyr)
library(multifunc)

# get the BIODEPTH data
data("all_biodepth")

# check the downloaded data
head(all_biodepth)
summary(all_biodepth)

# make a vector of the relevant function data
all_vars <- c("biomassY3", "root3", "N.g.m2",  "light3", "N.Soil", "wood3", "cotton3")

# make an id variable with the function names
var_id <- which(names(all_biodepth) %in% all_vars)

# check the possible locations
unique(all_biodepth$location)

# subset out Sweden
swe_dat <- 
  all_biodepth |>
  dplyr::filter(location == "Sweden")

# which variables have > 2/3 of the values not NA?
swe_vars <- which(names(swe_dat) %in% multifunc::whichVars(swe_dat, all_vars, thresh = 0))

# What are the names of species in this dataset
# that have at least some values > 0?
swe_sp <- multifunc::relevantSp(swe_dat, 26:ncol(swe_dat))

# get the column ids with species that have some data that are not zero
sp_id <- which(names(swe_dat) %in% swe_sp)

# get the relevant columns
swe_dat <- swe_dat[,c(1:14, sp_id, swe_vars)]

# write this into a .rds file
saveRDS(object = swe_dat, file = "data/biodepth_swe_data.rds")


# subset out Portugal
prt_dat <- 
  all_biodepth |>
  dplyr::filter(location == "Portugal")

# which variables have > 2/3 of the values not NA?
prt_vars <- which(names(prt_dat) %in% multifunc::whichVars(prt_dat, all_vars, thresh = 0))

# What are the names of species in this dataset
# that have at least some values > 0?
prt_sp <- multifunc::relevantSp(prt_dat, 26:ncol(prt_dat))

# get the column ids with species that have some data that are not zero
sp_id <- which(names(prt_dat) %in% prt_sp)

# get the relevant columns
prt_dat <- prt_dat[,c(1:14, sp_id, prt_vars)]

# write this into a .rds file
saveRDS(object = prt_dat, file = "data/biodepth_prt_data.rds")

### END
