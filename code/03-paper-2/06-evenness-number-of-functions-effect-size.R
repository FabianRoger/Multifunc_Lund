#'
#' @title Why does the slope change with the number of functions for evenness
#' 
#' @description Simulate the case of considering more functions to illustrate
#' that the evenness among functions 
#'

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load relevant scripts
source("code/helper-plotting-theme.R")
source("code/03-paper-2/06-helper-gamfeldt-roger-functions.R")

# set a seed for reproducibility
set.seed(54807)

# set the number of species and the number of functions
specnum <- 15
funcnum <- 9

# specify the distribution used to draw species' traits
distribution = "runif"

# get the species function matrix
FuncMat <- FunctionValue(specnum, funcnum, distribution, min = 0, max = 1)

# get function name and species name vectors
func_names <- as.character(unique(FuncMat$Functions))
spec_names <- as.character(unique(FuncMat$Species))

# set the maximum number of replicates of N functions and species combinations
maxrep <- 50

# set up the species combination matrix with the maximum replications
SpecMat <- SpeciesMatrix(specnum = specnum, maxrep = maxrep)

# set the method to calculate the function values in a plot based on species traits
method = "av"

# get the values of different functions for each species combination
AvFunc <- AverageFunction(SpecMat, FuncMat,
                          method = method,
                          CF = CF, 
                          compfunc = compfunc,
                          r = r)

# check the summary of the different variables
summary(AvFunc)

# empty dataframe to store results
RES_func <- tibble(div = numeric(),
                   MF_eff_Q1 = numeric(), 
                   MF_eff_Q2 = numeric(), 
                   EffN_Q1 = numeric(),
                   EffN_Q2 = numeric(), 
                   AV = numeric(),
                   nfunc = numeric(),
                   func_comb = numeric())

# loop over chosen subsets of all function of varying size
for (i in 2:9) { 
  
  # all possible combination of i out of funcnum functions
  func_comb <- combn(func_names, i)
  
  # sample 50 random function combinations if more than 50 possible combinations
  if(ncol(func_comb) > 50) {
    func_comb <- func_comb[, sample(c(1:ncol(func_comb)), 50)]
  }
  
  # loop over all function combinations of size i
  for ( k  in seq_len(ncol(func_comb))) { 
    
    # effective number of functions of q = 1 EMF metric
    MF_eff_Q1 <- multifunc::getMF_eff(AvFunc, func_comb[ ,k], q = 1)
    
    # effective number of functions of q = 1 EMF metric
    MF_eff_Q2 <- multifunc::getMF_eff(AvFunc, func_comb[ ,k], q = 2)
    
    # get effective number of functions q = 1
    EffN_Q1 <- multifunc::eff_num_func(AvFunc, func_comb[ ,k], q = 1, standardized = FALSE)
    
    # get effective number of functions q = 2
    EffN_Q2 <- multifunc::eff_num_func(AvFunc, func_comb[ ,k], q = 2, standardized = FALSE)
    
    # get average MF
    AV <- multifunc::getStdAndMeanFunctions(AvFunc, func_comb[ ,k])
    
    temp <- tibble(div = AvFunc$Richness,
                   MF_eff_Q1 = MF_eff_Q1, 
                   MF_eff_Q2 = MF_eff_Q2, 
                   EffN_Q1 = EffN_Q1,
                   EffN_Q2 = EffN_Q2, 
                   AV = AV$meanFunction,
                   nfunc = i,
                   func_comb = k)
    
    RES_func <- rbind(RES_func, temp)
    
  }
  
}

# check the results
head(RES_func)

RES_func %>% 
  mutate(EffN_Q1_st = EffN_Q1 / nfunc) %>% 
  pivot_longer(one_of(c("MF_eff_Q1", "EffN_Q1", "AV"))) %>% 
  ggplot(aes(x = div, y = value))+
  geom_point(size = 0.1, colour = "grey")+
  geom_smooth(method = "lm", se = F, size = 0.8, colour = "red")+
  # geom_smooth(se = F, size = 0.8, colour = "blue")+
  facet_grid(name~nfunc, scales = "free_y")+
  theme_bw()

# standardise the effective number of functions
RES_func <- 
  RES_func |>
  dplyr::mutate(EffN_Q1_st = EffN_Q1 / nfunc,
                EffN_Q2_st = EffN_Q2 / nfunc)
  
# plot the raw relationships  
RES_func |>
  ggplot(aes(x = div, y = AV, group = func_comb)) +
  geom_point(shape = 1, alpha = 0.1, colour = "grey") +
  geom_smooth(method = "lm", se = F, size = 0.1, colour = "red")+
  facet_wrap(~nfunc, nrow = 1, ncol = 8) +
  theme_meta()

RES_func |>
  tidyr::pivot_longer(one_of(c("MF_eff_Q2", "EffN_Q2", "AV"))) |>
  ggplot(aes(x = div, y = value, group = func_comb)) +
  geom_smooth(method = "lm", se = F, size = 0.1, colour = "red")+
  facet_grid(name~nfunc, scales = "free_y")+
  theme_meta()

# plot the range of evenness
RES_func |>
  dplyr::group_by(nfunc, func_comb) |>
  dplyr::summarise(min_EffN_Q1 = min(EffN_Q1),
                   max_EffN_Q1 = max(EffN_Q1), .groups = "drop") |>
  ggplot() +
  geom_segment(mapping = aes(x = nfunc, xend = nfunc, y = min_EffN_Q1, yend = max_EffN_Q1,
                             group = func_comb),
             position = position_jitter(width = 0.3), alpha = 0.1) +
  ylab("EFN range") +
  scale_x_continuous(breaks = c(1:9)) +
  theme_meta()

RES_func |>
  dplyr::group_by(nfunc, func_comb) |>
  dplyr::summarise(min_EffN_Q1 = min(EffN_Q1),
                   max_EffN_Q1 = max(EffN_Q1), .groups = "drop") |>
  ggplot(mapping = aes(x = nfunc, y = (max_EffN_Q1 - min_EffN_Q1 ))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("EFN range") +
  scale_x_continuous(breaks = c(1:9)) +
  theme_meta()
  
# plot the range of the average
RES_func |>
  dplyr::group_by(nfunc, func_comb) |>
  dplyr::summarise(min_AV = min(AV),
                   max_AV = max(AV), .groups = "drop") |>
  ggplot() +
  geom_segment(mapping = aes(x = nfunc, xend = nfunc, y = min_AV, yend = max_AV,
                             group = func_comb),
               position = position_jitter(width = 0.3), alpha = 0.1) +
  ylab("Average range") +
  scale_x_continuous(breaks = c(1:9)) +
  theme_meta()

RES_func |>
  dplyr::group_by(nfunc, func_comb) |>
  dplyr::summarise(min_AV = min(AV),
                   max_AV = max(AV), .groups = "drop") |>
  ggplot(mapping = aes(x = nfunc, y = (max_AV - min_AV ))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Average range") +
  scale_x_continuous(breaks = c(1:9)) +
  theme_meta()

# plot the range of the ENF-Q1 EMF
RES_func |>
  dplyr::group_by(nfunc, func_comb) |>
  dplyr::summarise(min_MF_Q1 = min(MF_eff_Q1),
                   max_MF_Q1 = max(MF_eff_Q1), .groups = "drop") |>
  ggplot() +
  geom_segment(mapping = aes(x = nfunc, xend = nfunc, y = min_MF_Q1, yend = max_MF_Q1,
                             group = func_comb),
               position = position_jitter(width = 0.3), alpha = 0.1) +
  ylab("ENF-Q1 EMF range") +
  scale_x_continuous(breaks = c(1:9)) +
  theme_meta()

### END
