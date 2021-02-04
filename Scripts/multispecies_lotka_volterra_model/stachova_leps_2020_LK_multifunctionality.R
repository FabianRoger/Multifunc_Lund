
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Add multifunctionality to the abundance data from the stachova and leps (2010) model

# select scripts to call functions from
library(here)
source(here("Scripts/multispecies_lotka_volterra_model/stachova_leps_2010_LK_model_function.R"))
source(here("Scripts/MF_functions_collated.R"))

# equal competition
# stabilising competition
# strong competition

# crossed with:
# high specialisation
# low specialisation
# mixed i.e. some functional generalists and some some specialists


# simulate a set of communities
df.1 <- s_l_function(lsp = c(5, 10, 15, 20, 25),
                     reps = 20,
                     rsp = 30,
                     t_steps = 500,
                     n0 = 20,
                     ext.thresh = 0.2,
                     a_mean = 0.15, a_sd = 0.05, a_min = 0, a_max = 0.75, 
                     a_spp = 1, sim.comp = "sym",
                     k_min = 20, k_max = 150,
                     r_min = 0.1, r_max = 0.5)

# check basic relationships
library(ggplot2)
library(dplyr)
library(tidyr)

df.sum <- 
  df.1[[1]] %>%
  filter(total_abundance > 0)

ggplot(data = df.sum,
       mapping = aes(x = richness, y = total_abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

ggplot(data = df.sum,
       mapping = aes(x = local.sp.pool, y = total_abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# check the species parameters
df.1[[3]]


# adding function values using the Weibull distribution

# get the species abundance data
df.spp <- df.1[[2]]

df.spp %>% 
  filter(abundance > 0) %>%
  pull(species) %>%
  unique() %>%
  length()
  

# make a species list
spp.list <- unique(df.spp$species)

# choose the number of functions
func.n <- 9

# define the distribution parameters: Weibull distribution

w.shape <- 3 # 0.5, 3
w.scale <- 0.5 # 0.25, 0.5

# check the distribution
l <- rweibull(n = 1000, shape = w.shape, scale = w.scale)
hist(l)
mean(l)

prob.neg <- 0.2

func.mat <- 
  lapply(spp.list, function(x){
    x <- rweibull(n = func.n, shape = w.shape, scale = w.scale)
    y <- x - quantile(x, prob.neg)
    z <- (x*y)
    round(z, digits = 4)
})

func.mat <- data.frame(do.call(rbind, func.mat))
func.mat <- cbind(spp.list, func.mat)

names(func.mat) <- c("species", paste("F", 1:func.n, sep = "_"))

# quantify species specialisation
vegan::diversity((func.mat[, -1] + (min(func.mat[, -1])*-1) ), index = "shannon")
mean(vegan::vegdist(x = abs(func.mat[, -1]), method = "bray"))

func.mat %>%
  pivot_longer(cols = starts_with(match = "F_"),
               names_to = "functioning",
               values_to = "val") %>%
  ggplot(data = .,
         mapping = aes(x = species, y = val, colour = functioning)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

multi.func <- full_join(df.spp, func.mat, by = "species")

# multiply abundance by these function values
multi.func <- 
  multi.func %>%
  mutate(across(.cols = starts_with(match = "F_"), ~(.*abundance) ))

# function to standardise and then translate the functions
func.std.trans <- function(x){
  y <- (x - mean(x))/sd(x)
  return(y + abs(min(y)))
}

multi.func <- 
  multi.func %>%
  group_by(patch) %>%
  summarise(across(.cols = starts_with("F_"), sum), .groups = "drop") %>%
  mutate(across(.cols = starts_with("F_"), func.std.trans ))


# calculate the MF metrics on these datasets with different levels of correlation
func.names <- names(multi.func)
func.names <- func.names[grepl(pattern = "F_", func.names)]

x <- select(multi.func, func.names)

multi.func %>%
  mutate(`scal. MF` = MF_jing(adf = x, vars = func.names),
         `sum MF` = MF_sum(adf = x, vars = func.names),
         `ave. MF` = MF_av(adf = x, vars = func.names),
         `MESLI MF` = MF_mesli(adf = x, vars = func.names),
         `Pasari MF` = MF_pasari(adf = x, vars = func.names),
         `SAM MF` = MF_dooley(adf = x, vars = func.names),
         `ENF MF` = hill_multifunc(adf = x, vars = func.names, scale = 1, HILL = TRUE),
         `Simp. MF` = MF_simpsons_div(adf = x, vars = func.names),
         `Manning.30 MF` = manning_multifunc(adf = x, vars = func.names, thresh = 0.3),
         `Manning.50 MF` = manning_multifunc(adf = x, vars = func.names, thresh = 0.5),
         `Manning.70 MF` = manning_multifunc(adf = x, vars = func.names, thresh = 0.7),
         `thresh.30 MF` = single_threshold_mf(adf = x, vars = func.names, thresh = 0.3),
         `thresh.50 MF` = single_threshold_mf(adf = x, vars = func.names, thresh = 0.5),
         `thresh.70 MF` = single_threshold_mf(adf = x, vars = func.names, thresh = 0.7),
         `Slade.10.90 MF` = MF_slade(adf = x, vars = func.names, A_quant = 0.10, B_quant = 0.90),
         `Slade.40.60 MF` = MF_slade(adf = x, vars = func.names, A_quant = 0.40, B_quant = 0.60),
         `PCA MF` = pca_multifunc(adf = x, vars = func.names, standardise = FALSE) ) %>%
  View()
  

# MF metric groups
MF_groups <- 
  c(rep("sum/ave.", 4),
    rep("evenness", 4),
    rep("threshold", 8),
    rep("other", 1))














