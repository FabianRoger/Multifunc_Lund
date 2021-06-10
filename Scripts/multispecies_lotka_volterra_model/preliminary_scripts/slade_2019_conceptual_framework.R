
# Slade et al. (2019, Trends in Plant Science) framework

# select scripts to call functions from
library(here)
source(here("Scripts/multispecies_lotka_volterra_model/stachova_leps_2010_LK_model_function.R"))
source(here("Scripts/MF_functions_collated.R"))

# libraries for data manipulation
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# fixed parameters

# lotka-volterra model
lsp = c(1, 2)
reps = 2
rsp = 2
t_steps = 500
n0 = 100
a_min = 0
a_max = 1
sim.comp = "sym"
a_scale = 1


# simulate a set of communities
df.x <- 
  s_l_function(lsp = lsp,
               reps = reps,
               rsp = rsp,
               t_steps = t_steps,
               n0 = n0,
               a_mean = 0.25, a_sd = 0, 
               a_min = a_min, a_max = a_max,
               a_spp = 0.75, sim.comp = sim.comp,
               k_min = 100, k_max = 200,
               r_min = 0.25, r_max = 0.5
  )

df.x$spp.info
df.x$data.summary
spp.list <- unique(df.x$data.raw$species)

func.n <- 2

# make a function name list
func.names <- paste("F_", 1:func.n, sep = "")

# generate a function matrix for the relationship between each species abundance and each function
func.mat <- 
  lapply(spp.list, function(x){
    x <- rweibull(n = 2, shape = 1, scale = 2)
    x[x > 10] <- rnorm(n = sum(x > 10), mean = 10, sd = 2) # constrain values to around 10
    y <- x - quantile(x, 0.1)
    z <- round(y, digits = 4)
  })

func.mat

func.mat <- data.frame(do.call(rbind, func.mat))
func.mat <- cbind(spp.list, func.mat)
func.mat
names(func.mat) <- c("species", func.names)

func.mat

df.sub <- 
  full_join(select(df.x$data.summary, patch, local.sp.pool),
            df.x$data.raw %>%
              pivot_wider(names_from = "species",
                          values_from = "abundance"), by = "patch") %>%
  mutate(local.sp.pool = as.character(local.sp.pool))

# plot the species abundance outcomes
ggplot(data = df.sub,
         mapping = aes(x = sp_1, y = sp_2, colour = local.sp.pool)) +
  geom_point(size = 2) +
  geom_line(data = filter(df.sub, local.sp.pool == 1),
            mapping = aes(x = sp_1, y = sp_2)) +
  theme_classic()



# get function values for species 1 and species 2
plot(df.x$data.summary$richness, df.x$data.summary$total_abundance)








