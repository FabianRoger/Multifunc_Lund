
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code to compare different metrics of multifunctionality given different function distributions

# what are the crucial characteristics?

# should be zero if all functions are zero
# should be maximised if all functions are their maximum
# should be positive always because what does a negative multifunctionality mean?

# must be defined for any possible functional distribution (i.e. cannot be undefined)

# debatable characteristics

# should be zero if there is only one positive function


# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(readr)

rm(list = ls())

# tell R where to call scripts from
source(here("Scripts/function_plotting_theme.R"))
source(here("Scripts/MF_functions_collated.R"))

# load the jena data
jena.dat <- read_csv(file = here("data/jena_data_Jochum_2020_clean.csv"))

# check the column names
col_names <- names(jena.dat)

# remove the individual species abundances
spp <- col_names[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", col_names)]
length(spp)
rm(col_names)

jena.dat <- 
  jena.dat %>%
  select(-all_of(spp))

# fit distributions to these functions
hist(jena.dat$biomass)

library(fitdistrplus)
Fit.dist <- function(func, dists = c("norm", "lnorm", "gamma", "unif", "weibull")) {
  
  x <- dists
  y <- vector(length = length(x))
  u <- vector("list", length(x))
  for (i in 1:length(x)) {
    
    z <- fitdist(func, x[i])
    u[[i]] <- z
    z.sum <- summary(z)
    y[[i]] <- z.sum$aic
    
  }
  
  u[which(y == min(y))][[1]]
  
}

# choose six functions and find the best fitting distributions
names(jena.dat)
func.names <- c("biomass",  "plantCN", "soilC",             
                "soilorgC", "herbi", "micBMC", "phoact",            
                "poll", "rootBM")

best.dist <- vector("list", length = length(func.names))
for(i in 1:length(func.names)) {
  
  best.dist[[i]] <- Fit.dist(func = jena.dat[[func.names[i] ]])
  
}

# choose distributions that represent a range of possible distributions
best.dist.sub <- best.dist[c(1, 3, 5, 7, 8)]

# set up the number of points to simulate
n.sim <- 1000

# weibull distribution
best.dist.sub[[1]]
x1 <- rweibull(n = n.sim, shape = best.dist.sub[[1]]$estimate[1], scale = best.dist.sub[[1]]$estimate[2])

# uniform distribution
best.dist.sub[[2]]
x2 <- runif(n = n.sim, best.dist.sub[[2]]$estimate[1],best.dist.sub[[2]]$estimate[2] )

# lnorm
best.dist.sub[[3]]
x3 <- rlnorm(n = n.sim, best.dist.sub[[3]]$estimate[1], best.dist.sub[[3]]$estimate[2])

# weibull
best.dist.sub[[4]]
x4 <- rweibull(n = n.sim, best.dist.sub[[4]]$estimate[1], best.dist.sub[[4]]$estimate[2])

# gamma
best.dist.sub[[5]]
x5 <- rgamma(n = n.sim, best.dist.sub[[5]]$estimate[1], best.dist.sub[[5]]$estimate[2])

standardisation <- "max"

x.sim <- 
  lapply(list(x1, x2, x3, x4, x5), function(x) {
  
    y <- c(0, x)
    
  if (standardisation == "z_score") {
    
    y <- scale(y)[,1]
    y <- y + abs(min(y))
    
  }  else if (standardisation == "max_0_1") {
    
    y <- ( y - min(y) )/( max(y) - min(y) )
    
  } else if (standardisation == "max") {
    
    y <- y/max(y)
    
  }
    
  z <- range(y)
  
  seq(z[1], z[2], length.out = 10)
  
})

x.sims <- do.call(expand.grid, x.sim)
names(x.sims) <- paste("F_", 1:length(x.sim), sep = "")
sim.dat <- as_tibble(x.sims)
rm(x.sims)

# write a function to generate these function distributions

# use the empirical distribution data
func.names <- names(sim.dat)

sim.metric <- 
  sim.dat %>%
  mutate(`scal. MF` = MF_jing(adf = sim.dat, vars = func.names),
         `sum MF` = MF_sum(adf = sim.dat, vars = func.names),
         `ave. MF` = MF_av(adf = sim.dat, vars = func.names),
         `sd. MF` = apply(sim.dat[, func.names], 1, sd),
         `MESLI MF` = MF_mesli(adf = sim.dat, vars = func.names),
         `Pasari MF` = MF_pasari(adf = sim.dat, vars = func.names),
         `SAM MF` = MF_dooley(adf = sim.dat, vars = func.names),
         `ENF MF` = hill_multifunc(adf = sim.dat, vars = func.names, scale = 1, HILL = TRUE),
         `Simp. MF` = MF_simpsons_div(adf = sim.dat, vars = func.names),
         `Manning.30 MF` = manning_multifunc(adf = sim.dat, vars = func.names, thresh = 0.3),
         `Manning.50 MF` = manning_multifunc(adf = sim.dat, vars = func.names, thresh = 0.5),
         `Manning.70 MF` = manning_multifunc(adf = sim.dat, vars = func.names, thresh = 0.7),
         `thresh.30 MF` = single_threshold_mf(adf = sim.dat, vars = func.names, thresh = 0.3),
         `thresh.50 MF` = single_threshold_mf(adf = sim.dat, vars = func.names, thresh = 0.5),
         `thresh.70 MF` = single_threshold_mf(adf = sim.dat, vars = func.names, thresh = 0.7),
         `Slade.10.90 MF` = MF_slade(adf = sim.dat, vars = func.names, A_quant = 0.10, B_quant = 0.90),
         `Slade.40.60 MF` = MF_slade(adf = sim.dat, vars = func.names, A_quant = 0.40, B_quant = 0.60),
         `PCA MF` = pca_multifunc(adf = sim.dat, vars = func.names, standardise = FALSE),
         one_positive = apply(sim.dat, 1, function(x) { sum(sum(ifelse(x > 0, 1, 0)) == 1) }),
         n_positive = apply(sim.dat, 1, function(x) { sum(ifelse(x > 0, 1, 0)) })
         )

rm(sim.dat)

# add a row.id column
sim.metric$row_id <- 1:nrow(sim.metric)

# what do we need?

# 1. all values at max
# 2. all values at zero
# 3. all values that are undefined (i.e. Inf, -Inf or NA)
# 4. the rest

# create the simulated dataset to plot
df <- sim.metric[sample(1:nrow(sim.metric), 1000), ]
df <- 
  df %>%
  mutate(fig.1.group = c(NA),
         group.colours = 1)

df.max <- 
  sim.metric %>%
  filter(across(.cols = starts_with("F_"), ~.x == max(.x, na.rm = TRUE)))

df.max <- 
  df.max %>%
  mutate(fig.1.group = c("a"),
         group.colours = 2)

df.min <- 
  sim.metric %>%
  filter(across(.cols = starts_with("F_"), ~.x == min(.x, na.rm = TRUE)))

df.min <- 
  df.min %>%
  mutate(fig.1.group = c("b"),
         group.colours = 3)

# remove duplicates in the df data.frame
df <- 
  df %>%
  filter(!(row_id %in% c(df.max$row_id, df.min$row_id)))

df.all <- 
  rbind(df.min, df, df.max)

# 



# get the median average MF and median

# make a test plot
library(ggrepel)
ggplot( mapping = aes(x = `ave. MF`, y = `sd. MF` ) ) +
  geom_point(data = df, alpha = 0.1) +
  geom_point(data = df.max, colour = "red", size = 2, shape = 17) +
  geom_point(data = df.min, colour = "green", size = 3, shape = 18) +
  geom_label_repel(data = df.label, mapping = aes(label = fig.1.group ),
                   segment.colour = "black", min.segment.length = 0.2,
                   fill = "grey", colour = "white") +
  ylab("SD among functions") +
  xlab("Mean among functions") +
  scale_colour_viridis_c() + 
  theme_meta()

# plot each highlighted point in the test plot
df.label.long <- 
  df.label %>%
  dplyr::select(fig.1.group, starts_with("F_") ) %>%
  pivot_longer(cols = starts_with("F_"),
               names_to = "function_id", 
               values_to = "function_value")

ggplot(data = df.label.long %>% filter(fig.1.group == "f"),
       mapping = aes(x = function_id, y = function_value)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2)) +
  theme_meta()


# for each metric, what do we want to plot separately?
# set the metric
names(df.all)
metric <- "Simp. MF"

df.metric <- 
  df.all %>%
  filter(!is.na(get(metric))) %>%
  filter(get(metric) != Inf) %>%
  filter(get(metric) != -Inf) %>%
  mutate(metric_negative = as.character(if_else(get(metric) < 0, 1, 0)) )

df.metric.max <- 
  df.metric %>%
  filter(get(metric) != Inf) %>%
  filter(get(metric) == max(get(metric)))

df.metric.min <- 
  df.metric %>%
  filter(get(metric) != -Inf) %>%
  filter(get(metric) == min(get(metric)))

df.metric.undefined <- 
  df.all %>%
  filter(is.na(get(metric)) | get(metric) == Inf | get(metric) == -Inf )

ggplot() +
  geom_point(data = df.metric,
             mapping = aes(x = `ave. MF`, y = `sd. MF`, shape = metric_negative, colour = get(metric) ),
             size = 3, alpha = 0.5, position = position_jitter(width = 0.005)) +
  scale_colour_viridis_c(option = "D") +
  scale_shape_manual(values=c(16, 15), guide = FALSE)+
  geom_hline(yintercept = df.min[["sd. MF"]], linetype = "dashed", colour = "black") +
  geom_vline(xintercept = df.max[["ave. MF"]], linetype = "dashed", colour = "black") +
  geom_point(data = df.min, 
             mapping = aes(x = `ave. MF`, y = `sd. MF`), 
             fill = "white", colour = "black", shape = 24, alpha = 1, size = 3.5, 
             position = position_nudge(x = -0.025), stroke = 1.25 ) +
  geom_point(data = df.max,
             mapping = aes(x = `ave. MF`, y = `sd. MF`),
             fill = "white", colour = "black", shape = 21, alpha = 1, size = 3.5,
             position = position_nudge(x = -0.025), stroke = 1.25 ) +
  geom_point(data = df.metric.max,
             mapping = aes(x = `ave. MF`, y = `sd. MF`),
             colour = "black", fill = "#FF6600", shape = 21, alpha = 1, size = 3.5, stroke = 1.25) +
  geom_point(data = df.metric.min,
             mapping = aes(x = `ave. MF`, y = `sd. MF`),
             colour = "black", fill =  "#FF6600", shape = 24, alpha = 1, size = 3.5, stroke = 1.25) +
  geom_point(data = df.metric.undefined,
             mapping = aes(x = `ave. MF`, y = `sd. MF`), 
             shape = 4, size = 3, position = position_nudge(x = 0.025) ) +
  ylab("SD among functions") +
  xlab("Mean among functions") +
  ggtitle(metric) +
  # scale_y_continuous(limits = c(-0.05, 0.55)) +
  # scale_x_continuous(limits = c(-0.05, 1.05)) +
  theme_meta() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 9)) +
  guides( colour = guide_colourbar(barwidth = 0.5, barheight = 10,
                                   frame.colour = "black", ticks = FALSE ) )



# not all the metrics are maximised when all the functions are at their maximum
# for a given set of plots

# SAM metric
# PCA metric
# Simpson index

sim.metric %>%
  filter(F_1 == max(F_1), 
         F_2 == max(F_2),
         F_3 == max(F_3),
         F_4 == max(F_4),
         F_5 == max(F_5)) %>%
  View()

lapply(sim.metric, range, na.rm = TRUE)

cor(sim.metric[sample(1:nrow(sim.metric), 10000), -c(1, 2, 3, 4, 5)])


x <- rnorm(n = 1000, 100, 10)
y <- runif(n = 1000, min = 0, max = 10)
z <- rpois(n = 1000, 50)

lapply(list(x, y, z), function(x){
  
  scale(x)[,1] %>% range()
  
})



# 0 - 1 sims

univariate_explorer <- function(funcnum = 5, grain = 0.1, error = 0.01) {
  
  # generate n functions between 0 and 1
  x <- replicate(funcnum, seq(0, 1, grain), simplify = FALSE)
  
  # put these functions into a matrix
  y <- do.call("expand.grid", x)
  
  # add a small amount of error and set to zero if zero
  nm <- nrow(y)*ncol(y)
  y <- y + rnorm(n = nm, 0, error)
  y[y < 0] <- 0
  
  names(y) <- paste("F_", 1:funcnum, sep = "")
  
  # convert y to a tibble
  df <- tibble(y)
  
  # only get unique rows
  df <- distinct(df)
  
  # output the df
  df
  
}

# simulate data for five ecosystem functions using univariate simulator
# set.seed(43767)
# sim.dat <- univariate_explorer() 
# head(sim.dat)

# get a vector of names
# func.names <- names(sim.dat)



