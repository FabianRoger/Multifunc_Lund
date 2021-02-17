
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Characterise the simulated data and Jena data by BEF relationships, abundance distributions etc.

# load relevant libraries
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# scripts to draw functions from
source(here("Scripts/function_plotting_theme.R"))

# load the simulated data

# load the simulation parameters
par.dat <- read_csv(file = here("data/parameters_sim.csv"))
head(par.dat)

# load the raw multifunctionality data
sim.raw <- read_csv(file = here("data/multifunctionality_data.csv"))
head(sim.raw)

# load the raw abundance data
abun.raw <- read_csv(file = here("data/raw_abundance_data.csv"))
head(abun.raw)


# choose the sim.id's to plot
# choose six random patches (one from each simulation type)
sim.id.plots <- 
  par.dat %>%
  group_by(sim.group) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(sim.id)

# sim.id
sim.number <- sim.id.plots[4]
params[sim.number,]

# modify the plotting theme to make the text smaller
theme_meta2 <- function() {
  theme_meta() +
    theme(
      axis.title.x = element_text(colour ="black", size = 9.5, face = "plain", margin=margin(5,0,0,0,"pt")),
      axis.title.y = element_text(colour = "black", size = 9.5, face = "plain", margin=margin(0,5,0,0,"pt")),
      axis.text.x = element_text(colour = "black", size=8, face = "plain",  margin=margin(10,0,0,0,"pt")),
      axis.text.y = element_text(colour ="black", size=8, face = "plain", margin=margin(0,10,0,0,"pt"))
    )
}

plot.list <- vector("list", length = length(sim.id.plots))
for (i in 1:length(sim.id.plots)) {
  
  # plot rank-abundance distributions for four randomly chosen patches
  p.df <- 
    abun.raw %>%
    filter(sim.id == sim.id.plots[i]) %>%
    filter(abundance > 0) %>%
    mutate(patch = as.character(patch))
  
  # get four random patches
  ran.patch <- sample(x = unique(p.df$patch), size = 4)
  
  p.df.ran <- 
    p.df %>%
    filter(patch %in% ran.patch) %>%
    group_by(patch) %>%
    mutate(abundance_rank = rank(-abundance)) %>%
    ungroup()
  
  p1 <- 
    ggplot(data = p.df.ran,
           mapping = aes(x = (abundance_rank), y = abundance, colour = patch)) +
    stat_smooth(method='lm', formula = y~poly(x, 2), se = FALSE, size = 1) +
    geom_jitter(width = 0.1) +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    scale_x_continuous(breaks = c(0:max(p.df.ran$abundance_rank))) +
    ylab("abundance") +
    xlab("rank") +
    theme_meta2() +
    theme(legend.position = "none") # +
  # facet_wrap(~patch, scales = "free") +
  # theme(legend.position = "bottom",
  # legend.key = element_rect(fill = NA),
  # strip.background =element_rect(fill = "white", colour = "black"))
  
  
  # plot richness-total abundance relationship
  mf.df <- 
    sim.raw %>%
    filter(sim.id == sim.id.plots[i])
  
  p2 <- 
    ggplot(data = mf.df,
           mapping = aes(x = richness, y = total_abundance)) +
    geom_smooth(method = "lm", se = TRUE, colour = "black", alpha = 0.2, size = 0.75) +
    geom_jitter(width = 0.1) +
    xlab("sp. richness") +
    ylab("tot. abundance") +
    theme_meta2() # +
  # theme(plot.margin = unit(c(5.5,5.5,57,5.5), "pt"))
  
  
  # plot richness-function relationships for the nine functions
  p3 <- 
    mf.df %>%
    pivot_longer(cols = starts_with("F_"),
                 names_to = "EF",
                 values_to = "value") %>%
    ggplot(data = .,
           mapping = aes(x = richness, y = value, colour = EF)) +
    geom_jitter(width = 0.1, alpha = 0.5, shape = 16) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    ylab("function") +
    xlab("sp. richness") +
    theme_meta2() +
    theme(legend.position = "none") # +
  # theme(legend.position = "bottom",
  # legend.key = element_rect(fill = NA),
  # strip.background =element_rect(fill = "white", colour = "black"),
  # legend.text = element_text(size = 9))
  
  # combine packages using patchwork
  
  # install.packages("devtools")
  # devtools::install_github("thomasp85/patchwork")
  library(patchwork)
  
  # combine the plots
  p.c <- (p1 + p2 + p3)
  
  # add the combined plot into a list
  plot.list[[i]] <- p.c
  
}

# high functional specialisation
par.dat[sim.id.plots,]
plot.list[[1]]/plot.list[[2]]/plot.list[[3]]

# low functional specialisation
par.dat[sim.id.plots,]
plot.list[[4]]/plot.list[[5]]/plot.list[[6]]






