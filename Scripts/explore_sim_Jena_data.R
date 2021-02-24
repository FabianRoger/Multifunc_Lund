
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
source(here("Scripts/MF_functions_collated.R"))


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

# modify the plotting theme to make the text smaller
theme_meta2 <- function() {
  theme_meta() +
    theme(
      axis.title.x = element_text(colour ="black", size = 9.5, face = "plain", margin=margin(5,0,0,0,"pt")),
      axis.title.y = element_text(colour = "black", size = 9.5, face = "plain", margin=margin(0,5,0,0,"pt")),
      axis.text.x = element_text(colour = "black", size=8, face = "plain",  margin=margin(10,0,0,0,"pt")),
      axis.text.y = element_text(colour ="black", size=8, face = "plain", margin=margin(0,10,0,0,"pt")),
      title = element_text(colour = "black", size = 8, face = "plain")
    )
}



# set-up plot titles
plot.titles <- 
  c("weak comp. + specialisation",
  "strong comp. + specialisation",
  "equal comp. + specialisation",
  "weak comp. + generalism",
  "strong comp. + generalism",
  "equal comp. + generalism")


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
    ggtitle(plot.titles[i]) +
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
    ggtitle("") +
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
    geom_smooth(method = "lm", se = FALSE, size = 0.75) +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    ylab("function") +
    xlab("sp. richness") +
    ggtitle("") +
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
h.s <- plot.list[[1]]/plot.list[[2]]/plot.list[[3]]
h.s

ggsave(filename = here("Figures/sim.exp1.png"), plot = h.s,
       width = 17.3, height = 16, units = "cm", dpi = 450)

# low functional specialisation
par.dat[sim.id.plots,]
l.s <- (plot.list[[4]]/plot.list[[5]]/plot.list[[6]])
l.s

ggsave(filename = here("Figures/sim.exp2.png"), plot = l.s,
       width = 17.3, height = 16, units = "cm", dpi = 450)


# create the same plots for the Jena data

# load the cleaned Jena data
jena.raw <- read_csv(file = here("data/jena_data_cleaned.csv"))

# define variable groups
var.names <- names(jena.raw)

# (1) get species names
spp.p <- ( grepl("+\\.+", var.names) & nchar(var.names) == 7 )
spp.names <- var.names[spp.p]

# (2) get site identifiers
site.id <- c("year", "sowndiv", "plotcode", "realised_diversity")

# (3) get function names
jena.func.names <- var.names[!(var.names %in% ( c(spp.names, site.id) ))  ]

# subset the site.id and function names
jena.dat <- 
  jena.raw %>%
  select(all_of( c(site.id, jena.func.names) ))

# standardise the different ecosystem functions
jena.dat <- 
  jena.dat %>%
  mutate(across(.cols = all_of(jena.func.names), standardise))

# plot richness versus community biomass
j1 <- 
  ggplot(data = jena.dat,
       mapping = aes(x = realised_diversity, y = BM_targ_DW)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", alpha = 0.2, size = 0.75) +
  geom_jitter(width = 0.1) +
  xlab("sp. richness") +
  ylab("community biomass") +
  theme_meta2() +
  theme(plot.margin = unit(c(5.5,5.5,57,5.5), "pt"))


# plot richness versus the other functions
j2 <- 
  jena.dat %>%
  pivot_longer(cols = jena.func.names[jena.func.names != "BM_targ_DW"],
               names_to = "EF",
               values_to = "value") %>%
  ggplot(data = .,
         mapping = aes(x = realised_diversity, y = value, colour = EF)) +
  geom_jitter(width = 0.1, alpha = 0.5, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, size = 0.75) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  ylab("function") +
  xlab("sp. richness") +
  theme_meta2() +
  theme(legend.position = "bottom",
        legend.key = element_rect(fill = NA),
        strip.background =element_rect(fill = "white", colour = "black"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7) )

library(patchwork)
dev.off()
j.c <- (j1 + j2) + plot_layout(guides = 'collect') &
  theme(legend.position='bottom') &
  theme(plot.margin = unit(c(5.5,7.5,5.5,7.5), "pt"))
j.c

ggsave(filename = here("Figures/jena.exp.png"), plot = j.c,
       width = 11, height = 7.5, units = "cm", dpi = 450)

### END

