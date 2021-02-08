
# Project: Mutifunctionality workshop (Lund 2019)

# Title: Cluster analysis of multifunctionality metrics based on simulated data

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)

# define scripts to pull functions from
source(here("Scripts/function_plotting_theme.R"))

# read in the simulated data
params <- read_csv(file = here("data/parameters_sim.csv"))
head(params)

mult.dat <- read_csv(file = here("data/multifunctionality_data.csv"))
head(mult.dat)

# create a simulation category variable
# reps <- length(unique(params$rep.id))
# id <- unique(paste(params$a_mean, params$w.shape, params$w.scale, sep = "_"))
# id <- LETTERS[1:length(id)]

# add this to the params data
# params$sim.group <- rep(id, each = reps)

# add this variable to the mult.dat data
mult.dat <- 
  full_join(select(params, rep.id, sim.group, sim.id), 
            mult.dat,
            by = "sim.id")

# examine the relationship between richness and average MF for each simulation
ggplot(data = mult.dat,
       mapping = aes(x = richness, y = `ave. MF`, colour = as.character(sim.id) )) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "none")

ggplot(data = mult.dat,
       mapping = aes(x = richness, y = `thresh.70 MF`, colour = as.character(sim.id) )) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~sim.id, scales = "free") +
  theme_classic() +
  theme(legend.position = "none")


### perform the cluster analysis
library(vegan)

# MF metric groups
MF_groups <- 
  c(rep("sum/ave.", 4),
    rep("evenness", 4),
    rep("threshold", 8),
    rep("other", 1))

# Run an nMDS to cluster the different multifunctionality metrics
nmds.figs <- 
  lapply(split(mult.dat, mult.dat$sim.group), function(data){
  
  mds_out <- 
    data %>%
    select(ends_with(" MF")) %>% 
    mutate(across(.cols = everything(), ~scale(.) ) ) %>% 
    as.matrix(.) %>% 
    t() %>% 
    dist(x = ., method = "euclidean") %>% 
    metaMDS(comm = .) 
  
  # check the nmds stress and fit
  stressplot(mds_out)
  print(mds_out$stress)
  
  mds_out$points %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(var = "MF_metrics") %>% 
    mutate(group = factor(as.factor(MF_groups), 
                          levels = c("sum/ave.", "evenness", "threshold", "other")) ) %>%
    ggplot(data = .,
           mapping = aes(x = MDS1, y = MDS2, colour = group)) +
    geom_point() +
    ggrepel::geom_label_repel(mapping = aes(label = MF_metrics),
                              show.legend = FALSE,
                              size = 3,
                              segment.alpha = 0.5,
                              label.size = NA, fill = "white") +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    guides(label = FALSE,
           colour = guide_legend(override.aes = list(size = 3))) +
    theme_meta() +
    theme(legend.position = "bottom",
          legend.key = element_blank(),
          legend.text = element_text(colour = "black", size = 14, face = "plain"),
          legend.title = element_blank())
  
} )


# combine these nmds plots
p1 <- 
  ggarrange(plotlist = nmds.figs, 
            labels = letters[1:length(nmds.figs)],
            font.label = list(size = 12, color = "black", face = "plain", family = NULL),
            common.legend = TRUE,
            legend = "bottom")

    
ggsave(filename = here("Figures/pca_clust_fig.png"), plot = p1,
       width = 19, height = 20, units = "cm", dpi = 300)

















