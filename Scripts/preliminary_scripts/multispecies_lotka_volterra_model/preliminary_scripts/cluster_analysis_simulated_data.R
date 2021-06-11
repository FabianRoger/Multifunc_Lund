
# Project: Mutifunctionality workshop (Lund 2019)

# Title: Cluster analysis of multifunctionality metrics based on simulated data

# load relevant libraries
library(dplyr)
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

# add this variable to the mult.dat data
mult.dat <- 
  full_join(select(params, rep.id, sim.group, sim.id), 
            mult.dat,
            by = "sim.id")

# examine the relationship between richness and average MF for each simulation
d1 <- 
  ggplot(data = mult.dat,
       mapping = aes(x = richness, y = ave._MF, colour = as.character(sim.id) )) +
  geom_jitter(width = 0.1, alpha = 0.1, shape = 16) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  ylab("average MF") +
  xlab("sp. richness") +
  theme_meta() +
  theme(legend.position = "none")

d2 <- 
  ggplot(data = mult.dat,
         mapping = aes(x = richness, y = Pasari_MF, colour = as.character(sim.id) )) +
  geom_jitter(width = 0.1, alpha = 0.1, shape = 16) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  ylab("Pasari MF") +
  xlab("sp. richness") +
  theme_meta() +
  theme(legend.position = "none")

d3 <- 
  ggplot(data = mult.dat,
         mapping = aes(x = richness, y = thresh.50_MF, colour = as.character(sim.id) )) +
  geom_jitter(width = 0.1, alpha = 0.1, shape = 16) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  ylab("threshold MF (70%)") +
  xlab("sp. richness") +
  theme_meta() +
  theme(legend.position = "none")


# plot a correlation matrix
library(ellipse)

# calculate a correlation matrix
mult.cor <- 
  mult.dat %>%
  select(ends_with("_MF")) %>%
  cor(.)

# order the correlation matrix
plotcorr(mult.cor , mar=c(1,1,1,1), type = c("lower"), yaxt="none")
d4 <- corrplot::corrplot(corr = mult.cor, method = c("ellipse"),
                   type = "lower",
                   tl.pos = "d", tl.col = "black", tl.cex = 0.75,
                   pos = 4, cl.align.text = "c",
                   mar = c(0, 0, 0, 0))


# perform the cluster analysis
library(vegan)

# MF metric groups
MF_groups <- 
  c(rep("sum/ave.", 3),
    rep("evenness", 5),
    rep("threshold", 8),
    rep("other", 1))

# Run an nMDS to cluster the different multifunctionality metrics

# split the data.frame into simulation types and append this to the full dataset
dat.list <- c(list(mult.dat), split(mult.dat, mult.dat$sim.group))

nmds.figs <- 
  lapply(dat.list, function(data){
  
  mds_out <- 
    data %>%
    select(ends_with("MF")) %>% 
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
                              size = 2.25,
                              segment.alpha = 0.5,
                              label.size = NA, fill = "white") +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    guides(label = FALSE,
           colour = guide_legend(override.aes = list(size = 3))) +
    theme_meta() +
    theme(legend.position = "right",
          legend.key = element_blank(),
          legend.title = element_blank())
  
} )

# plot the full data.set
nmds.figs[[1]]

# combine this plot with example plots using patchwork
library(patchwork)
e.g <- (d1 + d2 + d3)

ggsave(filename = here("Figures/pca_clust_fig.png"), plot = nmds.figs[[1]],
       width = 13, height = 8, units = "cm", dpi = 450)

# ggsave(filename = here("Figures/pca_clust_examples.png"), plot = e.g,
       # width = 17, height = 8, units = "cm", dpi = 450)


# combine the remaining plots using patchwork
library(patchwork)

comb.nmds <- 
  nmds.figs[[2]] + nmds.figs[[3]] + nmds.figs[[4]] +
  nmds.figs[[5]] + nmds.figs[[6]] + nmds.figs[[7]] +
  plot_layout(guides = 'collect') &
  theme(legend.position='bottom') &
  theme(axis.title.x = element_text(colour ="black", size = 9.5, face = "plain", margin=margin(5,0,0,0,"pt")),
        axis.title.y = element_text(colour = "black", size = 9.5, face = "plain", margin=margin(0,5,0,0,"pt")),
        axis.text.x = element_text(colour = "black", size=8, face = "plain",  margin=margin(10,0,0,0,"pt")),
        axis.text.y = element_text(colour ="black", size=8, face = "plain", margin=margin(0,10,0,0,"pt")),
        )
  
ggsave(filename = here("Figures/pca_clust_fig_all.png"), plot = comb.nmds,
       width = 20, height = 19, units = "cm", dpi = 450)

# combine these nmds plots
# p1 <- 
  # ggarrange(plotlist = nmds.figs, 
            # labels = letters[1:length(nmds.figs)],
            # font.label = list(size = 12, color = "black", face = "plain", family = NULL),
            # common.legend = TRUE,
            # legend = "bottom")

### END