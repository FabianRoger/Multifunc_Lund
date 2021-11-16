
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
library(vegan)

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
  dplyr::select(-all_of(spp))

# wrap this into a function
calculate_MF <- function(data, func.names) {
  
  data %>%
    mutate(`scal. MF` = MF_jing(adf = data, vars = func.names),
           `sum MF` = MF_sum(adf = data, vars = func.names),
           `ave. MF` = MF_av(adf = data, vars = func.names, stand_method = "z_score_abs"),
           `sd. MF` = MF_sd(adf = data, vars = func.names, stand_method = "z_score_abs"),
           `MESLI MF` = MF_mesli(adf = data, vars = func.names),
           `Pasari MF` = MF_pasari(adf = data, vars = func.names),
           `SAM MF` = MF_dooley(adf = data, vars = func.names),
           `ENF MF` = hill_multifunc(adf = data, vars = func.names, scale = 1, HILL = TRUE),
           `Simp. MF` = MF_simpsons_div(adf = data, vars = func.names),
           `Manning.30 MF` = manning_multifunc(adf = data, vars = func.names, thresh = 0.3),
           `Manning.50 MF` = manning_multifunc(adf = data, vars = func.names, thresh = 0.5),
           `Manning.70 MF` = manning_multifunc(adf = data, vars = func.names, thresh = 0.7),
           `thresh.30 MF` = single_threshold_mf(adf = data, vars = func.names, thresh = 0.3),
           `thresh.50 MF` = single_threshold_mf(adf = data, vars = func.names, thresh = 0.5),
           `thresh.70 MF` = single_threshold_mf(adf = data, vars = func.names, thresh = 0.7),
           `Slade.10.90 MF` = MF_slade(adf = data, vars = func.names, A_quant = 0.10, B_quant = 0.90),
           `Slade.40.60 MF` = MF_slade(adf = data, vars = func.names, A_quant = 0.40, B_quant = 0.60),
           `PCA MF` = pca_multifunc(adf = data, vars = func.names)
    )
  
}

### perform a cluster analysis on the jena data
names(jena.dat)

# set up the func.names
func.names <- c("biomass",  "plantCN", "soilC",             
                "soilorgC", "herbi", "micBMC", "phoact",            
                "poll", "rootBM")

jena.mf <- 
  calculate_MF(data = jena.dat, func.names = func.names)
names(jena.mf)

# perform the cluster analysis based on this subset of points
jena.clust <- 
  jena.mf %>% 
  dplyr::select(ends_with(" MF"), -starts_with("sd"))

mds_out <- 
  jena.clust %>% 
  mutate(across(.cols = everything(), ~(.-mean(., na.rm = TRUE))/sd(., na.rm = TRUE) ) ) %>% 
  as.matrix(.) %>% 
  t() %>% 
  dist(x = ., method = "euclidean") %>% 
  metaMDS(comm = .) 

mds_out$points %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column(var = "MF_metrics") %>% 
  ggplot(data = .,
         mapping = aes(x = MDS1, y = MDS2, colour = MF_metrics)) +
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

jena.clust %>%
  cor(., use = "na.or.complete") %>%
  corrplot::corrplot(method = "ellipse")


### simulate data to test the performance of the metrics

# fit distributions to the functions from the Jena data
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

# choose five functions and find the best fitting distributions
func.names

set.seed(548375)
best.dist <- vector("list", length = length(func.names))
for(i in 1:length(func.names)) {
  
  best.dist[[i]] <- Fit.dist(func = jena.dat[[func.names[i] ]])
  
}

# choose distributions that represent a range of possible distributions
best.dist.sub <- best.dist[c(1, 3, 5, 7, 8)]

# set up the number of points to simulate
n.sim <- 1001
n.sim <- as.integer(n.sim)

if ( (class(n.sim) != "integer") | (n.sim %% 2 == 0) ) {
  warning("code is only robust if odd simulations are used because it is based on quantile distributions")
}

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

x6 <- runif(n = n.sim, 0, 100)

x.sim <- 
  lapply(list(x6, x6, x6, x6, x6), function(x) {
  
  # translate the distribution so that it starts at zero
  x <- (x - min(x)) 
  
  s1 <- seq(0, 0.005, 0.001)
  s2 <- seq(0, 1, length = 12)
  s3 <- seq(0.995, 1, 0.001)
  
  q1 <- round(quantile(1:length(x), c(s1, s2, s3)) , 0)
          
  sort(x)[q1]
  
  # seq(z[1], z[2], length.out = 10)
  
})

# create a grid of means and sd's and sample points around that grid
# grid can be based on quantiles of the max-standardised functions

x.sims <- do.call(expand.grid, x.sim)
names(x.sims) <- paste("F_", 1:length(x.sim), sep = "")
x.sims <- as_tibble(x.sims)

# create the simulated dataset to plot
v <- apply(x.sims, 1, sum)
v <- which(v == min(v) | v == max(v))

sim.metric <- x.sims[c(v[1], sample(1:nrow( x.sims), 1000), v[length(v)]), ]
# rm(x.sims)

# calculate the multifunctionality metrics on these simulated data using function
sim.metric <- calculate_MF(data = sim.metric, func.names = func.names)

# add a row.id column
sim.metric$row_id <- 1:nrow(sim.metric)
# View(sim.metric)

# how to standardise the raw function data for plotting?
stand <- "max"
sim.metric <- 
  sim.metric %>%
  mutate(ave_MF = MF_av(adf = sim.metric, vars = func.names, stand_method = stand),
         sd_MF = MF_sd(adf = sim.metric, vars = func.names, stand_method = stand),
         across(.cols = starts_with("F_"), ~standardise_functions(.x, method = stand ) ) )

plot(sim.metric$ave_MF, sim.metric$sd_MF)

# remove duplicates in the df data.frame
df.max <- 
  sim.metric %>%
  filter(across(.cols = starts_with("F_"), ~.x == max(.x, na.rm = TRUE)))

df.min <- 
  sim.metric %>%
  filter(across(.cols = starts_with("F_"), ~.x == min(.x, na.rm = TRUE)))

df.all <- 
  sim.metric

# get four random points
ave.q <- quantile(df.all$ave_MF, c(0.2, 0.45, 0.55,  0.9))
sd.q <- quantile(df.all$sd_MF, c(0.1, 0.45, 0.55, 0.9))

df.plot.1 <- 
  df.all %>%
  filter(ave_MF < ave.q[1] | ave_MF > ave.q[4]) %>%
  filter(sd_MF > sd.q[2], sd_MF < sd.q[3] ) %>%
  mutate(ave_low_high = if_else(ave_MF < ave.q[1], "a", "b" )) %>%
  group_by(ave_low_high) %>%
  sample_n(size = 1, replace = TRUE) %>%
  ungroup() %>%
  dplyr::select(row_id, ave_MF, sd_MF)

df.plot.2 <- 
  df.all %>%
  filter(sd_MF < sd.q[1] | sd_MF > sd.q[4]) %>%
  filter(ave_MF > ave.q[2], ave_MF < ave.q[3]  ) %>%
  mutate(sd_low_high = if_else( sd_MF < sd.q[1], "a", "b" ) ) %>%
  group_by(sd_low_high) %>%
  sample_n(size = 1, replace = TRUE) %>%
  ungroup() %>%
  dplyr::select(row_id, ave_MF, sd_MF)

df.plot <- rbind(df.plot.1, df.plot.2)

row_id_in <- c(df.min$row_id, df.max$row_id, df.plot$row_id) 

df.label <- 
  df.all %>%
  filter(row_id %in% row_id_in) %>%
  dplyr::select(row_id, starts_with("F_"), ave_MF, sd_MF) %>%
  mutate(label = letters[1:length(row_id_in)])


# make a test plot
library(ggrepel)
ggplot( mapping = aes(x = ave_MF, y = sd_MF ) ) +
  geom_point(data = df.all %>% 
               filter(ave_MF != max(ave_MF)) %>%
               filter(ave_MF != 0), 
             size = 3, alpha = 0.2, position = position_jitter(width = 0.005)) +
  geom_hline(yintercept = df.min[["sd_MF"]], linetype = "dashed", colour = "black") +
  geom_vline(xintercept = df.max[["ave_MF"]], linetype = "dashed", colour = "black") +
  geom_point(data = df.min, 
             fill = "white", colour = "black", shape = 24, alpha = 1, size = 5, stroke = 1.25 ) +
  geom_point(data = df.max,
             fill = "white", colour = "black", shape = 21, alpha = 1, size = 5, stroke = 1.25 ) +
  geom_point(data = df.plot,
             fill = "#99CCFF", colour = "black", shape = 21, alpha = 1, size = 5, stroke = 1.25) +
  geom_text(data = df.label, 
            mapping = aes(label = label), size = 4, nudge_y = 0.0025) +
  ylab("SD among functions") +
  xlab("Mean among functions") +
  scale_colour_viridis_c() + 
  theme_meta()

# plot each highlighted point in the test plot
df.label.long <- 
  df.label %>%
  pivot_longer(cols = starts_with("F_"),
               names_to = "function_id", 
               values_to = "function_value")

# iterate this
ggplot(data = df.label.long %>% filter(label == "b"),
       mapping = aes(x = function_id, y = function_value)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.2)) +
  xlab("Function ID") +
  ylab("Function value") +
  theme_meta()


# iterate this
# for each metric, what do we want to plot separately?
# set the metric
names(df.all)
metric <- "SAM MF"

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
             mapping = aes(x = ave_MF, y = sd_MF, shape = metric_negative, colour = get(metric) ),
             size = 3, alpha = 0.5, position = position_jitter(width = 0.005)) +
  scale_colour_viridis_c(option = "D") +
  scale_shape_manual(values=c(16, 15), guide = FALSE)+
  geom_hline(yintercept = df.min[["sd_MF"]], linetype = "dashed", colour = "black") +
  geom_vline(xintercept = df.max[["ave_MF"]], linetype = "dashed", colour = "black") +
  geom_point(data = df.min, 
             mapping = aes(x = ave_MF, y = sd_MF), 
             fill = "white", colour = "black", shape = 24, alpha = 1, size = 3.5, 
             position = position_nudge(x = -0.025), stroke = 1.25 ) +
  geom_point(data = df.max,
             mapping = aes(x = ave_MF, y = sd_MF),
             fill = "white", colour = "black", shape = 21, alpha = 1, size = 3.5,
             position = position_nudge(x = -0.025), stroke = 1.25 ) +
  geom_point(data = df.metric.max,
             mapping = aes(x = ave_MF, y = sd_MF),
             colour = "black", fill = "#FF6600", shape = 21, alpha = 1, size = 3.5, stroke = 1.25) +
  geom_point(data = df.metric.min,
             mapping = aes(x = ave_MF, y = sd_MF),
             colour = "black", fill =  "#FF6600", shape = 24, alpha = 1, size = 3.5, stroke = 1.25) +
  geom_point(data = df.metric.undefined,
             mapping = aes(x = ave_MF, y = sd_MF), 
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

### END
