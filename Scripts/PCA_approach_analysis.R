
# Project: Multifunctionality workshop

# Title: PCA approach analysis

library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(corrplot)
library(grid)
library(gridExtra)

source(here("Scripts/Multifunctionality-Simulations/Multifunc_simulations_functions.R"))
source(here("Scripts/MF_functions_collated.R"))
source(here("Scripts/function_plotting_theme.R"))


# set number of runs to do
n <- 200

# set up a blank list to fill
pca_func <- vector("list", length = n)

# run the loop to generate 10 datasets with different seeds
for (i in (1:n) ) {
  
  set.seed(1555 + i)
  
  # number of species
  specnum <- 10
  
  # number of functions
  funcnum <- 9
  
  # distribution from which to draw function values
  # argument is irrelevant here but needs to be set
  distribution = "runif"
  
  # choose pairwise correlation strength
  COR <- 0
  
  # make correlation matrix (strictly speaking a covariance matrix but for these simulations it does not matter)
  Sigma <- matrix(COR, ncol = funcnum, nrow = funcnum)
  
  # make three 'cluster' of correlated functions
  # Sigma[1:2,1:2] <- 0.7
  # Sigma[4:6,4:6] <- 0.5
  # Sigma[7:9,7:9] <- 0.7
  
  diag(Sigma) <- 1
  
  Sigma
  
  # create function matrix
  FuncMat <- FunctionValue(specnum, funcnum, distribution, min = 0, max = 1)
  FuncMat
  
  # replace function values with correlated function values
  FuncMat_long <- 
    FuncMat %>% 
    pivot_wider(names_from = Functions, 
                values_from = Funcval)
  
  # draw correlated functions (with mean 0)
  corF <- mvrnorm(n = specnum, mu = rep(0, funcnum), Sigma = Sigma)
  
  # shift to positive
  corF <- apply(corF, 2, function(x){ x + abs(min(x)) })
  
  # backtransform to long dataset
  FuncMat_long[, 2:(funcnum+1) ] <- corF
  
  FuncMat <-
    FuncMat_long %>% 
    pivot_longer(cols = starts_with("F"),
                 names_to = "Functions",
                 values_to = "Funcval") 
  
  # extract function names
  func.names <- as.character( unique( FuncMat$Functions))
  
  
  # use these function values to generate data from a biodiversity experiment
  maxrep <- choose(specnum, floor(specnum/2))
  
  # simulate plot x species matrix
  SpecMat <- SpeciesMatrix(specnum = specnum, maxrep = maxrep)
  
  # select method and additional parameters if appropriate by setting the `if` statement to `TRUE`
  if (TRUE) {method = "av"}
  
  if (FALSE) {method = "comp"
  CF = 3
  compfunc = c("F 1", "F 6")
  r = 0.25}
  
  # Average function
  AvFunc <- AverageFunction(SpecMat, FuncMat,
                            method = method, 
                            CF = CF, 
                            compfunc = compfunc,
                            r = r)
  
  # add small normal error
  AvFunc <- 
    AvFunc %>% 
    mutate_at(func.names, function(x) {x + runif(n(), 0, 0.05)})
  
  # standardize by z-score
  AvFunc[,func.names] <- 
    
    apply(AvFunc[,func.names], 2, function(x) {
      
      y <- (x - mean(x))/sd(x)
      
      return(y+abs(min(y)))
    }
    
    )
  
  # add this output into an object
  pca_func[[i]] <- AvFunc
  
}


# calculate the MF metrics on these datasets with different levels of correlation

pca_mf <- 
  
  lapply(pca_func, function(x) {
    
    mutate(x, 
           Meyer_mf = pca_multifunc(adf = x, vars = func.names, standardise = FALSE),
           Av_mf = MF_av(adf = x, vars = func.names),
           Pasari_mf = MF_pasari(adf = x, vars = func.names) )
    
  }
  
  )

pca_mf <- bind_rows(pca_mf, .id = "run")

# rename the functions for plotting
pca_mf_plot <- 
  pca_mf %>%
  rename(`average MF` = Av_mf,
         `PCA MF` = Meyer_mf,
         `Pasari MF` = Pasari_mf)


# plot the relationship between species richness and each function for one random run
row.id <- sample(x = length(unique(pca_mf_plot$run)), size = 1)
row.id

S1a <- 
  pca_mf_plot %>%
  filter(run == row.id) %>%
  pivot_longer(cols = starts_with("F "),
               names_to = "function_id",
               values_to = "function_val") %>%
  ggplot(data = .,
         mapping = aes(x = Richness, y = function_val, colour = function_id)) +
  geom_jitter(alpha = 0.2, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, size = 0.75) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  ylab("function value") +
  xlab("species richness") +
  theme_meta() +
  theme(legend.position = "none")

S1b <- 
  pca_mf_plot %>%
  filter(run == row.id) %>%
  select(starts_with("F ")) %>%
  GGally::ggcorr(data = ., legend.position = "bottom")

S1 <- 
  ggarrange(plotlist = list(S1a, S1b), 
            labels = letters[1:length(list(S1a, S1b))],
            font.label = list(size = 12, color = "black", face = "plain", family = NULL))

ggsave(filename = here("Figures/fig_S1.png"), plot = S1,
       width = 19, height = 10, units = "cm", dpi = 300)

# plot a histogram of slopes of the relationship between average...
hist_slopes <- 
  lapply(split(pca_mf, pca_mf$run), function(x) {
  
  z <- 
    x %>%
    mutate(across(.cols = c("Av_mf", "Pasari_mf", "Meyer_mf"),
                  ~as.numeric( scale(., center = TRUE, scale = TRUE) )) )
  
  lm_1 <- lm(Pasari_mf ~ Av_mf, data = z)
  
  lm_2 <- lm(Meyer_mf ~ Av_mf, data = z)
  
  c("pasari_mf" = lm_1$coefficients["Av_mf"], 
    "pca_mf" = lm_2$coefficients["Av_mf"])
  
}
)

hist_slopes.df <- 
  bind_rows(hist_slopes) %>%
  rename(`Pasari MF` = pasari_mf.Av_mf,
         `PCA MF` = pca_mf.Av_mf) %>%
  pivot_longer(cols = everything(),
               names_to = "metric",
               values_to = "est.")
  
# plot the histograms
S3 <- 
  ggplot(data = hist_slopes.df,
       mapping = aes(x = est.)) +
  geom_histogram(colour = "white", alpha = 0.8, fill = "grey") +
  facet_wrap(~metric) +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed", size = 1) +
  theme_meta()
S3

ggsave(filename = here("Figures/fig_S3.png"), plot = S3,
       width = 19, height = 10, units = "cm", dpi = 300)


# plot average multifunctionality versus pca multifunctionality
S2 <- 
  pca_mf_plot %>%
  filter(run %in% sample((1:n), size = 50 ) ) %>%
  pivot_longer(cols = c(`PCA MF`, `Pasari MF`),
               names_to = "metric",
               values_to = "MF") %>%
  ggplot(data = .,
       mapping = aes(x = `average MF`, y = MF, colour = run)) +
  geom_point(alpha = 0.2, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, size = 0.75) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  facet_wrap(~metric, scales = "free") +
  theme_meta() +
  theme(legend.position = "none")
S2

ggsave(filename = here("Figures/fig_S2.png"), plot = S2,
       width = 19, height = 10, units = "cm", dpi = 300)


### potential additional figures


# plot richness-function plots for average and pca multifunctionality
p3 <- 
  pca_mf_plot %>%
  filter(run %in% sample((1:n), size = 50 ) ) %>%
  pivot_longer(cols = c(`average MF`, `PCA MF`, `Pasari MF`),
               names_to = "metric",
               values_to = "MF") %>%
  ggplot(data = .,
       mapping = aes(x = Richness, y = MF, colour = run)) +
  geom_point(alpha = 0.2, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, size = 0.75) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  facet_wrap(~metric, scales = "free") +
  theme_meta() +
  theme(legend.position = "none")

ggsave(filename = here("Figures/pca_fig_3.png"), plot = p3,
       width = 19, height = 10, units = "cm", dpi = 300)


p4 <- 
  ggplot(data = pca_mf_plot,
         mapping = aes(x = Richness, y = `PCA MF`, colour = `average MF`)) +
  geom_jitter(width = 0.5, alpha = 0.5, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, size = 0.75) +
  scale_colour_viridis_c(option = "C", end = 0.9) +
  facet_wrap(~run) +
  theme_meta()

ggsave(filename = here("Figures/pca_fig_4.png"), plot = p4,
       width = 19, height = 19, units = "cm", dpi = 300)

