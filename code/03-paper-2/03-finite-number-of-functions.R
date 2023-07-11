
# simulate a biodiversity gradient

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# set a seed for reproducibility
set.seed(495074209)

# set the number of functions
n_func <- 20

# set the number of plots
n <- 100

# get the diversity of these different plots
div <- rpois(n = n, lambda = 10)
range(div)
hist(div)

# simulate a set of coefficients for the effect of biodiversity on functioning
fslope <- rnorm(n = n_func, mean = 0.1, 0.1)
mean(fslope)
range(fslope)
hist(fslope)

# simulate a set of functions
fsim_list <-
  
  lapply(fslope, function(x) {
  
  y <- (div*x) + rnorm(n = length(div), mean = 0, sd = 0.5)
  z <- y + abs( min(y) )
  return(z/max(z))
  
})

# pull the functions into a data.frame
fsim <- do.call("cbind", fsim_list)

# fsim plot
fsim_plot <- bind_cols(tibble(plot = 1:n, div = div), 
                       bind_cols(fsim_list))
names(fsim_plot) <- c("plot", "div", paste0("F", 1:n_func))

# pull into the long format
fsim_plot <- 
  fsim_plot %>%
  pivot_longer(cols = contains("F"),
               names_to = "Function", 
               values_to = "Value")

p1 <- 
  ggplot(data = fsim_plot,
       mapping = aes(x = div, y = Value, group = Function)) +
  geom_line(stat = "smooth", method = "lm", size = 0.5, alpha = 0.3) +
  ylab("Function value (0-1)") +
  xlab("Species richness") +
  theme_test() +
  theme(axis.text = element_text(colour = "black"))

# plot the relationship between average EMF and diversity
fsim_mu <- tibble(div = div, ave_EMF = apply(fsim, 1, mean))

p2 <- 
  ggplot(data = fsim_mu,
       mapping = aes(x = div, y = ave_EMF)) +
  geom_jitter(width = 0.1, shape = 1, alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", size = 0.5, colour = "black") +
  ylab("True average EMF") +
  xlab("Species richness") +
  theme_test()+
  theme(axis.text = element_text(colour = "black"))

p12 <- 
  cowplot::plot_grid(p1, p2, labels = c("a", "b"),
                     label_fontface = "plain", label_size = 11
                     )
plot(p12)

ggsave(filename = "github/multifunc_sims/figures/fig_x.png", p12,
       units = "cm", width = 15, height = 7.5)

# calculate the true average multifunctionality slope

# get the population slope
mu <- lm(apply(fsim, 1, mean) ~ div)
mu <- summary(mu)
true_B <- mu$coefficients[2,][["Estimate"]]
print(true_B)
true_SE <- mu$coefficients[2,][["Std. Error"]]
print(true_SE)

# get a vector of the number of functions to sample
sample_funcs <- rep(c(2, 6, 10, 14, 18), each = 20)

# get the individual slope of each function
samp_B <- vector(length = length(sample_funcs) )
samp_SE <- vector(length = length(sample_funcs) )
for(i in 1:length(sample_funcs)) {
  
  funcs <- sample(1:ncol(fsim), sample_funcs[i])
  
  x <- lm(apply(fsim[,funcs], 1, mean) ~ div)
  y <- summary(x)
  samp_B[i] <- y$coefficients[2,][["Estimate"]]
  samp_SE[i] <- y$coefficients[2,][["Std. Error"]]
  
}

# pull into a data.frame
samp_df <- tibble(id = 1:length(samp_B),
                  n_funcs = sample_funcs,
                  slope = samp_B,
                  slope_SE = samp_SE)

# add some random noise to the n_funcs variable
samp_df$n_funcs_dodge <- samp_df$n_funcs + rnorm(n = nrow(samp_df), 0, 0.2)

# pull the true effect into a data.frame
true_df <- tibble(n_funcs = 20,
                  slope = true_B,
                  slope_SE = true_SE)

p3 <- 
  ggplot() +
  geom_hline(yintercept = true_df$slope, linetype = "dashed", colour = "red") +
  geom_point(data = samp_df,
             mapping = aes(x = n_funcs_dodge, y = slope),
             shape = 16, size = 1.5, alpha = 0.3) +
  geom_errorbar(data = samp_df,
                mapping = aes(x = n_funcs_dodge, 
                              ymin = slope-slope_SE, ymax = slope+slope_SE),
                width = 0, alpha = 0.25) +
  geom_point(data = true_df,
             mapping = aes(x = n_funcs, y = slope),
             shape = 18, size = 3, alpha = 0.75, colour = "red") +
  geom_errorbar(data = true_df,
                mapping = aes(x = n_funcs, 
                              ymin = slope-slope_SE, ymax = slope+slope_SE),
                width = 0, alpha = 0.75, colour = "red") +
  scale_x_continuous(limits = c(0, 21), breaks = c(unique(sample_funcs), 20) ) +
  xlab("Number of functions measured") +
  ylab("Slope est. (+-SE)") +
  theme_test() +
  theme(axis.text = element_text(colour = "black"))
plot(p3)

ggsave(filename = "github/multifunc_sims/figures/fig_y.png", p3,
       units = "cm", width = 7.5, height = 7.5)








