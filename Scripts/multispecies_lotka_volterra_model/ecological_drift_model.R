
# lotka-volterra model
lsp = c(1, 2, 4, 6)
reps = 5
rsp = 12
t_steps = 1000
n0 = 100
a_min = 0
a_max = 1
sim.comp = "sym"
a_scale = 1

# set the average level of interspecific competition
a_mean <- 1

# set the standard deviation of intraspecific competition
a_sd <- 0.1

# set the intraspecific competition value
a_spp <- 1

# set the min and max k-values
k_min = 150
k_max = 500

# set the min and max r-values
r_min = 0.1
r_max = 0.5

df <- s_l_function(lsp = lsp,
             reps = reps,
             rsp = rsp,
             t_steps = t_steps,
             n0 = n0,
             a_mean = a_mean, a_sd = a_sd, 
             a_min = a_min, a_max = a_max,
             a_spp = a_spp, sim.comp = sim.comp,
             k_min = k_min, k_max = k_max,
             r_min = r_min, r_max = r_max)
df$data.summary

library(ggplot2)
ggplot(data = df$data.summary,
       mapping = aes(x = local.sp.pool, y = richness)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()

ggplot(data = df$data.summary,
       mapping = aes(x = local.sp.pool, y = total_abundance)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


# write an ecological drift model (sensu Hubbell 2001)

# define a truncated poisson distribution
trunc_pois <- function(n, lambda) {
  x <- rpois(n = n, lambda = lambda)
  y <- ifelse(x == 0, 1, x)
  return(y)
}

# set-up model parameters
lsp = c(1, 2, 4, 6)
reps = 5
rsp = 12
t_steps = 10
n0 = 20

prop_change = 0.05

prob_local <- 0.9

# define the number of patches
patches <- reps*length(lsp)
patches

# assign a local species spool size to each patch
lsp.p <- rep(lsp, each = reps)[sample(x = (1:patches), size = patches, replace = FALSE)]
lsp.p

# put this into a data.frame
df.p <- data.frame(loc.sp.p = lsp.p)
df.p <- split(df.p, 1:nrow(df.p))
df.p

# for this patch
patch <- df.p[[4]]
patch

# get the local species pool size
lsp.size <- patch$loc.sp.p
lsp.size

# get a set of species from the regional species pool
lsp.patch <- sample(x = (1:rsp), size = lsp.size, replace = FALSE)
lsp.patch

# create a sack of individuals of each species
patch.sack <- rep(lsp.patch, each = round((n0/lsp.size), 0))
patch.sack

# create an output list of species abundances for each time point
n_t <- vector("list", length = t_steps)

# fill the first time point with starting abundances
n_t[[1]] <- patch.sack
n_t[[1]]

# round the proportion change variable
n_change <- round(n0*prop_change, 0)

# for each time point m
for(m in seq(from = 2, to = t_steps, by = 1)){
  
  # kill n individuals drawn from a poisson distribution
  post_death <- n_t[[m-1]][-sample(x = 1:length(n_t[[m-1]]), size = trunc_pois(n = 1, n_change))]
  
  # draw new recruits
  
  # determine to draw from local community or regional pool
  source_pool <- sample(x = c("LC", "RP"), size = 1, prob = c(prob_local, (1-prob_local)))
  
  # draw new recruits either from the local community or the regional pool
  if(source_pool == "LC") {
    
    z <- trunc_pois(n = 1, n_change)
    new_recruits <- post_death[sample(x = 1:length(post_death), size = z, replace = FALSE)]
    
  } else {
    
    z <- trunc_pois(n = 1, n_change)
    new_recruits <- sample(x = lsp.patch, size = z, replace = TRUE, prob = rep(1/lsp.size, lsp.size))
    
  }
  
  # join the post_death and new_recruits and write to second time-step
  n_t[[m]] <- c(post_death, new_recruits)
  
}












# code a nested for loop: for each time and for each species

# create a vector of starting values for each species in the regional species pool
n_vals <- rep(0, times = rsp)
n_vals

# add starting values for species in the lsp
n_vals[lsp.patch] <- (n0/lsp.size)
n_vals

# create an output list of species abundances for each time point
n_t <- vector("list", length = t_steps)

# fill the first time point with starting abundances
n_t[[1]] <- n_vals
n_t[[1]]

# for each time point m
for(m in seq(from = 2, to = t_steps, by = 1)){
  
  # for each species g
  for (g in 1:rsp ) {
    
    
    
  }
  
}










