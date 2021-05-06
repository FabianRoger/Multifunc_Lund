
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code for an ecological drift model (sensu Hubbell 2001) but without dispersal or speciation

# Next steps:

# fix neutral model when abundances drop very low
# sample abundances rather than individuals to fix this...


# parameter definitions

# lsp: gradient of local species pools (i.e. initial seeded diversity)
# reps: number of replicates for each local species pool size
# rsp: number of species in the regional species pool
# t_steps: number of time steps to model
# n0: initial abundance of all species (i.e. all species together)

# prop_change: expected proportion (Poisson distribution) of individuals that turnover at each time-step

# n_repeats: how many times to run the model with current parameters

# write an ecological drift model (sensu Hubbell 2001)
drift_model <- function(lsp = c(1, 2, 4, 6),
                        reps = 5,
                        rsp = 12,
                        t_steps = 1000,
                        n0 = 500,
                        prop_change = 0.05,
                        n_repeats = 5) {
  
  # load the dplyr library
  library(dplyr)
  
  # define a truncated poisson distribution
  trunc_pois <- function(n, lambda) {
    x <- rpois(n = n, lambda = lambda)
    y <- ifelse(x == 0, 1, x)
    return(y)
  }
  
  # set up the fixed model parameters
  
  # define the number of patches
  patches <- reps*length(lsp)
  
  # assign a local species spool size to each patch
  lsp.p <- rep(lsp, each = reps)
  
  # get a list of patches with a sack of species for each patch
  patch.sack <- vector("list", length = length(lsp.p))
  for (i in 1:length(lsp.p)) {
    
    x <- sample(x = (1:rsp), size = lsp.p[i], replace = FALSE)
    
    patch.sack[[i]] <- rep(x, each = round((n0/lsp.p[i]), 0))
    
  }
  
  
  # run the model n-times
  
  model_out <- vector("list", length = n_repeats)
  for (i in 1:n_repeats){
    
    neutral_list <- 
      lapply(patch.sack, function(patch) {
        
        # create an output list of species abundances for each time point
        n_t <- vector("list", length = t_steps)
        
        # fill the first time point with starting abundances
        n_t[[1]] <- patch
        
        # for each time point m
        for(m in seq(from = 2, to = t_steps, by = 1)){
          
          # round the proportion change variable
          n_change <- round(length(n_t[[m-1]])*prop_change, 0)
          n_change <- ifelse(n_change == 0, 1, n_change)
          
          # kill n individuals drawn from a poisson distribution
          post_death <- n_t[[m-1]][-sample(x = 1:length(n_t[[m-1]]), size = trunc_pois(n = 1, n_change))]
          
          # draw new recruits
          z <- trunc_pois(n = 1, n_change)
          new_recruits <- post_death[sample(x = 1:length(post_death), size = z, replace = FALSE)]
          
          # join the post_death and new_recruits and write to second time-step
          n_t[[m]] <- c(post_death, new_recruits)
          
        }
        
        # summarise model results into a data.frame
        df.patch <- lapply(n_t, function(x) {
          
          sp.abun <- vector(length = rsp)
          for (j in 1:rsp) {
            
            sp.abun[j] <- sum(j == x)
            
          }
          
          df <- data.frame(species = 1:rsp,
                           abundance = sp.abun)
          
          return(df)
          
        })
        
        df.patch <- bind_rows(df.patch, .id = "time")
        
        # add local species pool information
        df.patch$local_species_pool <- length(unique(patch))
        
        return(df.patch)
        
      })
    
    # bind this list into a data.frame
    neutral_data <- bind_rows(neutral_list, .id = "patch")
    
    # convert the time variable into a numeric variable
    neutral_data <- 
      neutral_data %>%
      mutate(time = as.numeric(time))
    
    model_out[[i]] <- neutral_data
    
  }
  
  return( bind_rows(model_out, .id = "run") )
  
}

# test the ecological drift model
x <- drift_model(lsp = c(1, 2, 4, 6),
                 reps = 5,
                 rsp = 9,
                 t_steps = 500,
                 n0 = 30,
                 prop_change = 0.05,
                 n_repeats = 1)

library(ggplot2)
head(x)
x %>%
  filter(abundance > 0) %>%
  mutate(species = as.character(species)) %>%
  group_by(local_species_pool, patch, time, species) %>%
  summarise(abundance = median(abundance)) %>%
  ggplot(data = .,
         mapping = aes(x = time, y = abundance, colour = species )) +
  geom_point() +
  scale_colour_viridis_d() +
  facet_wrap(~patch, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom")

x %>%
  filter(abundance > 0) %>%
  mutate(species = as.character(species)) %>%
  group_by(local_species_pool, patch, time, species) %>%
  summarise(abundance = median(abundance)) %>%
  filter(patch == 5) %>%
  filter(time > 450) %>%
  View()


x %>%
  filter(time == last(time)) %>%
  group_by(patch) %>%
  summarise(local_species_pool = mean(local_species_pool),
            realised_richness = sum(if_else(abundance > 0, 1, 0)),
            total_abundance = sum(abundance)) %>%
  ggplot(data = .,
         mapping = aes(x = local_species_pool, y = total_abundance)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(legend.position = "bottom")
