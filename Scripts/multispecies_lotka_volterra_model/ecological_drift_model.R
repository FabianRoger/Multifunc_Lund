
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code for an ecological drift model (sensu Hubbell 2001) but without dispersal or speciation

# parameter definitions

# lsp: gradient of local species pools (i.e. initial seeded diversity)
# reps: number of replicates for each local species pool size
# rsp: number of species in the regional species pool
# t_steps: number of time steps to model
# n0: initial abundance of all species (i.e. all species together)

# prop_change: expected proportion (Poisson distribution) of individuals that turnover at each time-step

# write an ecological drift model (sensu Hubbell 2001)

drift_model <- function(lsp = c(1, 2, 4, 6),
                        reps = 5,
                        rsp = 12,
                        t_steps = 1000,
                        n0 = 500,
                        prop_change = 0.3) {
  
  # define a truncated poisson distribution
  trunc_pois <- function(n, lambda) {
    x <- rpois(n = n, lambda = lambda)
    y <- ifelse(x == 0, 1, x)
    return(y)
  }
  
  # set up the model
  
  # define the number of patches
  patches <- reps*length(lsp)
  
  # assign a local species spool size to each patch
  lsp.p <- rep(lsp, each = reps)[sample(x = (1:patches), size = patches, replace = FALSE)]
  
  # put this into a data.frame
  df.p <- data.frame(loc.sp.p = lsp.p)
  df.p <- split(df.p, 1:nrow(df.p))
  
  neutral_list <- 
    lapply(df.p, function(patch) {
      
      # get the local species pool size
      lsp.size <- patch$loc.sp.p
      
      # get a set of species from the regional species pool
      lsp.patch <- sample(x = (1:rsp), size = lsp.size, replace = FALSE)
      
      # create a sack of individuals of each species
      patch.sack <- rep(lsp.patch, each = round((n0/lsp.size), 0))
      
      # create an output list of species abundances for each time point
      n_t <- vector("list", length = t_steps)
      
      # fill the first time point with starting abundances
      n_t[[1]] <- patch.sack
      
      # for each time point m
      for(m in seq(from = 2, to = t_steps, by = 1)){
        
        # round the proportion change variable
        n_change <- round(n_t[[m-1]]*prop_change, 0)
        
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
      df.patch$local_species_pool <- lsp.size
      
      return(df.patch)
      
    })
  
  # bind this list into a data.frame
  neutral_data <- bind_rows(neutral_list, .id = "patch")
  
  # convert the time variable into a numeric variable
  neutral_data <- 
    neutral_data %>%
    mutate(time = as.numeric(time))
  
  return(neutral_data)
  
}

# test the ecological drift model
x <- drift_model(lsp = c(1, 2, 4, 6, 9),
                 reps = 5,
                 rsp = 12,
                 t_steps = 500,
                 n0 = 300,
                 prop_change = 0.3)


x %>%
  filter(abundance > 0) %>%
  ggplot(data = .,
         mapping = aes(x = time, y = abundance, colour = as.character(species) )) +
  geom_point() +
  facet_wrap(~patch, scales = "free") +
  theme_classic()

x %>%
  filter(time == last(time)) %>%
  group_by(patch) %>%
  summarise(local_species_pool = mean(local_species_pool),
            realised_richness = sum(if_else(abundance > 0, 1, 0)),
            total_abundance = sum(abundance)) %>%
  ggplot(data = .,
         mapping = aes(x = local_species_pool, y = total_abundance)) +
  geom_point() +
  theme_classic()
