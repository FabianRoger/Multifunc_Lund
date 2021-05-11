
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code for an ecological drift model (sensu Hubbell 2001) but without dispersal or speciation

# parameter definitions

# lsp: gradient of local species pools (i.e. initial seeded diversity) but not monocultures (must be positive integers)
# mono: options for including monocultures:
# - "all" includes all possible monocultures from the regional species pool (i.e. rsp)
# - "random" draws monocultures randomly as per the other local species pool treatments
# - "none" excludes monocultures
# reps: number of replicates for each local species pool size
# technical_reps: how many repeats of each treatment exactly
# rsp: number of species in the regional species pool
# t_steps: number of time steps to model
# n0: initial abundance of all species (i.e. all species together)

# prop_change: expected proportion (Poisson distribution) of individuals that turnover at each time-step

# n_repeats: how many times to run the model with current parameters

# write an ecological drift model (sensu Hubbell 2001)
drift_model <- function(lsp = c(2, 4, 6),
                        mono = "all",
                        reps = 3,
                        technical_reps = 1,
                        rsp = 9,
                        t_steps = 1000,
                        n0 = 500,
                        prop_change = 0.05,
                        n_repeats = 1) {
  
  # check that the correct packages are installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires dplyr to be installed"
  )
  
  # load the dplyr library
  library(dplyr)
  
  # define a truncated poisson distribution
  trunc_pois <- function(n, lambda) {
    x <- rpois(n = n, lambda = lambda)
    y <- ifelse(x == 0, 1, x)
    return(y)
  }
  
  # perform error checking on arguments
  
  # lsp values cannot be less than or equal to 1
  if ( any(lsp <= 1) ) {
    stop("error! lsp cannot be one because monocultures are specified via mono argument")
  }
  
  # lsp values must be integers
  if ( any(lsp%%1 != 0 ) ) {
    stop("error! lsp must be integers")
  }
  
  # number of replicates (reps) must be at least 1 otherwise a default of 1 is assigned
  if (reps < 1) {
    stop("reps set to less than 1, choose appropriate number of replicates")
  }
  
  # make sure proportion change is between 0 and 1
  if ( (prop_change <= 0) | (prop_change >= 1) ) {
    stop("proportion change argument must be between 0 and 1")
  }
  
  # make sure that n0 is always greater than two times max lsp
  if (n0 < 2*max(lsp)) {
    stop("number of starting individuals (n0) must be at least twice as large as the maximum lsp")
  }
  
  # make sure the rsp is greater than the maximum lsp
  if (rsp < max(lsp)) {
    stop("regional species pool (rsp) must be at greater than or equal to maximum diversity lsp")
  }
  
  # make sure the number of time-steps is greater than or equal to 1
  if (t_steps < 2 ) {
    stop("time-steps less than 2, choose more time steps")
  }
  
  # make sure there is at least one model run
  if (n_repeats < 1) {
    stop("fewer than one model repeat select")
  }
  
  
  # set up the mixtures
  
  # assign a local species pool size to each mixture patch
  lsp.p <- rep(lsp, each = reps)
  
  # get a list of patches with a sack of species for each patch
  mix.sack <- vector("list", length = length(lsp.p))
  for (i in 1:length(lsp.p)) {
    
    x <- sample(x = (1:rsp), size = lsp.p[i], replace = FALSE)
    
    mix.sack[[i]] <- rep(x, each = round((n0/lsp.p[i]), 0))
    
  }
  
  # set up the monocultures
  if (mono == "all") {
    
    mono.sack <- vector("list", length = rsp)
    for (i in 1:rsp) {
      
      mono.sack[[i]] <- rep(i, each = round(n0, 0))
      
    }
    
  } else if (mono == "random") {
    
    mono.sack <- vector("list", length = reps)
    for (i in 1:reps) {
      
      x <- sample(x = (1:rsp), size = 1, replace = FALSE)
      
      mono.sack[[i]] <- rep(x, each = round(n0, 0))
      
    } 
    
  } else if (mono == "none") {
    
    mono.sack <- NULL
    
  } else {
    
    stop("error! choose appropriate monoculture option")
  }
  
  # add the monocultures to the mixtures
  mono_mix <- c(mono.sack, mix.sack)
  
  # set the number of technical replicates i.e. replicates of same exact composition
  if( technical_reps < 1) {
    
    stop("error! technical replicates cannot be less than one")
    
  } else if (technical_reps == 1) {
    
    patch.sack <- mono_mix  
    
  } else if (technical_reps > 1) {
    
    patch.sack <- c(mono_mix)
    for (i in 1:(technical_reps - 1)) {
      
      patch.sack <- c(patch.sack, mono_mix)
      
    }
    
  } else {
    
    stop("error! choose appropriate technical replicate")
    
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
        
        for(m in seq(from = 2, to = t_steps, by = 1)){
          
          # get total number of individuals in patch
          l <- length(n_t[[m-1]])
          
          # round the proportion change variable
          y <- round(l*prop_change, 0)
          n_change <- ifelse(y == 0, 1, y)
          
          # get the number of individuals to kill
          n_kill <- trunc_pois(n = 1, n_change)
          
          # if total length is less than or equal to or if number killed exceeds total population
          # assign an empty vector
          if ( (l <= 1) | (n_kill >= l) ) {
            
            n_t[[m]] <- vector()
            
          } else {
            
            # kill n individuals drawn from a poisson distribution
            post_death <- n_t[[m-1]][-sample(x = 1:l, size = n_kill)]
            
            # draw new recruits
            z <- trunc_pois(n = 1, n_change)
            
            if (z > length(post_death)) {
              
              n_t[[m]] <- c(post_death)
              
            } else {
              
              # draw new recruits
              new_recruits <- post_death[sample(x = 1:length(post_death), size = z, replace = FALSE)]
              
              # join the post_death and new_recruits and write to second time-step
              n_t[[m]] <- c(post_death, new_recruits)
              
            }
            
          }
          
        }
        
        # summarise model results into a data.frame
        df.patch <- lapply(n_t, function(x) {
          
          sp.abun <- vector(length = rsp)
          for (j in 1:rsp) {
            
            sp.abun[j] <- sum(j == x)
            
          }
          
          df <- data.frame(species = paste("sp_", 1:rsp, sep = ""),
                           abundance = sp.abun)
          
          return(df)
          
        })
        
        df.patch <- bind_rows(df.patch, .id = "time")
        
        # add local species pool information
        df.patch$local_species_pool <- length(unique(patch))
        
        # add species composition information
        df.patch$composition <- paste(unique(patch), collapse = "")
        
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
  
  data_out <- bind_rows(model_out, .id = "model_run")
  
  # reorder the columns
  data_out <- 
    data_out %>%
    dplyr::select(model_run, patch, time, local_species_pool, composition, species, abundance) %>%
    arrange(model_run, patch, time, local_species_pool, composition, species, abundance) %>%
    as_tibble()
  
  return( data_out )
  
}

# test the ecological drift model
x <- drift_model(lsp = c(2, 4, 6, 9),
                 mono = "all",
                 reps = 5,
                 technical_reps = 1,
                 rsp = 12,
                 t_steps = 1000,
                 n0 = 500,
                 prop_change = 0.05,
                 n_repeats = 5)

### END
