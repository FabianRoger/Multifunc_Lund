
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code the Stachova and Leps (2010) simulation model

# parameter definitions

# lsp: gradient of local species pools (i.e. initial seeded diversity)
# mono: options for including monocultures:
# - "all" includes all possible monocultures from the regional species pool (i.e. rsp)
# - "random" draws monocultures randomly as per the other local species pool treatments
# - "none" excludes monocultures
# reps: number of replicates for each local species pool size
# technical_reps: how many repeats of each treatment exactly
# rsp: number of species in the regional species pool
# t_steps: number of time steps to model
# n0: initial abundance of all species
# ext.thresh = threshold abundance of species extinction

# a_mean = mean interspecific competition of a truncated normal distribution
# a_sd = sd interspecific competition of a truncated normal distribution
# a_min = min interspecific competition of a truncated normal distribution
# a_max = max interspecific competition of a truncated normal distribution
# a_spp = intraspecific competition for each species
# a_scale = scaling coefficient of competition coefficients

# kmin = minimum carrying capacity (uniform distribution)
# kmax = maximum carrying capacity (unifrom distribution)

# rmin = minimum intrinsic growth rate (uniform distribution)
# rmax = maximum intrinsic growth rate (uniform distribution)

# n_repeats: how many times to run the model with current parameters

s_l_function <- function(lsp = c(5, 10, 15, 20, 25),
                         mono = "all",
                         reps = 10,
                         technical_reps = 1,
                         rsp = 50,
                         t_steps = 10,
                         n0 = 20,
                         a_mean = 0.8, a_sd = 0.2, a_min = 0.2,
                         a_max = 1.2, a_spp = 1, sim.comp = "sym", a_scale = 1,
                         k_min = 20, k_max = 150,
                         r_min = 0.01, r_max = 0.5,
                         n_repeats = 1){
  
  # check that the correct packages are installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires dplyr to be installed"
  )
  
  if(! "tidyr" %in% installed.packages()[,1]) stop(
    "this function requires tidyr to be installed"
  )
  
  if(! "truncnorm" %in% installed.packages()[,1]) stop(
    "this function requires truncnorm to be installed"
  )
  
  # install the pipe from dplyr to be used in further calculations
  `%>%` <- dplyr::`%>%`
  
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
  
  
  # set-up the model parameters for the species
  
  # generate matrix of competition coefficients
  alpha <- matrix( (truncnorm::rtruncnorm(n = rsp*rsp, 
                                          a = a_min, b = a_max, 
                                          mean = a_mean, 
                                          sd = a_sd)), rsp, rsp)
  
  # make the matrix symmetric if chosen
  if (sim.comp == "sym") {
    alpha[lower.tri(alpha)] = t(alpha)[lower.tri(alpha)]
  }
  
  # make all intraspecific competition coefficients equal to 1
  diag(alpha) <- rep(a_spp, rsp)
  
  # scale the competition-coefficient matrix
  alpha <- alpha*a_scale
  
  # generate the carrying capacities (K) for each species from the uniform distribution
  k <- runif(n = rsp, min = k_min, max = k_max)
  
  # generate growth rates values (r) for each species from the uniform distribution
  r <- runif(n = rsp, min = r_min, max = r_max)
  
  # set up the mixtures
  
  # assign a local species pool size to each mixture patch
  lsp.p <- rep(lsp, each = reps)
  
  mix.patch <- vector("list", length = length(lsp.p))
  for (i in 1:length(lsp.p)) {
    
    mix.patch[[i]] <- sample(x = (1:rsp), size = lsp.p[i], replace = FALSE)
    
  }
  
  # set up the monocultures
  
  if (mono == "all") {
    
    mono.patch <- vector("list", length = rsp)
    for (i in 1:rsp) {
      
      mono.patch[[i]] <- i
      
    }
    
  } else if (mono == "random") {
    
    mono.patch <- vector("list", length = reps)
    for (i in 1:reps) {
      
      mono.patch[[i]] <- sample(x = (1:rsp), size = 1, replace = FALSE)
      
    } 
    
  } else if (mono == "none") {
    
    mono.patch <- NULL
    
  } else {
    
    stop("error! choose appropriate monoculture option")
  }
  
  # add the monocultures to the mixtures
  mono_mix <- c(mono.patch, mix.patch)
  
  # set the number of technical replicates i.e. replicates of same exact composition
  if( technical_reps < 1) {
    
    stop("error! technical replicates cannot be less than one")
    
  } else if (technical_reps == 1) {
    
    patch.list <- mono_mix  
    
  } else if (technical_reps > 1) {
    
    patch.list <- c(mono_mix)
    for (i in 1:(technical_reps - 1)) {
      
      patch.list <- c(patch.list, mono_mix)
      
    }
    
  } else {
    
    stop("error! choose appropriate technical replicate")
    
  }
  
  
  # run model n number of times
  
  model_out <- vector("list", length = n_repeats)
  
  for (i in 1:n_repeats){
    
    # apply over all patches
    sim.out <- 
      lapply(patch.list, function(patch) {
        
        # get the local species pool size
        lsp.size <- length(patch)
        
        # code a nested for loop: for each time and for each species
        
        # create a vector of starting values for each species in the regional species pool
        n_vals <- rep(0, times = rsp)
        
        # add starting values for species in the lsp
        n_vals[patch] <- (n0/lsp.size)
        
        # create an output list of species abundances for each time point
        n_t <- vector("list", length = t_steps)
        
        # fill the first time point with starting abundances
        n_t[[1]] <- n_vals
        
        # for each time point m
        for(m in seq(from = 2, to = t_steps, by = 1)){
          
          # for each species g
          for (g in 1:rsp ) {
            
            # get the intrinsic growth rate (r) for each species (g)
            r.spp <- r[g]
            
            # get the carrying capacity (k) for each species (g)
            k.spp <- k[g]
            
            # define the first term
            t1 <- n_t[[m-1]][g] 
            
            # define the second term of the equation
            t2 <- (t1*r.spp)*(1 - ( sum( alpha[, g]*n_t[[m-1]] )/k.spp )   )
            
            # combine these terms
            t3 <- (t1 + t2)
            
            # use three terms in the full equation
            if (t3 > 0) {
              n_t[[m]][g] <- rpois(n = 1, lambda = t3)
            } else {
              n_t[[m]][g] <- 0
            }
            
          }
          
        }
        
        # collapse this into a dataframe
        df_n_t <- as.data.frame(do.call(rbind, n_t))
        names(df_n_t) <- paste("sp_", 1:rsp, sep = "")
        
        # add a column for the time-point
        df_n_t$time <- seq(from = 1, to = t_steps, by = 1)
        
        # pull this into two columns
        df_n_t <- 
          df_n_t %>%
          tidyr::pivot_longer(cols = -time,
                              names_to = "species",
                              values_to = "abundance") %>%
          dplyr::arrange(time, species)
        
        # add local species pool information
        df_n_t$local_species_pool <- lsp.size
        
        # add species composition information
        df_n_t$composition <- paste(unique(patch), collapse = "")
        
        return(df_n_t)
        
      }
      )
    
    library(dplyr)
    
    # bind the rows together
    bef.res <- bind_rows(sim.out, .id = "patch")
    
    # reorder the columns
    bef.res <- 
      bef.res %>%
      dplyr::select(patch, time, local_species_pool, composition, species, abundance)
    
    model_out[[i]] <- bef.res
    
  }
  
  data_out <- bind_rows(model_out, .id = "model_run")
  
  data_out <- 
    data_out %>%
    arrange(model_run, patch, time, local_species_pool, composition, species, abundance)
  
  # write these outputs into a list
  output.list <- list(data.raw = data_out,
                      spp.info = list(competition.coefficients = alpha,
                                      r.vals = r,
                                      k.vals = k))
  
  return(output.list)
    
}
