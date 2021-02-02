
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code the Stachova and Leps (2010) simulation model

# parameter definitions

# lsp: gradient of local species pools (i.e. initial seeded diversity)
# reps: number of replicates for each local species pool size
# rsp: number of species in the regional species pool
# t_steps: number of time steps to model
# n0: initial abundance of each species
# ext.thresh = threshold abundance of species extinction

# a_mean = mean interspecific competition of a truncated normal distribution
# a_sd = sd interspecific competition of a truncated normal distribution
# a_min = min interspecific competition of a truncated normal distribution
# a_max = max interspecific competition of a truncated normal distribution
# a_spp = intraspecific competition for each species
# sim.comp = if "sym" then competition is symmetric between species

# kmin = minimum carrying capacity (uniform distribution)
# kmax = maximum carrying capacity (unifrom distribution)

# rmin = minimum intrinsic growth rate (uniform distribution)
# rmax = maximum intrinsic growth rate (uniform distribution)

s_l_function <- function(lsp = c(5, 10, 15, 20, 25),
                         reps = 10,
                         rsp = 50,
                         t_steps = 10,
                         n0 = 3,
                         ext.thresh = 0.2,
                         a_mean = 0.8, a_sd = 0.2, a_min = 0.2,
                         a_max = 1.2, a_spp = 1, sim.comp = "sym",
                         k_min = 3, k_max = 150,
                         r_min = 0.01, r_max = 0.5){
  
  # check that the correct packages are installed
  if(! "dplyr" %in% installed.packages()[,1]) stop(
    "this function requires dplyr to be installed"
  )
  
  if(! "tidyr" %in% installed.packages()[,1]) stop(
    "this function requires tidyr to be installed"
  )
  
  if(! "gtools" %in% installed.packages()[,1]) stop(
    "this function requires gtools to be installed"
  )
  
  if(! "truncnorm" %in% installed.packages()[,1]) stop(
    "this function requires truncnorm to be installed"
  )
  
  # install the pipe from dplyr to be used in further calculations
  `%>%` <- dplyr::`%>%`
  
  # set-up the model
  
  # define the number of patches
  patches <- reps*length(lsp)
  
  # generate matrix of competition coefficients
  alpha <- matrix(truncnorm::rtruncnorm(n = rsp*rsp, 
                                        a = a_min, b = a_max, 
                                        mean = rnorm(n = 1, mean = a_mean, sd = 0.01), 
                                        sd = a_sd), rsp, rsp)
  
  # make the matrix symmetric if chosen
  if (sim.comp == "sym") {
    alpha[lower.tri(alpha)] = t(alpha)[lower.tri(alpha)]
  }
  
  # make all intraspecific competition coefficients equal to 1
  diag(alpha) <- rep(a_spp, rsp)
  
  # generate the carrying capacities (K) for each species from the uniform distribution
  k <- runif(n = rsp, min = k_min, max = k_max)
  
  # generate growth rates values (r) for each species from the uniform distribution
  r <- runif(n = rsp, min = r_min, max = r_max)
  
  # assign a local species spool size to each patch
  lsp.p <- rep(lsp, each = reps)[sample(x = (1:patches), size = patches, replace = FALSE)]
  
  # put this into a data.frame
  df.p <- data.frame(loc.sp.p = lsp.p)
  df.p <- split(df.p, 1:nrow(df.p))
  
  # apply over all patches
  sim.out <- 
    lapply(df.p, function(patch) {
      
      # get the local species pool size
      lsp.size <- patch$loc.sp.p
      
      # get a set of species from the regional species pool
      lsp.patch <- sample(x = (1:rsp), size = lsp.size, replace = FALSE)
      
      # code a nested for loop: for each time and for each species
      
      # create a vector of starting values for each species in the regional species pool
      n_vals <- rep(0, times = rsp)
      
      # add starting values for species in the lsp
      n_vals[lsp.patch] <- n0
      
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
          
          # if a species abundance drops below 0.2 it is considered extinct
          if (n_t[[m]][g] < ext.thresh) { n_t[[m]][g] <- 0 }
          
        }
        
      }
      
      # collapse this into a dataframe
      df_n_t <- as.data.frame(do.call(rbind, n_t))
      names(df_n_t) <- paste0("sp_", 1:rsp)
      
      # add a column for the time-point
      df_n_t$time <- seq(from = 1, to = t_steps, by = 1)
      
      # pull this into two columns
      df_n_t <- 
        df_n_t %>%
        tidyr::pivot_longer(cols = starts_with(match = "sp_"),
                            names_to = "species",
                            values_to = "abundance") %>%
        dplyr::arrange(time, species)
      
      return(df_n_t)
      
    }
    )
  
  library(dplyr)
  
  # bind the rows together
  bef.res <- bind_rows(sim.out, .id = "patch")
  
  # add a column for the local species pool
  lsp.patch <- 
    bef.res %>%
    filter(time == min(time)) %>%
    group_by(patch) %>%
    summarise(local.sp.pool = sum(if_else(abundance > 0, 1, 0)), .groups = "drop")
  
  bef.func <- 
    bef.res %>%
    filter(time == max(time)) %>%
    group_by(patch, time) %>%
    summarise(richness = sum(if_else(abundance > 0, 1, 0)),
              total_abundance = sum(abundance), .groups = "drop")
  
  # join these data.frames together
  sim.proc <- full_join(bef.func, lsp.patch, by = "patch")
  
  # individual species abundances
  sim.spp <- 
    bef.res %>%
    filter(time == max(time))
  
  # write these outputs into a list
  output.list <- list(data.summary = sim.proc,
                      data.raw = sim.spp)
  
  return(output.list)
    
}

# test the function
# df <- s_l_function()

# df[[1]]
# df[[2]]

# plot(df$richness, df$functioning)
# plot(df$local.sp.pool, df$functioning)

# x <- df[[2]]
# hist(x[x$patch == "1", ]$abundance)




