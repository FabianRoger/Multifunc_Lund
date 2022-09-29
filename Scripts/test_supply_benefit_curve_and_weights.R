
# load relevant libraries
library(dplyr)
library(here)

# tell R where to call scripts from
source(here("Scripts/01_general_functions/function_plotting_theme.R"))
source(here("Scripts/01_general_functions/MF_functions_collated.R"))

# simulate different functions
nf <- 5
reps <- 100

mf_out <- vector("list", length = nf)
for (i in 1:nf) {
  
  # draw a function from a normal distribution
  x <- rnorm(n = reps, mean = 10, sd = 2)
  
  # write it into a list
  mf_out[[i]] <- x
  
}

mf_out <- bind_cols(mf_out)
names(mf_out) <- paste0("F_", 1:nf)


# (1) define a supply-benefit curve for each function

# threshold supply-benefit
thresh <- function(x, thresh = 0.5) {
  
  y <- dplyr::if_else(x > (thresh*max(x)), 1, 0)
  return(y)
  
}

# linear supply-benefit
linear <- function(x) { 
  
  return(x) 
  
  }

# stepwise linear supply-benefit
step_linear <- function(x, thresh = 0.10) {
  
  y <- dplyr::if_else(x < (thresh*max(x)), 0, x)
  return(y)
  
}

# exponential supply-benefit
exponential <- function(x) { 
  
  y <- exp(x)
  return(y) 
  
  }

# choose supply-benefit curves for each function
sb_curves <- c("linear", "linear", "thresh", "step_linear", "exponential")
sb_args <- list(NA, NA, 0.5, 0.3, NA)

# apply the different supply-benefit functions
for(i in 1:ncol(mf_out)) {
  
  if(!is.na(sb_args[i])) {
    
    args_in <- c(list(mf_out[[i]]), (sb_args[i]))
    
  } else {
    
    args_in <- list(mf_out[[i]])
    
  }
  
  mf_out[[i]] <-  do.call(sb_curves[i], args_in)
  
}


# (2) transform to a common scale or not?

# transform Y-N
transform <- "Y"

# choose method
method = "max"

# apply the standardisation

if (transform == "Y") {
  
  mf_out <- 
    apply(mf_out, 2, function(x) round(standardise_functions(x = x, method = method), 3) ) %>%
    as_tibble()
  
}

# (3) apply weights to each function

# set-up the weights
ws <- c(1, 0.1, 0.5, 3, 1)
ws <- ( (ws/sum(ws))*nf )

mf <- apply(mf_out, 1, sum)

x <- sweep(mf_out, MARGIN = 2, ws, `*`)
x1 <- apply(x, 1, sum)

plot(mf, x1)








