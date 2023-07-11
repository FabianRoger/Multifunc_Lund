#'
#' @title Test whether slope changes with number of functions
#' 
#' @description Here we simulate relationships between biodiversity and
#' ecosystem multifunctionality (EMF) and then we test whether the relationship
#' between biodiversity and EMF changes by considering different numbers
#' of functions.
#'

# simulate purely linear relationships between biodiversity and functioning

# do the tests on these relationships

# simulate a mix of linear and non-linear relationships between biodiversity and functioning

# 

# set the number of functions
n_func <- 10

# set the number of plots
n <- 100

# get the diversity of these different plots
div <- rpois(n = n, lambda = 10)
range(div)
hist(div)

# simulate functions from a power function

fsim_list <- 
  lapply(1:n_func, function(x) {
  
  # draw the beta slope parameter
  b <- runif(n = 1, min = -2, max = 2)
  print(b)
  
  # simulate the function value
  x <- (div^b)
  
  # translate to positive values
  y <- x + abs( min(x) )
  
  # standardise by the maximum value and add random noise
  z <- (y/max(y)) + rnorm(n = length(y), mean = 0, sd = 0.1)
  
  # return the maximum standardised function value
  return(z)
  
})

# pull the functions into a data.frame
fsim <- do.call("cbind", fsim_list)

# bind the simulation data for efficient plotting
fsim_plot <- dplyr::bind_cols( dplyr::tibble(plot = 1:n, div = div), 
                               dplyr::bind_cols(fsim_list))

# rename the columns
names(fsim_plot) <- c("plot", "div", paste0("F", 1:n_func))

# pull into the long format
fsim_plot <- 
  fsim_plot |>
  tidyr::pivot_longer(cols = contains("F"),
                      names_to = "Function", 
                      values_to = "Value")

p1 <- 
  ggplot(data = fsim_plot,
         mapping = aes(x = div, y = Value, group = Function)) +
  geom_smooth(size = 0.5, colour = "grey", se = FALSE) +
  ylab("Function value (0-1)") +
  xlab("Species richness") +
  theme_test() +
  theme(axis.text = element_text(colour = "black"))
plot(p1)





