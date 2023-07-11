
# simulate the change in BEF slope with number of functions

# calculate the true threshold multifunctionality slope

# get the population slope
mu <- lm(apply(fsim, 1, function(x) sum(x > 0.5))/n_func ~ div)
mu <- coef(mu)[2]
print(mu)

# get the individual slope of each function
samp <- vector(length = ncol(fsim))
for(i in 1:ncol(fsim)) {
  
  x <- lm( ( fsim[,i] > 0.5) ~ div)
  samp[i] <- coef(x)[2]
  
}
mean(samp)
