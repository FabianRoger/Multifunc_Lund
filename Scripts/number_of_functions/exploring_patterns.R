
sim.n.out <- read_csv(here("data/sim_n_functions.csv"))

df <- 
  sim.n.out %>%
  filter(mod_id %in% first(mod_id)) %>%
  select(model_run, function_matrix, mod_id, func.comb.id, n.func.id, number_of_functions, sd_funcs, cv_funcs, range_funcs, sd_cor, cv_cor, range_cor) %>%
  distinct()

df$mod_id

plot(df$number_of_functions, df$sd_cor)

# load the raw drift model data
drift_mod_dat <- read_csv(here("data/drift_model_n_functions_processed.csv"))
head(drift_mod_dat)

df2 <- 
  drift_mod_dat %>%
  filter(mod_id %in% df$mod_id)

df2 %>% View()

df_F <- 
  df2 %>%
  select(starts_with("F_"))

f.v <- c("F_1","F_2", "F_3", "F_4", "F_5")
f.v2 <- f.v[c(3, 5)]

x <- 
  df_F %>%
  select( all_of(f.v2) ) %>%
  mutate(across(.cols = all_of(f.v2), standardise)) %>%
  apply(X = ., MARGIN = 1, sd)

plot(df2$local_species_pool, x)
cor(df2$local_species_pool, x, method = "spearman")


n_i = rep(5, 5)
sd_i <- rep(2, 5)
sd_i <- 1:5
xi <- 1:5

xy <- vector("list", length = n )
for (i in 1:length(n_i)) {
  
  xy[[i]] <- data.frame(x = rep(xi[i], 10000), 
                        y =rnorm(n = 10000, mean = n_i[i], sd = sd_i[i] ) )
  
}

dfx <- bind_rows(xy)

lm.y <- lm(y ~ x, data = dfx)

summary(lm.y)


# load the raw drift model data
drift_mod_dat <- read_csv(here("data/drift_model_n_functions_processed.csv"))
head(drift_mod_dat)

df2 <- 
  drift_mod_dat %>%
  filter(mod_id %in% sample(drift_mod_dat$mod_id, 1) )

names(df2)

fcomb <- get.function.combinations(function.names = c("F_1", "F_2", "F_3", "F_4", "F_5"))

df2 <- 
  df2 %>%
  select(local_species_pool, starts_with("F_"))

df2 <- 
  df2 %>%
  mutate(across(.cols = starts_with("F_"), standardise))

f_mat <- 
  df2 %>%
  select(starts_with("F_"))

sd_out <- vector("list", length = length(fcomb))
for(i in 1:length(fcomb)){
  
  z <- 
    f_mat %>%
    select(all_of(fcomb[[i]]))
  
  lm.n <- lm(apply(z, 1, sd) ~ df2$local_species_pool)
  lm.m <- lm(apply(z, 1, mean) ~ df2$local_species_pool)
  lm.p <- lm( (apply(z, 1, mean) - apply(z, 1, sd))  ~ df2$local_species_pool)
  
  range_in <- apply(z, 1, sd) %>% range(.)
  
  sd_out[[i]] <- tibble(func_id = paste(fcomb[[i]], collapse = "."),
                        n_funcs = length(fcomb[[i]]),
                        sd_richness_slope = lm.n$coefficients[[2]],
                        mean_richness_slope = lm.m$coefficients[[2]],
                        pas_richness_slope = lm.p$coefficients[[2]],
                        min_sd = range_in[1],
                        max_sd = range_in[2],
                        range_sd = diff(range_in),
                        mean_sd = mean(apply(z, 1, sd)),
                        mean_mean = mean(apply(z, 1, mean))
                        
                        )
}  

bind_rows(sd_out) %>%
  ggplot(data = .,
         mapping = aes(n_funcs, sd_richness_slope)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm")

bind_rows(sd_out) %>%
  ggplot(data = .,
         mapping = aes(n_funcs, mean_richness_slope)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm")

bind_rows(sd_out) %>%
  ggplot(data = .,
         mapping = aes(n_funcs, pas_richness_slope)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm")

bind_rows(sd_out) %>%
  ggplot(data = .,
         mapping = aes(n_funcs, range_sd)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm")

bind_rows(sd_out) %>%
  ggplot(data = .,
         mapping = aes(n_funcs, mean_sd)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm")

bind_rows(sd_out) %>%
  ggplot(data = .,
         mapping = aes(n_funcs, mean_mean)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm")

# we need to show why that decrease in the relationship between sd ecosystem functions
# and species richness with number of functions matters for the multifunctionality calculations

# we need to show why this is the case... Why is that when we have many functions,
# the standard deviation decreases with species richness

# the rate at which standard deviation decreases with species richness is higher 
# when more functions are considered... but why is this the case...


reps <- 100
n1 <- 5

m <- matrix(nrow = reps, ncol = n1, data = rnorm(n = reps*n1))

n2 <- get.function.combinations(function.names = c(1:n1))
# n2 <- n2[lapply(n2, length) > 8]

sd_sim <- vector("list", length = length(n2))
for (i in 1:length(n2) ) {
  
  ar.x <- apply(m[, n2[[i]]], MARGIN = 1, sd)
  
  sd_sim[[i]] <- tibble(numbers = length(n2[[i]]),
                        sd_mean = mean(ar.x),
                        sd_sd = sd(ar.x),
                        sd_range = diff(range(ar.x)) )
                        
  
}
  
sd_sim %>% 
  bind_rows() %>%
  ggplot(data = .,
         mapping = aes(x = numbers, y = sd_range)) +
  geom_point()

sd_sim %>% 
  bind_rows() %>%
  ggplot(data = .,
         mapping = aes(x = numbers, y = sd_mean)) +
  geom_point() +
  geom_smooth()

# on average, standard deviation increases when we increase the number of functions
# why would that lead to a decreasing slope between richness and sd?

sd_sim %>% bind_rows() %>%
  ggplot(data = .,
         mapping = aes(x = numbers, y = sd_sd)) +
  geom_point()

# overall variation in the standard deviation among a set of samples
# decreases when we consider more functions

# does this explain why SD decreases faster with species richness when many functions
# are considered?

yn = 10000

x <- runif(n = yn)

y <- -1.5*x + rnorm(n = yn, sd = 1)

lm(y~x)

y1 <- -1.5*x + rnorm(n = yn, sd = 10)

lm(y1~x)

x <- rnorm(n = 2)
x
y <- get.function.combinations(function.names = c(1:n1))
y

m.out <- vector(length = length(y))
sd.out <- vector(length = length(y))
for(i in 1:length(y)) {
  
  l <- x[y[[i]]]
  
  m.out[i] <- mean(l, na.rm = TRUE)
  sd.out[i] <- sd(l, na.rm = TRUE)
  
  
}

tibble(n_funcs = unlist(lapply(y, length)),
       mean.out = m.out,
       sd.out = sd.out) %>%
  ggplot(data = .,
         mapping = aes(x = n_funcs, y = mean.out)) +
  geom_point() +
  geom_smooth(method = "lm")

tibble(n_funcs = unlist(lapply(y, length)),
       mean.out = m.out,
       sd.out = sd.out) %>%
  ggplot(data = .,
         mapping = aes(x = n_funcs, y = sd.out)) +
  geom_point() +
  geom_smooth(method = "lm")
  
# is it possible that standard deviation increases when considering more things?
# this shouldn't be the case... but for small enough samples... it could be true...



y <- lapply(x, function(y){ data.frame(v1 = y,
                                       v2 = rnorm(y) )})

z <- bind_rows(y)

head(z)

z.sum <- z %>%
  group_by(v1) %>%
  summarise(v2 = sd(v2))

ggplot(data = z.sum,
       mapping = aes(x = v1, y = v2)) +
  geom_point() +
  geom_smooth(method = "lm")







