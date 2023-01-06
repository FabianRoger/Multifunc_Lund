
# empirical example Fig. 1

# set-up the data.frame
df <- data.frame(F1 = c(0, 3),
                 F2 = c(10, 5))
print(df)

# step 1 - transform to common scale
df$F1 <- ifelse(df$F1<quantile(df$F1, 0.30), 0, 1)
df$F2 <- df$F2/max(df$F2)
print(df)

# step 2 - apply function weights
df <- sweep(df, MARGIN = 2, (c(0.7, 0.3)*2) , `*`)

# calculate Hill number when q = 0
df_p <- apply(df, 1, function(x) round(x/sum(x), 3) )
df_p <- as.data.frame(t(df_p))
print(df_p)

cor(df_p)

vegan::renyi(df_p, scales = 0, hill = TRUE)

vegan::renyi

apply(df_p, 1, function(x) {
  
  q <- 0
  
  sum(x^q)^(1/(1-q))
  
} )

x <- c(0, 1)
q <- 0

sum(x^q)^(1/(1-q))

hill_multifunc(df, c("F1", "F2"), scale = 0, stand_method = "none")
hill_multifunc


