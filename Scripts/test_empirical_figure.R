
# empirical exploration figure test

# load the multifunctionality calculations
source("Scripts/01_general_functions/MF_functions_collated.R")
source("Scripts/01_general_functions/function_plotting_theme.R")

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# make a test data.frame
df <- tibble(Com = c(rep(1, 5),
                     rep(2, 5),
                     rep(3, 5),
                     rep(4, 5),
                     rep(5, 5),
                     rep(6, 5)),
             Func = rep(paste0("F", 1:5), 6),
             Func_val = c(0, 0, 0, 0, 0,
                          0.4, 0, 0, 0, 0,
                          0.2, 0.2, 0.1, 0.2, 0.1,
                          0, 0, 0.1, 0.9, 0,
                          0.3, 0.5, 0.7, 0.4, 0.4,
                          1, 1, 1, 1, 1))
head(df)

# convert the Com column into a character
df$Com <- as.character(df$Com)

# plot the function distributions
p1 <- 
  ggplot(data = df, 
       mapping = aes(x = Func, y = Func_val, fill = Com)) +
  geom_col(width = 0.5) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  facet_wrap(~Com, ncol = 6, nrow = 1) +
  ylab("Stand. function value (0-1)") +
  xlab(NULL) +
  theme_meta() +
  theme(legend.position = "none")

ggsave(filename = "Figures/test_metric1.png", p1, dpi = 400,
       units = "cm", width = 20, height = 5)  

# calculate the different metrics

# write a function with the relevant metrics
calculate_MF <- function(data, func.names) {
  
  data %>%
    mutate(`scal. MF` = MF_jing(adf = data, vars = func.names),
           `sum MF` = MF_sum(adf = data, vars = func.names, stand_method = "none"),
           `ave. MF` = MF_av(adf = data, vars = func.names, stand_method = "none"),
           `MESLI MF` = MF_mesli(adf = data, vars = func.names, stand_method = "none"),
           `Pasari MF` = MF_pasari(adf = data, vars = func.names, stand_method = "none"),
           `SAM MF` = MF_dooley(adf = data, vars = func.names,  stand_method = "none"),
           `Simp. MF` = MF_simpsons_div(adf = data, vars = func.names, stand_method = "none"),
           `Shannon MF` = MF_shannon_div(adf = data, vars = func.names, stand_method = "none"),
           `ENF.Q0 MF` = hill_multifunc(adf = data, vars = func.names, scale = 1, HILL = TRUE, stand_method = "none"),
           `ENF.Q1 MF` = hill_multifunc(adf = data, vars = func.names, scale = 2, HILL = TRUE, stand_method = "none"),
           `Cluster.30 MF` = manning_multifunc(adf = data, vars = func.names, thresh = 0.3),
           `Cluster.70 MF` = manning_multifunc(adf = data, vars = func.names, thresh = 0.7),
           `thresh.30 MF` = single_threshold_mf(adf = data, vars = func.names, thresh = 0.3),
           `thresh.70 MF` = single_threshold_mf(adf = data, vars = func.names, thresh = 0.7),
           `PCA MF` = pca_multifunc(adf = data, vars = func.names)
    )
  
}

# process the data into the correct format
df_w <- 
  df %>%
  pivot_wider(id_cols = "Com", 
              names_from = "Func", 
              values_from = "Func_val")

# calculate the different multifunctionality metrics
df_mf <- calculate_MF(data = df_w, func.names = paste0("F", 1:5))
df_mf <- 
  df_mf %>%
  select(Com, contains("MF"))

# convert the data to the long format
MF_names <- names(df_mf)[-1]
df_mf <- 
  df_mf %>%
  pivot_longer(cols = all_of(MF_names),
               names_to = "MF_metric",
               values_to = "MF")  %>%
  mutate(MF_metric = factor(MF_metric, levels = MF_names)) %>%
  arrange(MF_metric, Com) %>%
  mutate(MF = ifelse(is.infinite(MF) | is.na(MF), NA, MF))

df_undef <- 
  df_mf %>%
  filter(is.na(MF))
df_undef$MF <- "NA"
df_undef$y <- c(1, 0.5, 0.5)

# plot the results
p2 <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_line(data = df_mf,
            mapping = aes(x = as.integer(Com), y = MF),
            colour = "black", alpha = 0.3) +
  geom_point(data = df_mf,
             mapping = aes(x = as.integer(Com), y = MF, colour = Com),
             size = 2.5) +
  geom_text(data = df_undef,
            mapping = aes(x = as.integer(Com), y = y, label = MF),
            size = 2.5) +
  facet_wrap(~MF_metric, scales = "free_y", ncol = 5, nrow = 3) +
  theme_meta() +
  xlab("Community (1-6)") +
  ylab("Multifunctionality") +
  scale_x_continuous(breaks = c(1:6), limits = c(0.5, 6.5)) +
  theme(legend.position = "none")

ggsave(filename = "Figures/test_metric2.png", p2, dpi = 400,
       units = "cm", width = 20, height = 14)  



