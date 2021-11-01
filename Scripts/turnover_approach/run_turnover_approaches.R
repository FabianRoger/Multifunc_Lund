
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Run the turnover approaches on the Jena data

# load relevant libraries
library(here)
library(readr)
library(dplyr)
library(ggplot2)

# link to scripts with the relevant functions
source(here("Scripts/turnover_approach/turnover_approach_functions.R"))
source(here("Scripts/function_plotting_theme.R"))
source(here("Scripts/process_model_data.R"))

# read in the Jena data
jena.dat <- read_csv(file = here("data/jena_data_Jochum_2020_clean.csv"))

# get a vector of species names
col_names <- names(jena.dat)

spp <- col_names[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", col_names)]
length(spp)
rm(col_names)

# check if all species are present in at least one plot
any(colSums(jena.dat[, spp]) == 0)
min(colSums(jena.dat[, spp]))

# convert the species abundances into presence-absence data
jena.dat <- 
  jena.dat %>%
  mutate(across(.cols = all_of(spp), ~if_else(. > 0, 1, 0)))

# get a vector of function names
names(jena.dat)
func.names <- c("biomass", "plantCN", "soilC", "soilorgC", "herbi",
                "micBMC", "phoact", "poll", "rootBM")



# test the AIC_sp function
df2 <- SES_score(data = jena.dat, function_names = func.names, species_names = spp, 
          n_ran = 10)

df2

df <- AIC_sp(data = jena.dat, function_names = func.names[2], species_names = spp)

x <- apply(df2[, -1], 1, function(x) { ifelse(any(x > 0), TRUE, FALSE) })

df[x, ]$species

# for each dataset in mod.list and for each function matrix in func.list, run the AIC and SES turnover approaches

# set the number of reps for the SES method
ses.r <- 100

mod.turnover.test <- vector("list", length = length(mod.list))
for (i in 1:length(mod.list)) {
  
  func.reps <- vector("list", length = length(func.list))
  for (j in 1:length(func.list)) {
    
    # process the simulated cluster using: process_sim_data
    df.proc <- 
      process_sim_data(model_data = mod.list[[i]], 
                       func.mat =  func.list[[j]], 
                       time_final = TRUE, 
                       species_abun = "pa")
    
    # run the different turnover approaches
    func.reps[[j]] <- turnover_tester(model_dat = df.proc, function_matrix = func.list[[j]], ses_reps = ses.r) 
    
  }
  
  mod.turnover.test[[i]] <- bind_rows(func.reps, .id = "function_matrix")
  
}

# assign names to this object
names(mod.turnover.test) <- names(mod.list)

# bind this into one large data.frame
mod.turnover.df <- bind_rows(mod.turnover.test, .id = "parameter_combination")

# pull this into a longer data.frame
library(ggplot2)

df.plot <- 
  mod.turnover.df %>%
  pivot_longer(cols = starts_with("F_"),
               names_to = "function_id",
               values_to = "metric") %>%
  mutate(mod_id = paste(parameter_combination, function_matrix, sep = "_"), .before = 1)

View(df.plot)

# df.plot$mod_id <- as.factor(df.plot$mod_id)
# levels(df.plot$mod_id) <- paste("m.", 1:length(unique(df.plot$mod_id)), sep = "")
# df.plot$mod_id <- as.character(df.plot$mod_id)

df.plot %>%
  group_by(mod_id, model_run, method, output_metric) %>%
  summarise(mean_metric = mean(metric, na.rm = TRUE), .groups = "drop") %>%
  filter(output_metric == "prop_incorrect") %>%
  ggplot(data = .,
         mapping = aes(x = method, y = mean_metric)) +
  geom_jitter() +
  facet_wrap(~mod_id)

df.plot %>%
  ggplot(data = ., 
       mapping = aes(x = method, y = metric, colour = mod_id)) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic() +
  facet_wrap(~output_metric, scales = "free") +
  scale_colour_viridis_d() +
  theme(legend.position = "bottom",
        axis.text = element_text(colour = "black"))

### END
