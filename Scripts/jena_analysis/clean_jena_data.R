
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code to clean the Jena data for analysis

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)

# certain libraries must also be installed
if(! "vegan" %in% installed.packages()[,1]) print(
  "this script requires vegan to be installed"
)

# read in the Jena function data from Jochum et al. (2020)
jena.func <- read_csv(here("data/J_avg_20190514.csv"))
head(jena.func)
names(jena.func)

# rename the plot column to plotcode
jena.func <- 
  jena.func %>%
  rename(plotcode = plot)

# read in the Jena community composition data
jena.com <- read_csv(here("data/Jena_community_02-08.csv"))
head(jena.com)
names(jena.com)

# how many years are there?
length(unique(jena.com$year))
unique(jena.com$year)

# check if the plot identities match-up between the two datasets
# they do not
unique(jena.func$plotcode)
unique(jena.com$plotcode)

# why do they not match up?
# the jena.com data includes control plots that were not sown with anything
absent_plots <- unique(jena.com$plotcode)[which( !(unique(jena.com$plotcode) %in% jena.func$plotcode) )]

jena.com %>%
  filter(plotcode %in% absent_plots) %>%
  View()
rm(absent_plots)

# get the plots used in the function data
jena.com <- 
  jena.com %>%
  filter(plotcode %in% unique(jena.func$plotcode) )

# test if the subsetting worked: if FALSE, all plotcodes are equal
any(sort(unique(jena.com$plotcode)) != sort(jena.func$plotcode) )

# subset out the relevant columns

# get a list of species names
col_names <- names(jena.com)
spp <- col_names[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", col_names)]
rm(col_names)

# there should be 60 species in the dataset
length(spp) == 60

# summarise the percentage cover values by species across years
jena.com %>%
  select(all_of(spp)) %>%
  View()

jena.com <- 
  jena.com %>%
  group_by(year, sowndiv, plotcode) %>%
  summarise(across(.cols = all_of(spp), ~mean(., na.rm = TRUE)), .groups = "drop") %>%
  arrange(year, sowndiv, plotcode) %>%
  group_by(sowndiv, plotcode) %>%
  summarise(across(.cols = all_of(spp), ~mean(., na.rm = TRUE)), .groups = "drop")

# convert the NAs to 0's to indicate absence
jena.com <- 
  jena.com %>%
  mutate(across(.cols = all_of(spp), ~if_else(is.na(.), 0, . )))

# join these data to the jena.func data
# sort these data
jena.func <- 
  jena.func %>%
  arrange(sowndiv, plotcode)

# check if jena.com and jena.func match: if FALSE then they match
any(jena.com$plotcode != jena.func$plotcode)

# does the realised diversity calculated from jena.com match jena.func
# yes, Pearson r is 0.99 so the community data are likely a good representatino
plot(jena.func$S, rowSums(vegan::decostand(jena.com[, spp], method = "pa")) )
cor.test(jena.func$S, rowSums(vegan::decostand(jena.com[, spp], method = "pa")) )

# add the calculated realised diversity to the jena.func data
jena.func$realised_richness2 <- rowSums(vegan::decostand(jena.com[, spp], method = "pa"))

# join the species abundance data to the jena.func data
jena.all <- full_join(jena.func, jena.com, by = c("sowndiv", "plotcode"))

# rename the S column
jena.all <- 
  jena.all %>%
  rename(realised_richness1 = S)

# reorder the columns
names(jena.all)

jena.all <- 
  jena.all %>%
  select(plotcode, block, sowndiv, realised_richness1, realised_richness2,
         all_of(spp),
         biomass, plantCN, soilC, soilorgC, herbi, micBMC, phoact, poll, rootBM)

# remove the 60 species plots
jena.all <- 
  jena.all %>%
  filter(sowndiv < 60)

length(names(jena.all)[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", names(jena.all))])

# check for NAs
summary(jena.all)
func.names <- c("biomass", "plantCN", "soilC", "soilorgC", "herbi", "micBMC", "phoact", "poll","rootBM")

# get only the complete cases
# we lose three plots
jena.all <- jena.all[complete.cases(jena.all[, func.names]), ]

# output this cleaned data file
write_csv(x = jena.all, file = here("data/jena_data_Jochum_2020_clean.csv"))


# prepare the dataset for Laura
sp_comp <- read_delim(file = here("data/plotinfo.csv"), delim = ";")
head(sp_comp)

# get the relevant plotcodes
sp_comp <- 
  sp_comp %>%
  filter(plotcode %in% unique(jena.all$plotcode) )

nrow(sp_comp) == nrow(jena.all)

sp.df <- 
  sapply(sp_comp$composition, function(x) {
  
  strsplit(x, split = "[|]")
  
}, USE.NAMES = FALSE)

# convert into a species list
sp.list <- unique(unlist(sp.df))
any( sort(sp.list) != sort(spp) )

df <- data.frame(plotcode = sp_comp$plotcode)

for (i in 1:length(sp.list)) {
  
  z <- 
    lapply(sp.df, function(x) {
      
      if (sp.list[i] %in% x) {
        
        y <- 1/length(x)
        
      } else {
        
        y <- 0
        
      }
      
      y
      
    })
  
  df[, sp.list[i]] <- unlist(z)
  
}

any(rowSums(df[, -1]) != 1)

# join this dataset to a reduced version of jena.all
jena.di <- 
  full_join(jena.all %>%
              select(-all_of(spp)),
            df, by = "plotcode")

names(jena.di)

# reorder the columns to match Laura RPubs
jena.di <- 
  jena.di %>%
  select(plotcode, block, sowndiv, realised_richness1, realised_richness2,
         all_of(spp), 
         all_of(func.names))

write_csv(x = jena.di, file = here("data/jena_data_Jochum_2020_clean_di.csv"))

### END
