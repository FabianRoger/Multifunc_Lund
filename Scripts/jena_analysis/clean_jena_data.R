
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Code to clean the Jena data for analysis

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

library(here)

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

# get the plots used in the function data
jena.com <- 
  jena.com %>%
  filter(plotcode %in% unique(jena.func$plotcode) )

# subset out the relevant columns
names(jena.com)
spp <- names(jena.com)[52:100]

# replace the NAs with zeros then summarise abundance across seasons
jena.com <- 
  jena.com %>%
  mutate(across(.cols = all_of(spp), ~if_else(is.na(.), 0, .))) %>%
  group_by(year, sowndiv, plotcode) %>%
  summarise(across(.cols = spp, ~mean(.)), .groups = "drop") %>%
  arrange(year, sowndiv, plotcode)
head(jena.com)

# set-up a species by site matrix
spp.site <- 
  jena.com %>%
  select(all_of(spp))

# subset out the id vars
site <- 
  jena.com %>%
  select(-all_of(spp) )

# convert the spp.site to relative abundance
spp.site <- spp.site/if_else(rowSums(spp.site) == 0, 1, rowSums(spp.site))

# add realised diversity to the site data
site$realised_diversity <- rowSums( vegan::decostand(spp.site, method = "pa") )

# add the species abundances back to the site
site <- bind_cols(site, spp.site)

# add the function data to the site data
jena.dat <- 
  full_join(site,
            jena.func,
            by = "plotcode")

# remove plots with zero realised diversity and the 60 species plots
jena.dat <- 
  jena.dat %>%
  filter(realised_diversity > 0) %>%
  filter(sowndiv < 60)
head(jena.dat)

# reflect the soil nutrient functions (i.e. low nutrients means high uptake)
jena.dat <- 
  jena.dat %>%
  mutate(across(.cols = starts_with("Soil"), ~(.*-1) ))

# view the dataset
View(jena.dat)