
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Prepare the Jena data (specifically the 2004 data)

# load the relevant libraries
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# read in the function data from 2004
jena.func <- read_delim(here("data/41559_2017_391_MOESM5_ESM.csv"), delim = ";")
head(jena.func)

# rename the Plot column to plotcode
jena.func <- 
  jena.func %>%
  rename(plotcode = Plot)

# get nine reasonable functions for illustrative purposes
jena.func <- 
  jena.func %>%
  select(plotcode,
         c("BM_targ_DW", "Soil_NH4_AfterGrowth", "Soil_NO3_AfterGrowth",
           "Vol_coarseroot", "N_soil.insect.larv", "N_seedb_targ", "LAI",
           "Growth_root", "BM_microbes"))


# read in the Jena community composition data
jena.com <- read_csv(here("data/Jena_community_02-08.csv"))
head(jena.com)

# subset out the 2004 data
jena.com <- 
  jena.com %>%
  filter(year == 2004)

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
  mutate(across(.cols = spp, ~if_else(is.na(.), 0, .))) %>%
  group_by(year, sowndiv, plotcode) %>%
  summarise(across(.cols = spp, ~mean(.)), .groups = "drop") %>%
  arrange(year, sowndiv, plotcode)
head(jena.com)

# set-up a species by site matrix
spp.site <- 
  jena.com %>%
  select(spp)

# subset out the id vars
site <- 
  jena.com %>%
  select(-spp)

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

View(jena.dat)

# output this into a .csv file for further analysis
write_csv(x = jena.dat, file = here("data/jena_data_cleaned.csv"))




















