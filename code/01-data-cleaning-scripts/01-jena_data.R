#'
#' @title Clean the raw Jena data
#' 
#' @description Load the Jena data taken from Jochum et al. (2020, Nature 
#' Ecology and Evolution), clean it and output a cleaned version for 
#' further analyses
#'

# load relevant libraries
library(dplyr)
library(readr)

# certain libraries must also be installed
if(! "vegan" %in% installed.packages()[,1]) print(
  "this script requires vegan to be installed"
)

# read in the Jena function data from Jochum et al. (2020)
jena_func <- readr::read_csv("data/J_avg_20190514.csv")
head(jena_func)
names(jena_func)

# rename the plot column to plotcode
jena_func <- 
  jena_func |>
  dplyr::rename(plotcode = plot)

# read in the Jena community composition data
jena_com <- readr::read_csv("data/Jena_Community_02-08.csv")
head(jena_com)
names(jena_com)

# how many years are there?
length(unique(jena_com$year))
unique(jena_com$year)

# check if the plot identities match-up between the two datasets: they do not
unique(jena_func$plotcode)
unique(jena_com$plotcode)

# why do they not match up?
# the jena_com data includes control plots that were not sown with anything
absent_plots <- unique(jena_com$plotcode)[which( !(unique(jena_com$plotcode) %in% jena_func$plotcode) )]

# check these plots that are abasent from the jena_com data
jena_com |>
  dplyr::filter(plotcode %in% absent_plots) |>
  View()

# get the plots used in the function data
jena_com <- 
  jena_com |>
  dplyr::filter(plotcode %in% unique(jena_func$plotcode) )

# test if the subsetting worked: if FALSE, all plotcodes are equal
any(sort(unique(jena_com$plotcode)) != sort(jena_func$plotcode) )

# subset out the relevant columns

# get a list of species names
col_names <- names(jena_com)
spp <- col_names[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", col_names)]

# there should be 60 species in the dataset
length(spp) == 60

# summarise the percentage cover values by species across years
jena_com |>
  dplyr::select(dplyr::all_of(spp)) |>
  View()

jena_com <- 
  jena_com |>
  dplyr::group_by(year, sowndiv, plotcode) |>
  dplyr::summarise( dplyr::across(.cols = dplyr::all_of(spp), 
                                  ~mean(., na.rm = TRUE)), .groups = "drop") |>
  dplyr::arrange(year, sowndiv, plotcode) |>
  dplyr::group_by(sowndiv, plotcode) |>
  dplyr::summarise( dplyr::across(.cols = dplyr::all_of(spp), 
                                  ~mean(., na.rm = TRUE)), .groups = "drop")

# convert the NAs to 0's to indicate absence
jena_com <- 
  jena_com |>
  dplyr::mutate( dplyr::across(.cols = dplyr::all_of(spp), 
                               ~dplyr::if_else(is.na(.), 0, . )))

# join these data to the jena_func data and sort
jena_func <- 
  jena_func |>
  dplyr::arrange(sowndiv, plotcode)

# check if jena.com and jena.func match: if FALSE then they match
any(jena_com$plotcode != jena_func$plotcode)

# does the realised diversity calculated from jena_com match jena_func
# yes, Pearson r is 0.99 so the community data are likely a good representation
plot(jena_func$S, rowSums(vegan::decostand(jena_com[, spp], method = "pa")) )
cor.test(jena_func$S, rowSums(vegan::decostand(jena_com[, spp], method = "pa")) )

# add the calculated realised diversity to the jena.func data
jena_func$realised_richness2 <- rowSums(vegan::decostand(jena_com[, spp], method = "pa"))

# join the species abundance data to the jena.func data
jena_all <- dplyr::full_join(jena_func, jena_com, by = c("sowndiv", "plotcode"))

# rename the S column
jena_all <- 
  jena_all |>
  rename(realised_richness1 = S)

# reorder the columns
names(jena_all)

jena_all <- 
  jena_all |>
  dplyr::select(plotcode, block, sowndiv, realised_richness1, realised_richness2,
                dplyr::all_of(spp),
                biomass, plantCN, soilC, soilorgC, herbi, micBMC, phoact, poll, rootBM)

# remove the 60 species plots
jena_all <- 
  jena_all |>
  filter(sowndiv < 60)

length(names(jena_all)[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", names(jena_all))])

# check for NAs
summary(jena_all)
func_names <- c("biomass", "plantCN", "soilC", "soilorgC", "herbi", "micBMC", "phoact", "poll","rootBM")

# get only the complete cases
# we lose three plots
jena_all <- jena_all[complete.cases(jena_all[, func_names]), ]

# output this cleaned data file
saveRDS(object = jena_all, file = "data/jena_data.rds")


# prepare the diversity-interactions model dataset
sp_comp <- readr::read_delim(file = "data/plotinfo.csv", delim = ";")
head(sp_comp)

# get the relevant plotcodes
sp_comp <- 
  sp_comp |>
  dplyr::filter(plotcode %in% unique(jena_all$plotcode) )

# check that the number of rows match
nrow(sp_comp) == nrow(jena_all)

sp_df <- 
  sapply(sp_comp$composition, function(x) {
  
  strsplit(x, split = "[|]")
  
}, USE.NAMES = FALSE)

# convert into a species list
sp_list <- unique(unlist(sp_df))
any( sort(sp_list) != sort(spp) )

# generate an output data.frame
df <- data.frame(plotcode = sp_comp$plotcode)
for (i in 1:length(sp_list)) {
  
  z <- 
    lapply(sp_df, function(x) {
      
      if (sp_list[i] %in% x) {
        
        y <- 1/length(x)
        
      } else {
        
        y <- 0
        
      }
      
      y
      
    })
  
  df[, sp_list[i]] <- unlist(z)
  
}

any(rowSums(df[, -1]) != 1)

# join this dataset to a reduced version of jena.all
jena_di <- dplyr::full_join( dplyr::select(jena_all, -dplyr::all_of(spp) ),
                             df, 
                             by = "plotcode")

# check the jena_di dataset
names(jena_di)
head(jena_di)

# reorder the columns to match Laura RPubs
jena_di <- 
  jena_di |>
  dplyr::select(plotcode, block, sowndiv, realised_richness1, realised_richness2,
                dplyr::all_of(spp), 
                dplyr::all_of(func_names))

# export the di format data
saveRDS(object = jena_di, file = "data/jena_data_di.rds")

### END
