#'
#' @title Example of the workflow
#' 
#' @description Reanalyse data from Manning et al. (2018) to show how we might
#' operationalise some of the ideas. The data comes from FunDivEurope
#' 

# load the data
mf_df <- readr::read_delim("data/manning_2018_fundiv_data.txt", delim = "\t")
head(mf_df)

# check the column names
names(mf_df)

# remove the threshold columns
mf_df <-
  mf_df |>
  dplyr::select(!contains(match = "threshold"))

# get the relevant function columns
mf_df <-
  mf_df |>
  dplyr::select(PlotID, Latitude, Longitude, Ownership, Altitude, Slope,
                Exposition, Rocks.and.boulder.cover, Bedrock.type, 
                Calcareous.bedrock, Sand, Silt, Clay, Soil.drainage,
                Soil.depth, Soil.type, Stand.origin, Current.management,
                Forest.structure, Age.distribution, Age.canopy.trees,
                SpeciesRichness_All, TrueShannonIndex_All, TrueSimpsonsIndex_All,
                evergreen, coniferous, Country, Composition, Annual.Mean.T,
                Annual.Precip,
                productivity, biomass, wood_decomposition, litter_decomposition_all,
                total_soil_C, total_soil_N, WUE_average, root_biomass)

# get the country with the most data
mf_df |>
  dplyr::group_by(Country) |>
  dplyr::summarise(n = n())

# subset out Poland
mf_po <- 
  mf_df |>
  dplyr::filter(Country == "Poland")




