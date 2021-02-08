
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Does the number of functions matter?

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)

# load the cleaned Jena data
jena.dat <- read_csv(file = here("data/jena_data_cleaned.csv"))

# define variable groups
var.names <- names(jena.dat)

# (1) get species names
spp.p <- ( grepl("+\\.+", var.names) & nchar(var.names) == 7 )
spp.names <- var.names[spp.p]

# (2) get site identifiers
site.id <- c("year", "sowndiv", "plotcode", "realised_diversity")

# (3) get function names
func.names <- var.names[!(var.names %in% ( c(spp.names, site.id) ))  ]


# BEF slope and n-functions

n.funcs <- 1:length(func.names)

function.combinations <- vector("list", length = length(func.names))
for (i in n.funcs){
  function.combinations[[i]] <- combn(x = func.names, m = i)
}
function.combinations[[9]]


