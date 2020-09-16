
# Manning et al. (2018): Clustering-based multifunctionality

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(NbClust)

source(here("Scripts","Multifunctionality-Simulations", "Multifunc_simulations_functions.R"))
source(here("Scripts", "MF_index_function.R"))
source(here("Scripts", "MF_Hill_function.R"))

### Simulate full diversity experiment
set.seed(777)

specnum <- 10
funcnum <- 10

distribution = "runif"

FuncMat <- FunctionValue(specnum,funcnum, distribution, min = 0.1, max = 0.9)

func.names <- as.character( unique( FuncMat$Functions))
spec.names <- as.character( unique( FuncMat$Species))


maxrep <- 100 # using the full replications is prohibitive

SpecMat <- SpeciesMatrix(specnum = specnum, maxrep = maxrep)

method = "av"

compfunc <- func.names[1:3]

AvFunc <- AverageFunction(SpecMat, FuncMat,
                          method = method, 
                          compfunc = compfunc)

set.seed(563)
errM <- matrix(rnorm(n = nrow(AvFunc)*funcnum, mean = 0, sd = 0.01), ncol = funcnum)

#add variance
AvFunc[,func.names] <- AvFunc[,func.names] + errM

# standardize functions 
AvFunc_func <- 
  AvFunc %>% 
  mutate_at(vars(one_of(func.names)), function(x) {(x) / max(x)})


# get a matrix of functions
func_matrix <- AvFunc_func[ , func.names]

# cluster analysis of the functions
func_matrix_inv <- as.data.frame(t((func_matrix[ , ])))

d <- dist(func_matrix_inv, method = "euclidean")

dendrogram <- hclust(d, method = "complete" )

grp <- cutree(dendrogram, k = 4)
grp

NbClust(diss = d, distance = NULL, min.nc = 2, max.nc = 9, 
        method = "complete", index = "mcclain")

NbClust(diss = d, distance = NULL, min.nc = 2, max.nc = 9, 
        method = "ward.D2", index = "cindex")

NbClust(diss = d, distance = NULL, min.nc = 2, max.nc = 9, 
        method = "ward.D2", index = "silhouette")

NbClust(diss = d, distance = NULL, min.nc = 2, max.nc = 9, 
        method = "ward.D2", index = "dunn")










