
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

# first, we transpose the data matrix because we want to examine similarity among functions
func_matrix_inv <- as.data.frame(t((func_matrix[ , ])))

# create a distance matrix of the functions
func_dist <- dist(func_matrix_inv, method = "euclidean")

# we can plot these in a dendrogram

# start the function

# can use any method except: "kmeans" as this requires all the data

# index = "frey" doesn't work for some reason


adf <- func_matrix

vars <- names(func_matrix)

ind <- c("mcclain", "cindex", "silhouette", "dunn")

met <- "ward.D2"

dis <- "euclidean"


# get functions from the adf matrix
adf_mat <- adf[, vars]

# transpose the adf_mat data
adf_mat_t <- as.data.frame(t((adf_mat[ , ])))

# generate a dissimilarity matrix from adf_mat_t
adf_dist <- dist(adf_mat_t, method = dis)

# use the dissimilarity matrix to determine optimum group number for cluster analysis

opt_clus <- vector(length = length(ind))
  
for(i in 1:length(ind)) {
  
  x <- NbClust(diss = adf_dist, 
               distance = NULL, 
               min.nc = 2, max.nc = (ncol(adf_mat) - 1), 
               method = met, index = ind[i])
  
  opt_clus[i] <- x$Best.nc[1]
  
}

clus_n <- round(x = median(opt_clus), digits = 0)

adf_den <- hclust(func_dist, method = met )

adf_groups <- cutree(adf_den, k = clus_n)
adf_groups

# function loadings
func_loads <- 
  unlist(lapply(adf_groups, FUN = function(x) { ( (1)/length(adf_groups[adf_groups == x]) ) } ))













