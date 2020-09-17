
# function to calculate the multifunctionlity index proposed by
# Manning et al 2018

# Manning, Peter, Fons Plas, Santiago Soliveres, Eric Allan, Fernando T Maestre,
# Georgina Mace, MARK J WHITTINGHAM, and Markus Fischer. 2018.
# “Redefining Ecosystem Multifunctionality.” 
# Nature Ecology & Evolution, February. Springer US, 1–10. 
# doi:10.1038/s41559-017-0461-7.

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector of functions to include which has 
# to correspond to column names

# Manning et al suggest to use cluster analysis to define cluster
# of functions that should be given equal weigt


# Function parameters:
# 
# cluster.method
# distance.metric
# Ellbow.detection.method


###### run to test function ######

source(here("Scripts","Multifunctionality-Simulations", "Multifunc_simulations_functions.R"))

specnum <- 9
funcnum <- 12
distribution = "rnorm"

FuncMat <- FunctionValue(specnum,funcnum, distribution, mean = 3, sd = 1)

vars <- func.names <- as.character( unique( FuncMat$Functions))

FuncMat[FuncMat$Functions %in% func.names[1:2],]$Funcval <- rnorm(specnum * 2, mean = 4, sd = 1)
FuncMat[FuncMat$Functions %in% func.names[3:4],]$Funcval <- rnorm(specnum * 2, mean = 5, sd = 1)

FuncMat %>% 
  spread(Functions, Funcval) %>% 
  select(-Species) %>% 
  cor %>% 
  corrplot

spec.names <- as.character( unique( FuncMat$Species))

maxrep <- 100 #using the full replications is prohibitive

SpecMat <- SpeciesMatrix(specnum = specnum, maxrep = maxrep)

method = "species_complementarity"

spec_comp <- SpecComp(specnum = specnum, funcnum = funcnum,
                      distribution = "rnorm", mean = 1, sd = 0.2,
                      spec_compfunc = func.names[1:3])

adf <- AvFunc <- AverageFunction(SpecMat, FuncMat,
                          method = method, 
                          spec_comp = spec_comp)

errM <- matrix(rnorm(n = nrow(AvFunc)*funcnum, mean = 0, sd = 0.01), ncol = funcnum)

#add variance
adf[,func.names] <- adf[,func.names] + errM

#####################################
# code from Manning:

#  Identify clustersof related functions
functions_matrix <-t(as.matrix(adf[adf$Richness == 1,vars]))
functions_matrix <-scale(functions_matrix)
d <-dist(functions_matrix, method = "euclidean")
dendrogram <-hclust(d, method = "complete" )
plot(dendrogram, cex = 0.6, hang = -1)

#use the dendrogram to determine the optimal number of clustersusing the Elbow method
factoextra::fviz_nbclust(functions_matrix, hcut, method = "wss", k.max = funcnum-1)


# In the final EF-multifunctionality,measure all clusters are weighted equally, 
# and this requires individual functions within a cluster to be given proportional weighing, 
# so that functions in a cluster sum to one. We coded this as:

loading_values <-c(0.167,1,0.333,0.333,1,0.167,0.333,0.167,0.167,0.167,0.167)


#recode the function values so that those which exceed the threshold are 
#assigned a value of 1, and those below are assigned a value of 0. 
mf_data_scaled <-matrix(nrow=nrow(mf_data),ncol=length(function_nrs2))
threshold <-0.5

for(i in 1:length(function_nrs2)){
  maximum <-mean(sort(mf_data[,function_nrs2[i]], decreasing=T)[c(1:10)])
  mf_data_scaled[,i] <-rep(0,length(mf_data_scaled[,i])) 
  high_values <-which(mf_data[,function_nrs2[i]] >= threshold * maximum) 
  mf_data_scaled[high_values,i] <-loading_values[i]}


# we sum the number of functions exceeding the 50% threshold and divide this 
# by the number of clustersand so maximum possible score

mf_data$multifunctionality<-rowSums(mf_data_scaled) / 7


######### code this into one function with automated cluster selection

hclust_multifunc <- function(adf, vars, scale = 1, HILL = TRUE){
  
  if(! "vegan" %in% installed.packages()[,1]) stop(
    "this function requires vegan to be installed"
  )
  
  if(length(scale) > 1) stop(
    "this function requires to choose a single order of diversity"
  )
  
  adf_mat <- adf[,vars]
  effN <- vegan::renyi(adf_mat, scales = scale, hill = TRUE)
  if(!HILL) effN <- effN/length(vars)
  meanFunc <- rowMeans(adf_mat) 
  
  adf$multifunc_effN <-  effN*meanFunc
  return(adf)
  
}
