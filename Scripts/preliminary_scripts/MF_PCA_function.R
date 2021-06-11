
# function to calculate the pca based multifunctionality index as suggested by 
# Meyer et al. (2017) Biodiversity-multifunctionality relationships depend on identity and number of measured functions #
# Nature Ecology & Evolution   

# adf, is dataframe with plots in rows, and functions in columns
# vars has to bee a named vector 

pca_multifunc <- function(adf, vars){
  
  if(! "vegan" %in% installed.packages()[,1]) stop(
    "this function requires vegan to be installed"
  )
  
  adf_mat <- adf[,vars]
  pca2<-vegan::rda(adf_mat)
  
  temp2 <- vegan::scores(pca2, choices=1:length(vars), display=c("sites")) 
  eig<-summary(pca2)$cont$importance[1,]
  for(i in 1:length(eig)) temp2[,i] <- temp2[,i] * eig[i]
  Index.wt <- rowSums(temp2)
  adf$multifunc_pca_ind <-  Index.wt
  return(adf)
  
}
