###########################################################################################################################
#                                                                                                                         #
# R-Code to updating the turnover approach to multifunctionality to included more stringent criteria for considering an	  #
# effect of an individual species informative. This allows to compensate for a high number of tests when many functions   #
# and many species are included in an analysis using the turnover approach. The code also includes analyses of simulated  #
# data to explore, based on permutation procedures, how large the potential bias in the results of the turnover approach  #
# are when using very large data sets.                       								  # 
#                                                                                                                         #
# Suppplementary material to 												  #
# Meyer et al. (2017) Biodiversity-multifunctionality relationships depend on identity and number of measured functions   #
# Nature Ecology & Evolution                                                                                              #
#                                                                                                                         #
###########################################################################################################################

###############################################################
#                                                             #
#   Define the functions to calculate the turnover approach   #
#                                                             #
###############################################################

#The functions from the multifunc package are used to calculate the turnover approach
install_github("multifunc", username="jebyrnes")
library(multifunc)

#One function from the multifunc package needs to be adapted to work with large numbers of functions when testing ALL combinations
#of functions becomes impossible
#get the function divNeeded to work...
divNeeded2 <- function(overData, type) {
  overData <- filterOverData2(overData, type = type) #Ignore for the moment as there are no negative effects
  dn <- sapply(1:nrow(overData), function(m) { #nrow(overData)
    print(m)
    if (m == 1) {
      rowSums(overData)
    }
    else {
      b <- NA
      for (i in 1:1000) {
        a <- length(which(colSums(overData[sample(1:nrow(overData), m, replace = FALSE, prob = NULL),])  > 0))
        b <- c(b, a)
      }
      b <- b[-1]                                                                                                    
      #combn(1:nrow(overData), m, function(x) length(which(colSums(overData[x, ]) > 0)))
    }
  })
  ret <- data.frame(nfunc = 1, div = dn[[1]])
  for (i in 2:length(dn)) ret <- rbind(ret, data.frame(nfunc = i, div = dn[[i]]))
  return(ret)
}

###############################################################################
# Update functions from the multifunc package to work with different k-values #
###############################################################################
#sAICFit does the business of fitting a model using a stepAIC approach
########
sAICFit2 <- function(response, species, data, method="lm", combine="+", k, ...)
{
  f <- as.formula(paste(response, "~" , paste(species, collapse="+")))
  fit <- eval(substitute(lm(f, data=data, ...)))
  obj <- stepAIC(fit, trace=0, k=k)
  obj
}

#sAICfun takes a dataset, response, and function, and then uses a stepAIC approach
#to determine the best model.  From that it extracts the species with a positive, 
#negative, and neutral effect on that function
#
#In this new defined version 2 an option to specify disired cutoff value as k hase been implemented
#########
sAICfun2 <- function(response, species, data, positive.desired=T, method="lm", combine="+", k=2, ...)
{
  #first fit the model
  obj<-sAICFit2(response, species, data, method, combine, k=k, ...)
  
  #now extract the important information about positive, negative, etc.
  
  #return that info in a list
  if(positive.desired) {
    pos.sp <- names(summary(obj)[[4]][,4][(summary(obj)[[4]][,1]>0) & names(summary(obj)[[4]][,4])!="(Intercept)"])
    neg.sp <- names(summary(obj)[[4]][,4][(summary(obj)[[4]][,1]<0) & names(summary(obj)[[4]][,4])!="(Intercept)"])
  }else{
    pos.sp <- names(summary(obj)[[4]][,4][(summary(obj)[[4]][,1]<0) & names(summary(obj)[[4]][,4])!="(Intercept)"])
    neg.sp <- names(summary(obj)[[4]][,4][(summary(obj)[[4]][,1]>0) & names(summary(obj)[[4]][,4])!="(Intercept)"])
  }
  neu.sp <- species[!(species %in% pos.sp) & !(species %in% neg.sp)]
  
  #make a vector of 1s and 0s
  effects<-rep(0, length(species))
  names(effects)<-species  
  effects[which(names(effects) %in% pos.sp)]<-1
  effects[which(names(effects) %in% neg.sp)]<- -1
  
  coefs<-c(effects,0)
  names(coefs)[length(coefs)]<-"(Intercept)"
  coefs[ match(names(coef(obj)), names(coefs)) ]<-coef(obj)
  
  return(list(pos.sp=pos.sp,neg.sp=neg.sp,neu.sp=neu.sp, functions=response, coefs=coefs, effects=effects))
}

getRedundancy2 <- function (vars, species, data, negVars = NA, method = "lm", combine = "+", k, output = "effect", ...) 
{
  res.list <- lapply(vars, function(x) {
    positive.desired <- T
    if (x %in% negVars) 
      positive.desired <- F
    sAICfun2(response = x, species = species, data = data, 
             positive.desired, method = method, combine = combine, k=k,
             ...)
  })
  if (output == "coef") {
    ret <- ldply(res.list, function(x) x$coefs)
  }
  else {
    ret <- ldply(res.list, function(x) x$effects)
  }
  rownames(ret) <- vars
  return(ret)
}

filterOverData2 <- function (overData, type = "positive") 
{
  neg <- which(overData < 0, arr.ind = T)
  pos <- which(overData > 0, arr.ind = T)
  if (type == "positive" & nrow(neg)>0) 
    apply(neg, 1, function(x) overData[x[1], x[2]] <<- 0)
  if (type == "negative" & nrow(pos)>0) 
    apply(pos, 1, function(x) overData[x[1], x[2]] <<- 0)
  if (type == "negative" & nrow(neg)>0) 
    apply(neg, 1, function(x) overData[x[1], x[2]] <<- 1)
  if (type == "all" & nrow(neg)>0) 
    apply(neg, 1, function(x) overData[x[1], x[2]] <<- 1)
  overData
}

########################################################
#                                                      #
#    Simulate the data to test the turnover approach   #
#                                                      #
########################################################

#Defining the function to performe the simulations of data
#N_plots: The number of plots in the simulated data. Needs to be a multiple of 4
#N_species: The number of species in the simulated species pool
#N_funct: The number of simulated functions
#PropPool: The proportion of the species pool that has an effect on functioning (= affecting species pool)
#PropSpecies: The proportion of the affecting species pool that has an effect on an individual function
#k: the cut of value to define what effects are considered informative
#SRatio: The signal to noise ratio, the random term will be devided by this number

#defining a functions that simulates the data and calculates the results pased on the functions for the turnover approach defined above
addPropLine <- function(N_plots=40, N_species=10, N_funct=20, PropPool=0.5, PropSpecies=0.2, k=2, SRatio=10)
{
  SmallSim <- as.data.frame(matrix(NA, nrow=N_plots, ncol=3+N_species+N_funct))
  names(SmallSim) <- c("plotcode", "div", "comp", paste("S", 1:N_species, sep=""), paste("F", 1:N_funct, sep=""))
  SmallSim$plotcode <- paste("P", 1:N_plots, sep="")
  SmallSim$div <- rep(c("1", "2", "4", "8"), each=N_plots/4)
  SmallSim[,(1+3):(N_species+3)] <- 0
  for (i in 1:N_plots) SmallSim[i, sample(1:N_species, SmallSim[i,"div"])+3] <- 1 #create information on presence of species
  
  #Define the affecting species pool
  AffectPool <- sample(1:N_species, PropPool*N_species, replace=FALSE)
  
  #Simulate the functions as the sum of presence/absance data of a given proportion of the species pool
  #with a random variable added that is scaled with a tenth of the number of species to keep relative size of the random
  #term independent of the size of species pool and the proportion of contributing species
  if (ceiling(PropSpecies*PropPool*N_species) > 1) 
  {
    for (i in 1:N_funct)
    SmallSim[,paste("F", i, sep="")] <- rowSums(SmallSim[,sample(paste("S", AffectPool, sep=""), ceiling(PropSpecies*PropPool*N_species), replace=TRUE)]) + rnorm(N_plots)*N_species*PropPool*PropSpecies/SRatio + rnorm(N_plots)/1000000 #Very small constant random number added so that the simulations work when no significant species effects occur 
  }
  if (ceiling(PropSpecies*PropPool*N_species) == 1) 
  {
    for (i in 1:N_funct)
    SmallSim[,paste("F", i, sep="")] <- SmallSim[,sample(paste("S", AffectPool, sep=""), ceiling(PropSpecies*PropPool*N_species), replace=TRUE)] + rnorm(N_plots)*N_species*PropPool*PropSpecies/SRatio + rnorm(N_plots)/1000000 #Very small constant random number added so that the simulations work when no significant species effects occur 
  }
  if (ceiling(PropSpecies*PropPool*N_species) == 0) 
  {
    for (i in 1:N_funct)
      SmallSim[,paste("F", i, sep="")] <- rnorm(N_plots)/1000000 #Very small constant random number added so that the simulations work when no significant species effects occur 
  }
  
  #Permutation of the data
  SmallSimRand <- SmallSim
  SmallSimRand[,(N_species+4):(N_species+N_funct+3)] <- as.data.frame(apply(SmallSimRand[,(N_species+4):(N_species+N_funct+3)], 2, function(x) sample(x, N_plots, replace=FALSE)))
  
  #calculate the turn-over approach for the simulated data with the functions from the multifunc package
  Data <- SmallSim
  vars <- names(Data[,(4+N_species):(3+N_species+N_funct)])
  SpeciesList <- names(Data[,4:(3+N_species)])
 
   
  res.list<-lapply(vars, function(x) sAICfun2(x, SpeciesList, Data, k=k))
  names(res.list)<-vars
  
  redund<-getRedundancy2(vars, SpeciesList, Data, k=k)
  coefs<-getRedundancy2(vars, SpeciesList, Data, output="coef", k=k)
  stdCoefs<-stdEffects(coefs, Data, vars, SpeciesList)
  
  #plot the number of functions by fraction of the species pool needed
  posCurve<-divNeeded2(redund, type="positive")
  posCurve2 <- posCurve
  posCurve2$div<-posCurve2$div/ncol(redund)
  
  #find number of functions where the mean proportion differes from 1
  M_nfunc <- aggregate(posCurve2$div, by=list("nfunc"=posCurve2$nfunc), FUN=mean)
  SD_nfunc <- aggregate(posCurve2$div, by=list("nfunc"=posCurve2$nfunc), FUN=sd)
  M_nfunc$Q <- NA 
  M_nfunc$Q <- qnorm(0.95, mean=M_nfunc$x, sd=SD_nfunc$x)
  
  #safe variables from normal data
  res.list_ref <- res.list 
  redund_ref <- redund 
  coefs_ref <- coefs 
  stdCoefs_ref <- stdCoefs 
  posCurve_ref <- posCurve 
  posCurve2_ref <- posCurve2 
  M_nfunc_ref <- M_nfunc
  
  #Calculate the turn-over approach for perumtated data with the functions from the multifunc package
  Data <- SmallSimRand
  vars <- names(Data[,(4+N_species):(3+N_species+N_funct)])
  SpeciesList <- names(Data[,4:(3+N_species)])
  
  res.list<-lapply(vars, function(x) sAICfun2(x, SpeciesList, Data, k=k))
  names(res.list)<-vars
  
  redund<-getRedundancy2(vars, SpeciesList, Data, k=k)
  coefs<-getRedundancy2(vars, SpeciesList, Data, output="coef", k=k)
  stdCoefs<-stdEffects(coefs, Data, vars, SpeciesList)
  
  #plot the number of functions by fraction of the species pool needed
  posCurve<-divNeeded2(redund, type="positive")
  posCurve2 <- posCurve
  posCurve2$div<-posCurve2$div/ncol(redund)
  
  #find number of functions where the mean proportion differes from 1
  M_nfunc <- aggregate(posCurve2$div, by=list("nfunc"=posCurve2$nfunc), FUN=mean)
  SD_nfunc <- aggregate(posCurve2$div, by=list("nfunc"=posCurve2$nfunc), FUN=sd)
  M_nfunc$Q <- NA 
  M_nfunc$Q <- qnorm(0.95, mean=M_nfunc$x, sd=SD_nfunc$x)
  
  #safe variables from randomized data
  res.list_rand <- res.list 
  redund_rand <- redund 
  coefs_rand <- coefs 
  stdCoefs_rand <- stdCoefs 
  posCurve_rand <- posCurve 
  posCurve2_rand <- posCurve2 
  M_nfunc_rand <- M_nfunc
  
  #Calculated the average proportions of the species pool with significant effects for different numbers of functions
  ref <- aggregate(posCurve2_ref$div, by=list("nfunc"=posCurve2_ref$nfunc), mean)
  names(ref)[2] <- "ref"
  rand <- aggregate(posCurve2_rand$div, by=list("nfunc"=posCurve2_rand$nfunc), mean)
  names(rand)[2] <- "rand"
  temp <- merge(ref, rand)
  
  #return the results
  return(temp)
}

##############################################################################################################
#                                                                                                            #
#   Test for effects of the different properties in the simulated data on results of the turnover approach   #
#                                                                                                            #
##############################################################################################################

#Test for the range of simulated results_____________________________________________________________________________
R1 <- addPropLine()
R2 <- addPropLine()
R3 <- addPropLine()
R4 <- addPropLine()
R5 <- addPropLine()
R6 <- addPropLine()
R7 <- addPropLine()
R8 <- addPropLine()
R9 <- addPropLine()
R10 <- addPropLine()

#Plot the results
x11(height=7,width=14)
par(mai=c(2, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
par(mfrow=c(1,2))
plot(c(1,10), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nsimulated data"))
for (i in c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")) 
lines(ref ~ nfunc, data=get(i), lwd=2)

plot(c(1,10), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\npermutated data"))
for (i in c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")) 
  lines(rand ~ nfunc, data=get(i), lwd=2)

savePlot("TenRunsWithStandarSetting", type="wmf")



#Test for the effect of number of plots on the simulated functions_____________________________________________
N20 <- addPropLine(N_plots=20)
N40 <- addPropLine(N_plots=40)
N60 <- addPropLine(N_plots=60)
N80 <- addPropLine(N_plots=80)
N100 <- addPropLine(N_plots=100)
N150 <- addPropLine(N_plots=152)
N200 <- addPropLine(N_plots=200)

#Plot the results
x11(height=7,width=14)
  par(mai=c(1, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
  m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
  layout(mat = m ,heights = c(0.85,0.15))
  
  #par(mfrow=c(1,2))
  plot(c(1,10), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nsimulated data"))
  lines(ref ~ nfunc, data=N20, lwd=1)
  lines(ref ~ nfunc, data=N40, lwd=2)
  lines(ref ~ nfunc, data=N60, lwd=3)
  lines(ref ~ nfunc, data=N80, lwd=4)
  lines(ref ~ nfunc, data=N100, lwd=5)
  lines(ref ~ nfunc, data=N150, lwd=6)
  lines(ref ~ nfunc, data=N200, lwd=7)
  
  plot(c(1,10), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\npermutated data"))
  lines(rand ~ nfunc, data=N20, lwd=1)
  lines(rand ~ nfunc, data=N40, lwd=2)
  lines(rand ~ nfunc, data=N60, lwd=3)
  lines(rand ~ nfunc, data=N80, lwd=4)
  lines(rand ~ nfunc, data=N100, lwd=5)
  lines(rand ~ nfunc, data=N150, lwd=6)
  lines(rand ~ nfunc, data=N200, lwd=7)
  
  par(mai=c(0.1, 2, 0.1, 0.1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top", inset = 0, legend=c("20", "40", "60", "80", "100", "152", "200"), lwd=c(1,2,3,4,5,6,7), horiz=TRUE, title="Number of plots")
savePlot("EffectOfInreasingNumberOfPlots", type="wmf")

#Test for the effect of number of species_________________________________________________________________________________________
S10 <- addPropLine(N_species=10, N_funct=40)
S15 <- addPropLine(N_species=15, N_funct=40)
S20 <- addPropLine(N_species=20, N_funct=40)
S40 <- addPropLine(N_species=40, N_funct=40)
S60 <- addPropLine(N_species=60, N_funct=40)

#Plot the results
x11(height=7,width=14)
  par(mai=c(1, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
  m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
  layout(mat = m ,heights = c(0.85,0.15))

  plot(c(1,40), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nreal data"))
  lines(ref ~ nfunc, data=S10, lwd=1)
  lines(ref ~ nfunc, data=S15, lwd=2)
  lines(ref ~ nfunc, data=S20, lwd=3)
  lines(ref ~ nfunc, data=S40, lwd=4)
  lines(ref ~ nfunc, data=S60, lwd=5)

  plot(c(1,40), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nrandomized data"))
  lines(rand ~ nfunc, data=S10, lwd=1)
  lines(rand ~ nfunc, data=S15, lwd=2)
  lines(rand ~ nfunc, data=S20, lwd=3)
  lines(rand ~ nfunc, data=S40, lwd=4)
  lines(rand ~ nfunc, data=S60, lwd=5)

  par(mai=c(0.1, 2, 0.1, 0.1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top", inset = 0, legend=c("10", "15", "20", "40", "60"), lwd=c(1,2,3,4,5), horiz=TRUE, title="Number of species", xpd=TRUE)
savePlot("EffectOfInreasingNumberOfSpecies", type="wmf")

#Test for the effect of number of functions_________________________________________________________________________________________
F5  <- addPropLine(N_funct=5)
F10 <- addPropLine(N_funct=10)
F15 <- addPropLine(N_funct=15)
F20 <- addPropLine(N_funct=20)
F40 <- addPropLine(N_funct=40)
F60 <- addPropLine(N_funct=60)
F80 <- addPropLine(N_funct=80)

#Plot the results
x11(height=7,width=14)
  par(mai=c(1, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
  m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
  layout(mat = m ,heights = c(0.85,0.15))
  plot(c(1,80), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nreal data"))
  lines(ref ~ nfunc, data=F5, lwd=7)
  lines(ref ~ nfunc, data=F10, lwd=6)
  lines(ref ~ nfunc, data=F15, lwd=5)
  lines(ref ~ nfunc, data=F20, lwd=4)
  lines(ref ~ nfunc, data=F40, lwd=3)
  lines(ref ~ nfunc, data=F60, lwd=2)
  lines(ref ~ nfunc, data=F80, lwd=1)

  plot(c(1,40), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\npermutated data"))
  lines(rand ~ nfunc, data=F5, lwd=7)
  lines(rand ~ nfunc, data=F10, lwd=6)
  lines(rand ~ nfunc, data=F15, lwd=5)
  lines(rand ~ nfunc, data=F20, lwd=4)
  lines(rand ~ nfunc, data=F40, lwd=3)
  lines(rand ~ nfunc, data=F60, lwd=2)
  lines(rand ~ nfunc, data=F80, lwd=1)

  par(mai=c(0.1, 2, 0.1, 0.1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top", inset = 0, legend=c("5", "10", "15", "20", "40", "60", "80"), lwd=c(7,6,5,4,3,2,1), horiz=TRUE, title="Number of functions")
savePlot(paste("EffectOfInreasingNumberOfFunctions_PropOfSpecies",i,sep="_"), type="wmf")


#Test for the proportion of species with significant effects on the simulated functions_____________________________________________
P00 <- addPropLine(PropSpecies=0.0, PropPool=1)#, N_funct = 10)
P02 <- addPropLine(PropSpecies=0.2, PropPool=1)#, N_funct = 10)
P04 <- addPropLine(PropSpecies=0.4, PropPool=1)#, N_funct = 10)
P06 <- addPropLine(PropSpecies=0.6, PropPool=1)#, N_funct = 10)
P08 <- addPropLine(PropSpecies=0.8, PropPool=1)#, N_funct = 10)
P10 <- addPropLine(PropSpecies=1.0, PropPool=1)#, N_funct = 10)

#Plot the results
x11(height=7,width=14)
  par(mai=c(1, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
  m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
  layout(mat = m ,heights = c(0.85,0.15))
  plot(c(1,20), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nsimulated data"))
  lines(ref ~ nfunc, data=P00, lwd=1)
  lines(ref ~ nfunc, data=P02, lwd=2)
  lines(ref ~ nfunc, data=P04, lwd=3)
  lines(ref ~ nfunc, data=P06, lwd=4)
  lines(ref ~ nfunc, data=P08, lwd=5)
  lines(ref ~ nfunc, data=P10, lwd=6)

  plot(c(1,20), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\npermutateded data"))
  lines(rand ~ nfunc, data=P00, lwd=1)
  lines(rand ~ nfunc, data=P02, lwd=2)
  lines(rand ~ nfunc, data=P04, lwd=3)
  lines(rand ~ nfunc, data=P06, lwd=4)
  lines(rand ~ nfunc, data=P08, lwd=5)
  lines(rand ~ nfunc, data=P10, lwd=6)

  par(mai=c(0.1, 2, 0.1, 0.1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top", inset = 0, legend=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), lwd=c(1,2,3,4,5,6), horiz=TRUE, title="Proportion of species contributing to individual functions", xpd=TRUE)
savePlot("EffectsOfIncreasingProportionOfSpeciesWithEffects", type="wmf")


#Test for effects of the size of the affecting pool_____________________________________________
PP00 <- addPropLine(PropPool=0.0)#, N_funct = 10)
PP02 <- addPropLine(PropPool=0.2)#, N_funct = 10)
PP04 <- addPropLine(PropPool=0.4)#, N_funct = 10)
PP06 <- addPropLine(PropPool=0.6)#, N_funct = 10)
PP08 <- addPropLine(PropPool=0.8)#, N_funct = 10)
PP10 <- addPropLine(PropPool=1.0)#, N_funct = 10)

#Plot the results
x11(height=7,width=14)
  par(mai=c(1, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
  m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
  layout(mat = m ,heights = c(0.85,0.15))
  plot(c(1,20), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nsimulated data"))
  lines(ref ~ nfunc, data=PP00, lwd=1)
  lines(ref ~ nfunc, data=PP02, lwd=2)
  lines(ref ~ nfunc, data=PP04, lwd=3)
  lines(ref ~ nfunc, data=PP06, lwd=4)
  lines(ref ~ nfunc, data=PP08, lwd=5)
  lines(ref ~ nfunc, data=PP10, lwd=6)

  plot(c(1,20), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\npermutateded data"))
  lines(rand ~ nfunc, data=PP00, lwd=1)
  lines(rand ~ nfunc, data=PP02, lwd=2)
  lines(rand ~ nfunc, data=PP04, lwd=3)
  lines(rand ~ nfunc, data=PP06, lwd=4)
  lines(rand ~ nfunc, data=PP08, lwd=5)
  lines(rand ~ nfunc, data=PP10, lwd=6)

  par(mai=c(0.1, 2, 0.1, 0.1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top", inset = 0, legend=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), lwd=c(1,2,3,4,5,6), horiz=TRUE, title="Proportion of species pool with effects", xpd=TRUE)
savePlot("EffectsOfIncreasingProportionOfSpeciesWithEffects", type="wmf")


#Test for the effect of changes in the signal/noise ratio____________________________________________
for (i in c(qchisq(0.16, 1, lower.tail=F), qchisq(0.01, 1, lower.tail=F), qchisq(0.001, 1, lower.tail=F), qchisq(0.0001, 1, lower.tail=F)))
{
  C1 <- addPropLine(SRatio = 1000)
  C2 <- addPropLine(SRatio = 100)
  C3 <- addPropLine(SRatio = 10)
  C4 <- addPropLine(SRatio = 1)
  C5 <- addPropLine(SRatio = 0.1)
  C6 <- addPropLine(SRatio = 0.01)
  C7 <- addPropLine(SRatio = 0.001)
  
x11(height=7,width=14)
  par(mai=c(1, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
  m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
  layout(mat = m ,heights = c(0.85,0.15))
  plot(c(1,20), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nsimulated data"))
  lines(ref ~ nfunc, data=C1, lwd=1)
  lines(ref ~ nfunc, data=C2, lwd=2)
  lines(ref ~ nfunc, data=C3, lwd=3)
  lines(ref ~ nfunc, data=C4, lwd=4)
  lines(ref ~ nfunc, data=C5, lwd=5)
  lines(ref ~ nfunc, data=C6, lwd=6)
  lines(ref ~ nfunc, data=C7, lwd=7)
  
  plot(c(1,20), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\npermutateded data"))
  lines(rand ~ nfunc, data=C1, lwd=1)
  lines(rand ~ nfunc, data=C2, lwd=2)
  lines(rand ~ nfunc, data=C3, lwd=3)
  lines(rand ~ nfunc, data=C4, lwd=4)
  lines(rand ~ nfunc, data=C5, lwd=5)
  lines(rand ~ nfunc, data=C6, lwd=6)
  lines(rand ~ nfunc, data=C7, lwd=7)
  
  par(mai=c(0.1, 2, 0.1, 0.1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top", inset = 0, legend=c("1000", "100", "10", "1", "0.1", "0.01", "0.001"), lwd=c(1,2,3,4,5,6), horiz=TRUE, title="Factor by which noise in simulation is devided", xpd=TRUE)
savePlot("EffectsOfIncreasingNoiseEffectRatio", type="wmf")
  
  
  
#Test the effect of using stricter k-values when calculating the_________________________________________________________________  
PP00_0001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0, k=qchisq(0.0001, 1, lower.tail=F))
PP02_0001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.2, k=qchisq(0.0001, 1, lower.tail=F))
PP04_0001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.4, k=qchisq(0.0001, 1, lower.tail=F))
PP06_0001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.6, k=qchisq(0.0001, 1, lower.tail=F))
PP08_0001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.8, k=qchisq(0.0001, 1, lower.tail=F))
PP10_0001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=1.0, k=qchisq(0.0001, 1, lower.tail=F))

PP00_001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0, k=qchisq(0.001, 1, lower.tail=F))
PP02_001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.2, k=qchisq(0.001, 1, lower.tail=F))
PP04_001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.4, k=qchisq(0.001, 1, lower.tail=F))
PP06_001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.6, k=qchisq(0.001, 1, lower.tail=F))
PP08_001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.8, k=qchisq(0.001, 1, lower.tail=F))
PP10_001 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=1.0, k=qchisq(0.001, 1, lower.tail=F))

PP00_01 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0, k=qchisq(0.01, 1, lower.tail=F))
PP02_01 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.2, k=qchisq(0.01, 1, lower.tail=F))
PP04_01 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.4, k=qchisq(0.01, 1, lower.tail=F))
PP06_01 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.6, k=qchisq(0.01, 1, lower.tail=F))
PP08_01 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.8, k=qchisq(0.01, 1, lower.tail=F))
PP10_01 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=1.0, k=qchisq(0.01, 1, lower.tail=F))

PP00_05 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0, k=qchisq(0.05, 1, lower.tail=F))
PP02_05 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.2, k=qchisq(0.05, 1, lower.tail=F))
PP04_05 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.4, k=qchisq(0.05, 1, lower.tail=F))
PP06_05 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.6, k=qchisq(0.05, 1, lower.tail=F))
PP08_05 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.8, k=qchisq(0.05, 1, lower.tail=F))
PP10_05 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=1.0, k=qchisq(0.05, 1, lower.tail=F))

#classic k-values of 2
PP00_16 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0, k=qchisq(0.16, 1, lower.tail=F))
PP02_16 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.2, k=qchisq(0.16, 1, lower.tail=F))
PP04_16 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.4, k=qchisq(0.16, 1, lower.tail=F))
PP06_16 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.6, k=qchisq(0.16, 1, lower.tail=F))
PP08_16 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=0.8, k=qchisq(0.16, 1, lower.tail=F))
PP10_16 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.5, N_funct=20, PropPool=1.0, k=qchisq(0.16, 1, lower.tail=F))

#Plot the results
#calles separately for the different significance levels
x11(height=7,width=14)
par(mai=c(2, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
par(mfrow=c(1,2))
plot(c(1,20), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nsimulated data"))
par(col="green")
lines(ref ~ nfunc, data=PP00_16, lwd=1)
lines(ref ~ nfunc, data=PP02_16, lwd=2)
lines(ref ~ nfunc, data=PP04_16, lwd=3)
lines(ref ~ nfunc, data=PP06_16, lwd=4)
lines(ref ~ nfunc, data=PP08_16, lwd=5)
lines(ref ~ nfunc, data=PP10_16, lwd=6)
par(col="yellow")
lines(ref ~ nfunc, data=PP00_05, lwd=1)
lines(ref ~ nfunc, data=PP02_05, lwd=2)
lines(ref ~ nfunc, data=PP04_05, lwd=3)
lines(ref ~ nfunc, data=PP06_05, lwd=4)
lines(ref ~ nfunc, data=PP08_05, lwd=5)
lines(ref ~ nfunc, data=PP10_05, lwd=6)
par(col="orange")
lines(ref ~ nfunc, data=PP00_01, lwd=1)
lines(ref ~ nfunc, data=PP02_01, lwd=2)
lines(ref ~ nfunc, data=PP04_01, lwd=3)
lines(ref ~ nfunc, data=PP06_01, lwd=4)
lines(ref ~ nfunc, data=PP08_01, lwd=5)
lines(ref ~ nfunc, data=PP10_01, lwd=6)
par(col="red")
lines(ref ~ nfunc, data=PP00_001, lwd=1)
lines(ref ~ nfunc, data=PP02_001, lwd=2)
lines(ref ~ nfunc, data=PP04_001, lwd=3)
lines(ref ~ nfunc, data=PP06_001, lwd=4)
lines(ref ~ nfunc, data=PP08_001, lwd=5)
lines(ref ~ nfunc, data=PP10_001, lwd=6)
par(col="pink")
lines(ref ~ nfunc, data=PP00_0001, lwd=1)
lines(ref ~ nfunc, data=PP02_0001, lwd=2)
lines(ref ~ nfunc, data=PP04_0001, lwd=3)
lines(ref ~ nfunc, data=PP06_0001, lwd=4)
lines(ref ~ nfunc, data=PP08_0001, lwd=5)
lines(ref ~ nfunc, data=PP10_0001, lwd=6)

par(col="black")
legend("bottom", inset=c(-1,-0.55), c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), lwd=c(1,2,3,4,5,6), horiz=TRUE, title="Proportion of species pool contributing", xpd=TRUE)

plot(c(1,20), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\npermutated simulated data"))
par(col="green")
lines(rand ~ nfunc, data=PP00_16, lwd=1)
lines(rand ~ nfunc, data=PP02_16, lwd=2)
lines(rand ~ nfunc, data=PP04_16, lwd=3)
lines(rand ~ nfunc, data=PP06_16, lwd=4)
lines(rand ~ nfunc, data=PP08_16, lwd=5)
lines(rand ~ nfunc, data=PP10_16, lwd=6)
par(col="yellow")
lines(rand ~ nfunc, data=PP00_05, lwd=1)
lines(rand ~ nfunc, data=PP02_05, lwd=2)
lines(rand ~ nfunc, data=PP04_05, lwd=3)
lines(rand ~ nfunc, data=PP06_05, lwd=4)
lines(rand ~ nfunc, data=PP08_05, lwd=5)
lines(rand ~ nfunc, data=PP10_05, lwd=6)
par(col="orange")
lines(rand ~ nfunc, data=PP00_01, lwd=1)
lines(rand ~ nfunc, data=PP02_01, lwd=2)
lines(rand ~ nfunc, data=PP04_01, lwd=3)
lines(rand ~ nfunc, data=PP06_01, lwd=4)
lines(rand ~ nfunc, data=PP08_01, lwd=5)
lines(rand ~ nfunc, data=PP10_01, lwd=6)
par(col="red")
lines(rand ~ nfunc, data=PP00_001, lwd=1)
lines(rand ~ nfunc, data=PP02_001, lwd=2)
lines(rand ~ nfunc, data=PP04_001, lwd=3)
lines(rand ~ nfunc, data=PP06_001, lwd=4)
lines(rand ~ nfunc, data=PP08_001, lwd=5)
lines(rand ~ nfunc, data=PP10_001, lwd=6)
par(col="pink")
lines(rand ~ nfunc, data=PP00_0001, lwd=1)
lines(rand ~ nfunc, data=PP02_0001, lwd=2)
lines(rand ~ nfunc, data=PP04_0001, lwd=3)
lines(rand ~ nfunc, data=PP06_0001, lwd=4)
lines(rand ~ nfunc, data=PP08_0001, lwd=5)
lines(rand ~ nfunc, data=PP10_0001, lwd=6)

par(col="black")
legend("bottom", inset=c(-1,-0.55), c("0.16", "0.05", "0.01", "0.001", "0.0001"), lwd=3, col=c("green", "yellow", "orange", "red", "pink"), horiz=TRUE, title="Significance level", xpd=TRUE)

savePlot("EffectsOfIncreasingSizeOAffectingPool_strictInclusion0001", type="wmf")


#Test for conditions like in the main analysis_________________________________________________________________________________________
Jena  <- addPropLine(N_plots=80, N_species=60, PropSpecies=0.2, N_funct=82)
x11(height=7,width=14)
par(mai=c(2, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
par(mfrow=c(1,2))
plot(c(1,82), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of  of species pool\nreal data"))
lines(ref ~ nfunc, data=Jena, lwd=1)
plot(c(1,82), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of  of species pool\nrandomized data"))
lines(rand ~ nfunc, data=Jena, lwd=1)


#Test for the effect of changes in the signal/noise ratio
#at different proportions of species contributing to functioning

for (i in c(qchisq(0.16, 1, lower.tail=F), qchisq(0.01, 1, lower.tail=F), qchisq(0.001, 1, lower.tail=F), qchisq(0.0001, 1, lower.tail=F)))
{
  C1 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.2, N_funct=40, PropPool = 0.5, k=i, SRatio = 100)
  C2 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.2, N_funct=40, PropPool = 0.5, k=i, SRatio = 10)
  C3 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.2, N_funct=40, PropPool = 0.5, k=i, SRatio = 5)
  C4 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.2, N_funct=40, PropPool = 0.5, k=i, SRatio = 3)
  C5 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.2, N_funct=40, PropPool = 0.5, k=i, SRatio = 2)
  C6 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.2, N_funct=40, PropPool = 0.5, k=i, SRatio = 1)
  C7 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.2, N_funct=40, PropPool = 0.5, k=i, SRatio = 0.5)
  C8 <- addPropLine(N_plots=40, N_species=20, PropSpecies=0.2, N_funct=40, PropPool = 0.5, k=i, SRatio = 0.2)
  
  #Plot the results
  x11(height=7,width=14)
  par(mai=c(2, 2, 0.1, 0.1), mgp=c(4,1.5,0), cex.lab=2, cex.axis=1.5, pch=21, cex=1, lwd=1.5, tck=0.015)
  par(mfrow=c(1,2))
  plot(c(1,40), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\nsimulated data"))
  lines(ref ~ nfunc, data=C1, lwd=1)
  lines(ref ~ nfunc, data=C2, lwd=2)
  lines(ref ~ nfunc, data=C3, lwd=3)
  lines(ref ~ nfunc, data=C4, lwd=4)
  lines(ref ~ nfunc, data=C5, lwd=1, lty=2)
  lines(ref ~ nfunc, data=C6, lwd=2, lty=2)
  lines(ref ~ nfunc, data=C7, lwd=3, lty=2)
  lines(ref ~ nfunc, data=C8, lwd=4, lty=2)
  legend("bottom", inset=c(1,-0.55), c("1/100", "1/10", "1/5", "1/3"), lwd=c(1,2,3,4), lty=c(1), horiz=TRUE, title="Signal to noise ratio", xpd=TRUE)
  
  
  plot(c(1,40), c(0,1), col="transparent", xlab="Number of functions", ylab=("Proportion of species\npermutated data"))
  lines(rand ~ nfunc, data=C1, lwd=1)
  lines(rand ~ nfunc, data=C2, lwd=2)
  lines(rand ~ nfunc, data=C3, lwd=3)
  lines(rand ~ nfunc, data=C4, lwd=4)
  lines(rand ~ nfunc, data=C5, lwd=1, lty=2)
  lines(rand ~ nfunc, data=C6, lwd=2, lty=2)
  lines(rand ~ nfunc, data=C7, lwd=3, lty=2)
  lines(rand ~ nfunc, data=C8, lwd=4, lty=2)
  legend("bottom", inset=c(1,-0.55), c("1/2", "1", "1/0.5", "1/0.2"), lwd=c(1,2,3,4), lty=c(2), horiz=TRUE, title="Signal to noise ratio", xpd=TRUE)
  
  savePlot(paste("EffectOfSignalNoiseRati_k",i,sep="_"), type="wmf")
}



