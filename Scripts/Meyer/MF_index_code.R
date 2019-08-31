###########################################################################################################################
#                                                                                                                         #
# R-Code to calculate an index of multifunctionality using a multivariate approach as presented in                        #
#                                                                                                                         #
# Meyer et al. (2017) Biodiversity-multifunctionality relationships depend on identity and number of measured functions   #
# Nature Ecology & Evolution                                                                                              #
#                                                                                                                         #
# The index is caclulated based on a principal componant analysis (PCA) using data on measurements of various ecosystem   #
# functions at a range of sites. To calculate the PCA the rda function from the vegan package is used. Therefore no       #
# missing values in any function are allowed.                                                                             #
#                                                                                                                         #
# For each variable the information if higher or lower values indicate higher functioning needs to be defined by hand.    #
# This information is used to give the PCA-axis biological meaning, so that for all axes higher values indicate higher    #
# functioning.                                                                                                            #
#                                                                                                                         #
# As described in the paper the index of multifunctionality is calculated by summing up all site-scores of the PCA after  #
# having weighted the axes by their eigenvalues.                                                                          # 
#                                                                                                                         #
# In the below code the data from the paper is used for illustration. That is 82 ecosystem variables that have been       #
# measured on 81 plots spanning a gradient of plant species richness from monocultures to 60-species mixtures             #
# in the Jena Experiment.                                                                                                 #
#                                                                                                                         #
###########################################################################################################################

#Version 1.0 by Sebastian T. Meyer (Sebastian.T.Meyer@TUM.de), 2017/10/07

#Files used

#PlotInformation.csv: Explanatory variables describing the plots of the Jena Experiment
#plotcode: name of the plot
#block: spatial blocks in which plots are arranged on the field site
#sowndiv: sown species richness of the plot
#numfg: number of functional groups sown in the plot
#num[...]: number of species from grasses, small herbs, tall herbs, and legumes sown in the plot
#[...].ef: effect of functional groups as presence (1) and absence (0) of grasses, small herbs, tall herbs, and legumes

#DATA_Processes.csv: the measured values for the 82 variables indicating functions in the 81 plots of the Jena Experiment

#VariableList.csv: Information for each variable whether higher or lower values are indication of higher functioning
#DirectionOfBetter: Variable indicating that higher values are considered higher functioning (1), lower values are 
#                   considered higher functioning (-1) or no clear direction can be defined (0).

#requirements
library(vegan)
library(readr)
library(here)

#Read in data from the online Suppl. 

tmp = tempfile(fileext = ".csv")
download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41559-017-0391-4/MediaObjects/41559_2017_391_MOESM3_ESM.csv", destfile = tmp, mode="wb")
PlotInfo <- read_csv2(tmp)

PlotInfo$div.col<-rainbow(6, start=0, end=0.88)[6:1][as.numeric(factor(PlotInfo$sowndiv))] # color-code for diversity levels

tmp = tempfile(fileext = ".csv")
download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41559-017-0391-4/MediaObjects/41559_2017_391_MOESM4_ESM.csv", destfile = tmp, mode="wb")
DataProcess <- read_delim(tmp, col_types = paste(c("c", rep("d", 82)), collapse = ""), delim = ";")

#Center the variables around zero and standardize to variance of one standard deviation to homogenize scales  
DataProcessStand <- DataProcess
DataProcessStand[,2:ncol(DataProcessStand)] <- decostand(DataProcessStand[,2:ncol(DataProcessStand)], meth="stand")
depStand <- DataProcessStand[,2:ncol(DataProcessStand)]
row.names(depStand) <- DataProcessStand$Plot

tmp = tempfile(fileext = ".csv")
download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41559-017-0391-4/MediaObjects/41559_2017_391_MOESM6_ESM.csv", destfile = tmp, mode="wb")
VList <- read_csv2(tmp)

#Calculate the PCA
pca1<-rda(depStand)
summary(pca1)

#Here the orientation of all axes is homogenized based on their biological meaning
SiteSc <- as.data.frame(pca1$CA$u) 
SpeciesSc <- as.data.frame(pca1$CA$v) 

for(i in 1:length(SiteSc[1,])) 
  {
  IdentV <- row.names(SpeciesSc)[which(abs(SpeciesSc[,i])==max(abs(SpeciesSc[,i])))] #identifies the variable with highest loading on an axis
  Soll <- as.numeric(unique(VList[VList$Variable==IdentV,"DirectionOfBetter"]))
  if (Soll == 0) Soll <- 1 #Variables for which no clear direction of better was defined are kept as they are
  Ist <- if (SpeciesSc[IdentV,i]<0) -1 else 1
  if (Soll != Ist) 
  {
    SiteSc[,i] <- SiteSc[,i]*-1
    SpeciesSc[,i] <- SpeciesSc[,i]*-1
  }
}
# Check of new orientation
# for (i in 1: ncol(pca1$CA$u)) {plot(pca1$CA$u[,i] ~ SiteSc[,i], main=i)
# readline()}

#create new pca-pbject with biologically orientated axis scores
pca2 <- pca1
pca2$CA$u <- as.matrix(SiteSc)
pca2$CA$v <- as.matrix(SpeciesSc)

#plot the first two PCA-axes
png(filename=here("Figures", "Fig1A_V2.png"), type="cairo", units="px", pointsize=20, width=2000, height=2000, res=200)
par(mai=c(2, 2, 0.2, 0.2), mgp=c(2,0.5,0), cex.lab=1.5, cex.axis=1.0, pch=21, cex=1.5, lwd=5, tck=0.015)
pl <- biplot(pca2, scal=3, choices=1:2, col="grey", xlab="PCA-axis 1", ylab="PCA-axis 2")
#points in gray
points(pca2, col="grey", pch=19, scal=3, choices=1:2)
#points with color according to diversity gradient
points(pca2, col=PlotInfo$div.col, pch=19, scal=3, choices=1:2)
dev.off()

#Calculate the mutivariate index of multifunctionality by summing the site scores of the PCA
#To calculate the index the axes are weighted by their eigenvalues
temp2 <- scores(pca2, choices=1:80, display=c("sites")) 
eig<-summary(pca2)$cont$importance[1,]
for(i in 1:length(eig)) temp2[,i] <- temp2[,i] * eig[i]
Index.wt <- rowSums(temp2)
Index.wt <- as.data.frame(Index.wt)
Index.wt$plotcode <- row.names(Index.wt)
Index.wt <- merge(PlotInfo, Index.wt)

#plot the multifunctionality index together with prediction of linear model for the multifunctionality-diversity relationship
png(filename=here("Figures", "Fig1B.png"), type="cairo", units="px", pointsize=20, width=2000, height=2000, res=200)
par(mai=c(2, 2, 0.2, 0.2), mgp=c(2,0.5,0), cex.lab=1.5, cex.axis=1.0, pch=21, cex=1.5, lwd=5, tck=0.015)
color <- rainbow(6, start=0, end=0.88)[6:1][as.numeric(factor(Index.wt$sowndiv))]
plot(Index.wt ~ log2(sowndiv), data=Index.wt, ylab="Index of multifunctionality", xlab="Plant species richness (log2)", pch=19, col=color)
mod.wt <- lm(Index.wt ~ block+log2(sowndiv), data=Index.wt)
new=seq(1,60, by=0.5); lines(log2(new), rowMeans(cbind(predict(mod.wt, data.frame(block="B1", sowndiv=new)), predict(mod.wt, data.frame(block="B2", sowndiv=new)), predict(mod.wt, data.frame(block="B3", sowndiv=new)), predict(mod.wt, data.frame(block="B4", sowndiv=new)))))
lines(log2(new), rowMeans(cbind(
  predict(mod.wt, data.frame(block="B1", sowndiv=new), interval ="confidence", level=0.95)[,2], 
  predict(mod.wt, data.frame(block="B2", sowndiv=new), interval ="confidence", level=0.95)[,2], 
  predict(mod.wt, data.frame(block="B3", sowndiv=new), interval ="confidence", level=0.95)[,2], 
  predict(mod.wt, data.frame(block="B4", sowndiv=new), interval ="confidence", level=0.95)[,2])), lty=2)
lines(log2(new), rowMeans(cbind(
  predict(mod.wt, data.frame(block="B1", sowndiv=new), interval ="confidence", level=0.95)[,3], 
  predict(mod.wt, data.frame(block="B2", sowndiv=new), interval ="confidence", level=0.95)[,3], 
  predict(mod.wt, data.frame(block="B3", sowndiv=new), interval ="confidence", level=0.95)[,3], 
  predict(mod.wt, data.frame(block="B4", sowndiv=new), interval ="confidence", level=0.95)[,3])), lty=2)
dev.off()

#Effect of plant species richness on multifunctionality
anova(mod.wt)

