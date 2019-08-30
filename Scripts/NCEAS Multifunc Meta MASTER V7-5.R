#######################################################################################################
#                                                                                                     #
#                 META-ANALYSIS OF BIODIVERSITY AND ECOSYSTEM MULTIFUNCTIONALITY                      #
#                                                                                                     #
#######################################################################################################

#Author: Jon Lefcheck & Jarrett Byrnes

#Last updated: V7-5.2014-10-30

#######################################################################################################
#                                        TABLE OF CONTENTS                                            #
#   Line ##: Required libraries                                                                       #
#   Line ##: Importing and formatting the data                                                        #
#   Line ##: Data exploration                                                                         #
#   Line ##: Averaging approach                                                                       #
#   Line ##: Multiple threshold approach                                                              #
#   Line ##: Turnover approach                                                                        #
#   Line ##: Multiplicative approach                                                                  #
#                                                                                                     #
#######################################################################################################

library(ggplot2) #Calls: ggplot
library(gridExtra) #Calls: grid.arrange
library(MASS) #Calls: glmmPQL
library(nlme) #Calls: lmeControl
library(plotrix) #Calls: std.error
library(plyr) #Calls: ddply, rbind.fill
library(reshape2) #Calls: melt

setwd("C:/Users/Jon/Dropbox/nceas_bdef/multifunctionality/Mutifunc meta-analysis/Analysis")
load("2014-08-20 Multifunc Meta.RData")

#######################################################################################################
#                                IMPORTING AND FORMATTING THE DATA                                    #
#######################################################################################################

#Import from file: Monoculture meta-master ALL DATA.xlxs
multifunc=read.csv("Multifunctionality Meta MASTER V7.5 (Final).csv")

#Remove large grassland studies for sensitivity analysis
#multifunc=subset(multifunc,Reference!="BIODEPTH (Hector et al, Science 1999)")
#multifunc=subset(multifunc,Reference!="Cedar Creek (E120)")
#multifunc=subset(multifunc,Reference!="BioCON (E140)")
#multifunc=subset(multifunc,Reference!="Jena (Allan et al. 2013)")

#Remove the rows where Direction!="Positive" or !="Negative"
multifunc=droplevels(subset(multifunc,multifunc$Direction=="Positive" | multifunc$Direction=="Negative"))

#Convert all response means, sample sizes, and standard deviations to numeric
#First, extract names of columns for response means, N, and SD
Y.colnames=colnames(multifunc)[ grep("Y",colnames(multifunc))[grep("Y",colnames(multifunc))>=27] ]
N.colnames=colnames(multifunc)[ grep("N",colnames(multifunc))[grep("N",colnames(multifunc))>=27] ]
SD.colnames=colnames(multifunc)[ grep("SD",colnames(multifunc))[grep("SD",colnames(multifunc))>=27] ]
#Convert response values to numeric
multifunc[,c(Y.colnames,N.colnames,SD.colnames)]=apply(multifunc[,c(Y.colnames,N.colnames,SD.colnames)],2,function(x) as.numeric(as.character(x)) )

#Check recorded number of species (Smax) against actual number of species in maximum polyculture treatment
#First, retrieve the column name of the last column with actual values
poly.colnames=apply(multifunc[,Y.colnames],1,function(x) { y=rev(x[is.finite(x)])[1]; names(y)[length(y)] })
#Use regular exprsesions to grab the number in the column name and covert them to a numeric vector
Smax.colnames=as.numeric(gsub("X([0-9]+).*","\\1",poly.colnames))
cbind(as.character(multifunc$Reference),Smax.colnames,multifunc$Smax,Smax.colnames==multifunc$Smax)

#Check to see if any experiments report only one function
byexpt=ddply(multifunc,c("Reference","Study","Expt"),nrow)
byexpt[byexpt$V1==1,c("Reference","Study","Expt")]
   
#If Direction=="Negative", then transform based on Byrnes et al. 2014 MEE: x = -x + max(x)
multifunc[,Y.colnames]=ddply(multifunc,1,function(x) 
  if(x$Direction=="Negative") -x[Y.colnames]+max(x[Y.colnames],na.rm=T) else x[Y.colnames] )[,-1]

#If any responses are negative, scale so that they are all >0
multifunc[,Y.colnames]=ddply(multifunc,1,function(x) 
  if(any(x[Y.colnames]<0,na.rm=T)) x[Y.colnames]+max(abs(x[Y.colnames]),na.rm=T) else x[Y.colnames] )[,-1] 

#Create a dataset with values scaled by the maximum value (for averaging approach)
multifunc.scaled=multifunc
multifunc.scaled[,SD.colnames]=ddply(multifunc,1,function(x) x[SD.colnames]/max(abs(x[Y.colnames]),na.rm=T) )[,-1]
multifunc.scaled[,Y.colnames]=ddply(multifunc,1,function(x) x[Y.colnames]/max(abs(x[Y.colnames]),na.rm=T) )[,-1]
  
#######################################################################################################
#                                      DATA EXPLORATION                                               #
#######################################################################################################

#Number of studies
nrow(ddply(multifunc,"Study",nrow))
#Number of experiments
nrow(ddply(multifunc,c("Study","Expt"),nrow))
#Total number of functions
nrow(multifunc)
#Number of experiments for each function
table(ddply(multifunc,c("Study","Expt"),nrow)$V1)

#Number of habitats
count(ddply(multifunc,c("Study","Expt","Sys1"),nrow),vars="Sys1")
#Number of trophic levels
count(ddply(multifunc,c("Study","Expt","FTG"),nrow),vars="FTG")
#And both
count(ddply(multifunc,c("Study","Expt","Sys1","FTG"),nrow),vars=c("FTG","Sys1"))

#Level of richness within an experiment
hist(multifunc$Smax)
median(multifunc$Smax)
range(multifunc$Smax)
#Level of richness within an experiment by habitat
ddply(multifunc,c("Sys1"),summarize,median=median(Smax))
ddply(multifunc,c("Sys1"),function(x) data.frame(min=range(x$Smax)[1],max=range(x$Smax)[2]))
#Level of richness within an experiment by trophic level
ddply(multifunc,c("FTG"),summarize,median=median(Smax))
ddply(multifunc,c("FTG"),function(x) data.frame(min=range(x$Smax)[1],max=range(x$Smax)[2]))

#Number of functions per experiment
median(ddply(multifunc,c("Study","Expt"),nrow)$V1)
range(ddply(multifunc,c("Study","Expt"),nrow)$V1)
#Number of functions per experiment by habitat
ddply(ddply(multifunc,c("Study","Expt","Sys1"),nrow),"Sys1",function(x) data.frame(median=median(x$V1)))
ddply(ddply(multifunc,c("Study","Expt","Sys1"),nrow),"Sys1",function(x) data.frame(min=range(x$V1)[1],max=range(x$V1)[2]))
#Number of functions per experiment by trophic level
ddply(ddply(multifunc,c("Study","Expt","FTG"),nrow),"FTG",function(x) data.frame(median=median(x$V1)))
ddply(ddply(multifunc,c("Study","Expt","FTG"),nrow),"FTG",function(x) data.frame(min=range(x$V1)[1],max=range(x$V1)[2]))

#Look at average pairwise correlation between all functions within a study
pairwisecor.df=ddply(multifunc,c("Reference","Study","Expt","Sys1","FTG"),function(x) {
  Smax=unique(x$Smax)
  x=x[,Y.colnames]
  cormat=cor(t(x),use="complete.obs",method=c("kendall"))
  data.frame(
    Smax=Smax,
    no.fn=nrow(x),
    avg.cor=mean(cormat[lower.tri(cormat)]) ) } ); pairwisecor.df
#And across all studies
mean(pairwisecor.df$avg.cor); std.error(pairwisecor.df$avg.cor)
#By habitat
ddply(pairwisecor.df,"Sys1",summarize,mean(avg.cor),std.error(avg.cor))
#By trophic level
ddply(pairwisecor.df,"FTG",summarize,mean(avg.cor),std.error(avg.cor))
#Plot as a function of number of functions
ggplot(pairwisecor.df,aes(x=no.fn,y=avg.cor))+
  geom_hline(yintercept=0,lwd=0.8,lty=1,col="grey30")+
  geom_point(size=3)+
  scale_x_continuous(breaks=seq(2,12,2))+
  labs(x="Number of functions",y="Average pairwise correlation")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#Look at number of positive / negative relationships 
#First cast data.frame longways
multifunc.long=melt(cbind(multifunc[,c(2:4,6:7,9,12)],multifunc[,Y.colnames]),id.vars=c(1:7),measure.vars=c(8:112))
multifunc.long$richness=suppressWarnings(
  ifelse(grepl("mono",multifunc.long$variable),1,as.numeric(gsub("X([0-9]+).*","\\1",multifunc.long$variable))) )

#Next calculate correlation between richness and functioning for each function
diversitycor.df=ddply(multifunc.long,c("Reference","Study","Expt","Sys1","FTG","Ydesc"),function(x)
    data.frame(cor=cor(x$richness,x$value,use="complete.obs",method="kendall"),
               p.value=cor.test(x$richness,x$value,na.action=na.omit,method="kendall")$p.value) )

#Calculate number of positive/negative functions for each study
propcor.df=ddply(diversitycor.df,c("Study","Expt","Reference"),function(x) 
  cbind(no.fn=length(unique(x$Ydesc)),
        sig.neg=sum(x[x$p.value<=0.05,"cor"]<0),
        sig.pos=sum(x[x$p.value<=0.05,"cor"]>0),
        neutral=length(x[x$p.value>0.05,"cor"])) )

propcor.df=melt(propcor.df,id.vars=c(1:4),measure.vars=c(5:7))
propcor.df$variable=factor(propcor.df$variable,levels=c("sig.pos","sig.neg","neutral"))
levels(propcor.df$variable)=c("Positive","Negative","Neutral")
#Set symbols based on references in simulation, below
# propcor.df=adply(propcor.df,1,function(x) {
# if(x$Reference=="Wardle et al. 2003") "diamond" else 
#   if(x$Reference=="Cedar Creek (E120)") "triangle" else "none" } )

#Plot results
ggplot(propcor.df,aes(x=no.fn,y=value/no.fn,group=variable,col=variable,shape=variable))+#,shape=V1))+
  geom_point(size=4,alpha=0.5,position="jitter")+
  scale_x_continuous(breaks=c(2,4,6,8,10,12))+
  scale_color_manual(values=c("red","blue","grey20"),name="")+
  scale_shape_manual(values=15:17,name="")+
#   scale_shape_manual(values=c(15,1,17),guide="none")+
  stat_smooth(method="glm",family=binomial(),aes(lty=variable),lwd=2,se=F)+
  scale_linetype(guide="none")+
  labs(x="Total number of functions",y="Proportion of total functions")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.direction="horizontal",legend.position="bottom")

#######################################################################################################
#                                       AVERAGING APPROACH                                            #
#######################################################################################################

#Calculate average level of functioning across all functions for each treatment, for each experiment
multifunc.avg=ddply(multifunc.scaled,c("Reference","Study","Expt","FTG","Sys1","Sys2"),function(x) {
  z=data.frame(
    richness=colnames(x[,Y.colnames]),
    no.fn=length(unique(x$Ydesc)),
    avg.fn=colMeans(x[,Y.colnames]),
    avg.fn.SD=sqrt(
      colSums((x[,N.colnames]-1)*(x[,SD.colnames]^2),na.rm=T)/
      (colSums(x[,N.colnames],na.rm=T)-nrow(x[,N.colnames])) ),
    avg.fn.N=colSums(x[,N.colnames]) ) 
  #Remove rows where there is no response (i.e., avg.fn==NA)
  z=z[!is.na(z$avg.fn),]
  #Set richness by splitting column names
  z$richness=suppressWarnings(ifelse(grepl("mono",z$richness),1,as.numeric(gsub("X([0-9]+).*","\\1",z$richness))))
  return(z) } )

#Investigate proper functional form to use
#Group data for random effects
multifunc.avg.grouped=groupedData(avg.fn~richness|Study,data=multifunc.avg)
#Fit different functional forms using non-linear mixed models
Null=nlme(avg.fn~a,fixed=a~1,random=~a~1,start=c(a=0.2),data=multifunc.avg.grouped)
Linear=nlme(avg.fn~a+b*richness,fixed=a+b~1,random=~a+b~1,start=c(a=1.5,b=1),data=multifunc.avg.grouped)
Logarithmic=nlme(avg.fn~a+b*log(richness),fixed=a+b~1,random=~a+b~1,start=c(a=1,b=1),data=multifunc.avg.grouped)
Power=nlme(avg.fn~a*richness^b,fixed=a+b~1,random=~a+b~1,start=c(a=0.2,b=2),data=multifunc.avg.grouped)
Saturating=nlme(avg.fn~richness/(k+richness),fixed=k~1,random=k~1,start=c(k=1),data=multifunc.avg.grouped)
#Compare models using AIC
AIC(Null,Linear,Logarithmic,Power,Saturating)

#Fit log relationship using linear mixed effects model, allowing slopes and intercepts to vary by Study
avgmods.list=lapply(c("unweighted","variance","sample.size"),function(i) {
  #Subset dataset to include only non-NA data points for each type of analysis
  if(i=="variance") { multifunc.avg=multifunc.avg[!is.na(multifunc.avg$avg.fn.SD),] 
  } else if(i=="sample.size") { multifunc.avg=multifunc.avg[!is.na(multifunc.avg$avg.fn.N),]
  } else { multifunc.avg }
  #Fit linear mixed effects model for each weighting scheme
  if(i=="unweighted") { 
    mod=glmmPQL(avg.fn~log(richness),random=~richness|Study,family=quasibinomial(link="identity"),
                start=c(0.5,0),data=multifunc.avg,verbose=F)
  } else if(i=="variance") { 
    mod=glmmPQL(avg.fn~log(richness),random=~richness|Study,weights=1/((multifunc.avg$avg.fn.SD^2)+0.01),
                family=quasibinomial(link="identity"),start=c(0.5,0),data=multifunc.avg,verbose=F)
  } else { 
    mod=glmmPQL(avg.fn~log(richness),random=~richness|Study,weights=sqrt(multifunc.avg$avg.fn.N),
                family=quasibinomial(link="identity"),start=c(0.5,0),data=multifunc.avg,verbose=F) }
  #Return model
  return(mod)
} )
#Append reduced model (S <= 16)
avgmods.list=append(avgmods.list,list(update(avgmods.list[[1]],data=subset(multifunc.avg,richness<=16))))
#Look at output and diagnostic plots
lapply(avgmods.list,summary); lapply(avgmods.list,plot)

#Extract predicted fits for plot
pred.df.list=lapply(seq_along(avgmods.list),function(i) {
  if(i<4) multifunc.avg=multifunc.avg else multifunc.avg=subset(multifunc.avg,!paste(Study,Expt) %in%    
    unique(paste(subset(multifunc.avg,richness>16)$Study,subset(multifunc.avg,richness>16)$Expt)))
  #Modified from: http://glmm.wikidot.com/faq 
  #Create dataframe for predicted values for overall fit
  newdata=expand.grid(richness=1:max(multifunc.avg$richness),no.fn=2:max(multifunc.avg$no.fn),avg.fn=0)
  #Generate predicted values for overall trend
  newdata$avg.fn=predict(avgmods.list[[i]],newdata,type="response",level=0)
  #Obtain model matrix
  mm=model.matrix(terms(avgmods.list[[i]]),newdata)
  #Obtain estimate of SE based on fixed effects only
  newdata$fixedSE=sqrt(diag(mm %*% tcrossprod(vcov(avgmods.list[[i]]),mm)))
  
  #Create dataframe for predicted values for each study
  newdata2=ldply(unique(multifunc.avg$Study),function(j)
    data.frame(Study=j,
               richness=1:max(subset(multifunc.avg,Study==j)$richness),
               no.fn=max(subset(multifunc.avg,Study==j)$no.fn)) )
  #Generate predicted values for each study
  newdata2$avg.fn=predict(avgmods.list[[i]],newdata2,type="response",level=1)
  
  #Return dataframes in a list
  list(newdata,newdata2)  
} )     

#Plot predicted values and confidence bands across all studies on top of fitted values for each individual study
avgplots.list=lapply(1:3,function(i) {
  ggplot()+
    #Plot raw points
    geom_point(data=multifunc.avg,aes(x=richness,y=avg.fn),size=2,col="grey60",alpha=0.5)+
    #Plot curves for each study
    geom_line(data=pred.df.list[[i]][[2]],aes(x=richness,y=avg.fn,group=Study),col="black",lwd=0.8,alpha=0.7)+
    #Add confidence band for mixed model predictions based on fixed effects only
    geom_ribbon(data=pred.df.list[[i]][[1]],
                aes(x=richness,y=avg.fn,ymax=avg.fn+2*fixedSE,ymin=avg.fn-2*fixedSE),
                fill="red",alpha=0.4)+
    #Add line for mixed model predictions
    geom_line(data=pred.df.list[[i]][[1]],aes(x=richness,y=avg.fn),col="red",lwd=1.5,alpha=0.9)+
    scale_x_continuous(breaks=c(1,20,40,60))+
    scale_y_continuous(limits=c(0,1.1),breaks=c(0,0.5,1))+
    labs(x="Richness",y="Average multifunctionality")+
    theme_bw(base_size=18)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
} )#; avgplots.list

#Plot figure 1 with inset from list above (5 x 4.5")
p1=ggplot(data=subset(multifunc.avg,!paste(Study,Expt) %in%    
                     unique(paste(subset(multifunc.avg,richness>16)$Study,subset(multifunc.avg,richness>16)$Expt))),
       aes(x=richness,y=avg.fn))+
  #Plot raw points
  geom_point(size=2.5,col="grey60",alpha=0.3)+
  #Plot curves for each study
  geom_line(data=pred.df.list[[4]][[2]],aes(x=richness,y=avg.fn,group=Study),col="grey20",lwd=0.75,alpha=0.8)+
  #Add confidence band for mixed model predictions based on fixed effects only
  geom_ribbon(data=pred.df.list[[4]][[1]],
              aes(x=richness,y=avg.fn,ymax=avg.fn+2*fixedSE,ymin=avg.fn-2*fixedSE),
              fill="red",alpha=0.4)+
  #Add line for mixed model predictions
  geom_line(data=pred.df.list[[4]][[1]],aes(x=richness,y=avg.fn),col="red",lwd=2.5)+
  coord_cartesian(ylim=c(-0.05,1.1))+
  scale_x_continuous(breaks=c(1,4,8,12,16),labels=c("1","4","8","12","16"))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  labs(x="Richness",y="Average multifunctionality")+ 
  geom_text(data=data.frame(
    labels=letters[1]),
    aes(x=-Inf,y=Inf,label=labels),vjust=1.5,hjust=-1.5,col="black",fontface="bold",size=9)+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  annotation_custom(grob=ggplotGrob(avgplots.list[[1]]+labs(x="",y="")+theme(plot.margin=unit(c(0,0,0,0),"cm"))),
                    xmin=7,xmax=16,ymin=-0.075,ymax=0.385)

#######################################################################################################

#Look at diversity effects on single functions
#Rearrange scaled responses so they are in the same format as the averaged dataset
multifunc.single=data.frame(
  melt(multifunc.scaled,
       id.vars=c("Reference","Study","Expt","FTG","Sys1","Sys2","Ydesc"),
       measure.vars=c(Y.colnames)),
  fn.N=melt(multifunc.scaled,
            id.vars=c("Reference","Study","Expt","FTG","Sys1","Sys2","Ydesc"),
            measure.vars=c(N.colnames))[9],
  fn.SD=melt(multifunc.scaled,
             id.vars=c("Reference","Study","Expt","FTG","Sys1","Sys2","Ydesc"),
             measure.vars=c(SD.colnames))[9] )
names(multifunc.single)[9:11]=c("fn.mean","fn.N","fn.SD")
multifunc.single=multifunc.single[!is.na(multifunc.single$fn.mean),]
multifunc.single$richness=
  suppressWarnings(ifelse(grepl("mono",multifunc.single$variable),1,as.numeric(gsub("X([0-9]+).*","\\1",multifunc.single$variable))))

#Fit linear mixed effects model allowing slopes and intercepts to vary by Study
singlemods.list=lapply(c("unweighted","variance"),function(i) { #,"sample.size"),function(i) {
    #Subset dataset to include only non-NA data points for each type of analysis
    if(i=="variance") { multifunc.single=multifunc.single[!is.na(multifunc.single$fn.SD),] 
    } else if(i=="sample.size") { multifunc.single=multifunc.single[!is.na(multifunc.single$fn.N),]
    } else { multifunc.single }
    #Fit linear mixed effects model for each weighting scheme
    if(i=="unweighted") { 
      mod=glmmPQL(fn.mean~log(richness),random=~richness|Study/Ydesc,family=quasibinomial(link="identity"),
                  start=c(0.5,0),control=lmeControl(msTol=1e-5,opt="optim"),data=multifunc.single,verbose=F)
    } else if(i=="variance") { 
      mod=glmmPQL(fn.mean~log(richness),random=~richness|Study/Ydesc,weights=1/((multifunc.single$fn.SD^2)+0.01),
                  family=quasibinomial(link="identity"),start=c(0.5,0),control=lmeControl(msTol=1e-5,opt="optim"),data=multifunc.single,verbose=F)
    } else { 
      mod=glmmPQL(fn.mean~log(richness),random=~richness|Study/Ydesc,weights=sqrt(multifunc.single$fn.N),
                  family=quasibinomial(link="identity"),start=c(0.2,0),control=lmeControl(msTol=1e-5,opt="optim"),data=multifunc.single,verbose=F) 
    }
  #Return model
  return(mod)
} )
#Look at output and diagnostic plots
lapply(singlemods.list,summary); lapply(singlemods.list,plot)

#Extract predicted fits for plot
singlepred.df.list=lapply(seq_along(singlemods.list),function(i) {
  #Modified from: http://glmm.wikidot.com/faq 
  #Create dataframe for predicted values for overall fit
  newdata=expand.grid(richness=1:max(multifunc.single$richness),
                      fn.mean=0)
  #Generate predicted values for overall trend
  newdata$fn.mean=predict(singlemods.list[[i]],newdata,type="response",level=0)
  #Obtain model matrix
  mm=model.matrix(terms(singlemods.list[[i]]),newdata)
  #Obtain estimate of SE based on fixed effects only
  newdata$fixedSE=sqrt(diag(mm %*% tcrossprod(vcov(singlemods.list[[i]]),mm)))
  
  #Create dataframe for predicted values for each study
  newdata2=ldply(unique(multifunc.single$Study),function(j) {
    ldply(unique(subset(multifunc.single,Study==j)$Expt),function(k) {
    expand.grid(Study=j,Expt=k,
                Ydesc=unique(subset(multifunc.single,Study==j)$Ydesc),
                richness=1:max(subset(multifunc.single,Study==j)$richness)) } ) } )
  #Generate predicted values for each study
  newdata2$fn.mean=predict(singlemods.list[[i]],newdata2,type="response",level=2)
  newdata2$Ydesc=paste(newdata2$Ydesc,newdata2$Study,newdata2$Expt,sep=".")
  
  #Return dataframes in a list
  list(newdata,newdata2)  
} )     

#Plot predicted values and confidence bands across all studies on top of fitted values for each individual study
singleplots.list=lapply(1:2,function(i) {
  ggplot()+
    #Plot raw points
    #geom_point(data=multifunc.single,aes(x=richness,y=fn.mean),size=2,col="grey60",alpha=0.5)+
    #Plot curves for each function
#     geom_line(data=singlepred.df.list[[i]][[2]],aes(x=richness,y=fn.mean,group=Ydesc),
#               alpha=0.3,lwd=0.8)+
    #Add confidence band for mixed model predictions based on fixed effects only
    geom_ribbon(data=singlepred.df.list[[i]][[1]],
                aes(x=richness,y=fn.mean,ymax=fn.mean+2*fixedSE,ymin=fn.mean-2*fixedSE),
                fill="red",alpha=0.3)+
    #Add line for mixed model predictions
    geom_line(data=singlepred.df.list[[i]][[1]],aes(x=richness,y=fn.mean),col="red",lwd=1.5,alpha=0.9)+
    scale_y_continuous(limits=c(0,1.1),breaks=seq(0,1,0.2))+
    labs(x="Richness",y="Functioning (single functions pooled)")+
    theme_bw(base_size=18)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
} )

#Plot with average multifunctionality
p2=ggplot()+
  geom_point(data=multifunc.single,aes(x=richness,y=fn.mean),size=2,col="grey60",alpha=0.5)+
#   geom_ribbon(data=pred.df.list[[1]][[1]],
#               aes(x=richness,y=avg.fn,ymax=avg.fn+2*fixedSE,ymin=avg.fn-2*fixedSE),
#               fill="red",alpha=0.4)+
  geom_ribbon(data=singlepred.df.list[[1]][[1]],
              aes(x=richness,y=fn.mean,ymax=fn.mean+2*fixedSE,ymin=fn.mean-2*fixedSE),
              fill="blue",alpha=0.3)+
  #Add line for mixed model predictions
#   geom_line(data=pred.df.list[[1]][[1]],aes(x=richness,y=avg.fn),col="red",lwd=1.5,alpha=0.9)+
  geom_line(data=singlepred.df.list[[1]][[1]],aes(x=richness,y=fn.mean),col="blue",lwd=1.5,alpha=0.9)+
#   scale_fill_manual(values=c("red","blue"),name="")+
  coord_cartesian(ylim=c(-0.05,1.1),xlim=c(-2,62))+
  scale_x_continuous(breaks=c(1,20,40,60))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  geom_text(data=data.frame(
    labels=letters[2]),
    aes(x=-Inf,y=Inf,label=labels),vjust=1.5,hjust=-1.5,col="black",fontface="bold",size=9)+
  labs(x="Richness",y="Scaled single functions")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position=c(0.9,0.3))

#13" x 6"
grid.arrange(p1,p2,nrow=1)

#######################################################################################################

#Fit global mod using all predictors to see whether there are significant differences in the effect of richness
#on multifunctionality across systems or trophic groups
global.mod=glmmPQL(avg.fn~log(richness)*no.fn+log(richness)*FTG+log(richness)*Sys1,random=~richness|Study,family=quasibinomial(link="identity"),
                   start=c(0.5,rep(0,13)),data=multifunc.avg,verbose=F)
summary(global.mod); plot(global.mod)

#Extract predicted fits for plot
predglobal.df=ldply(unique(multifunc.avg$Sys1),function(i) {
  ldply(unique(multifunc.avg$FTG),function(j) {
    if(nrow(subset(multifunc.avg,Sys1==i & FTG==j))==0) {
      data.frame()
    } else {    
    #Generate dataframe for predicted values
    expand.grid(
      richness=1:max(subset(multifunc.avg,Sys1==i & FTG==j)$richness),
      no.fn=1:max(subset(multifunc.avg,Sys1==i & FTG==j)$no.fn),
      Sys1=i,
      FTG=j,
      avg.fn=NA)
} } ) } )
#Generate predict values for overall trends
predglobal.df$avg.fn=predict(global.mod,predglobal.df,type="response",level=0)

#Plot the results
ggplot()+
  #Plot raw points
  geom_point(data=multifunc.avg,aes(x=richness,y=avg.fn),size=2,col="grey60",alpha=0.5)+
  #Add confidence band for mixed model predictions based on fixed effects only
#   geom_ribbon(data=pred.df.list[[i]][[1]],
#               aes(x=richness,y=avg.fn,fill=no.fn,group=no.fn,ymax=avg.fn+2*fixedSE,ymin=avg.fn-2*fixedSE),
#               alpha=0.15)+
  #Add line for mixed model predictions
  geom_line(data=predglobal.df,aes(x=richness,y=avg.fn,col=no.fn,group=no.fn),lwd=1.5,alpha=0.8)+
  scale_color_gradient(high="red",low="blue",name="Number of\nfunctions")+
  #Add facets by Sys1 and FTG
  facet_grid(Sys1~FTG,scales="free")+#,ncol=2,nrow=5)+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  labs(x="Richness",y="Average multifunctionality")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),strip.text.x=element_text(size=12),
        legend.position=c(0.075,0.25))

#######################################################################################################

#Repeat, but include number of functions as an interaction
#Fit linear mixed effects model allowing slopes and intercepts to vary by Study
avgintmods.list=lapply(c("unweighted","variance","sample.size"),function(i) {
  #Subset dataset to include only non-NA data points for each type of analysis
  if(i=="variance") { multifunc.avg=multifunc.avg[!is.na(multifunc.avg$avg.fn.SD),] 
  } else if(i=="sample.size") { multifunc.avg=multifunc.avg[!is.na(multifunc.avg$avg.fn.N),]
  } else { multifunc.avg }
  #Fit linear mixed effects model for each weighting scheme
  if(i=="unweighted") { 
    mod=glmmPQL(avg.fn~log(richness)*no.fn,random=~richness|Study,family=quasibinomial(link="identity"),
                start=c(0.5,0,0,0),data=multifunc.avg,verbose=F)
  } else if(i=="variance") { 
    mod=glmmPQL(avg.fn~log(richness)*no.fn,random=~richness|Study,weights=1/((multifunc.avg$avg.fn.SD^2)+0.01),
                family=quasibinomial(link="identity"),start=c(0.5,0,0,0),data=multifunc.avg,verbose=F)
  } else { 
    mod=glmmPQL(avg.fn~log(richness)*no.fn,random=~richness|Study,weights=sqrt(multifunc.avg$avg.fn.N),
                family=quasibinomial(link="identity"),start=c(0.5,0,0,0),data=multifunc.avg,verbose=F) }
  #Return model
  return(mod)
} )
#Look at output and diagnostic plots
lapply(avgintmods.list,summary); lapply(avgintmods.list,plot)

#Extract predicted fits for plot
predint.df.list=lapply(seq_along(avgintmods.list),function(i) {
  #Modified from: http://glmm.wikidot.com/faq 
  #Create dataframe for predicted values for overall fit
  newdata=expand.grid(richness=1:max(multifunc.avg$richness),
                      no.fn=2:max(multifunc.avg$no.fn),
                      avg.fn=0)
  #Generate predicted values for overall trend
  newdata$avg.fn=predict(avgintmods.list[[i]],newdata,type="response",level=0)
  #Obtain model matrix
  mm=model.matrix(terms(avgintmods.list[[i]]),newdata)
  #Obtain estimate of SE based on fixed effects only
  newdata$fixedSE=sqrt(diag(mm %*% tcrossprod(vcov(avgintmods.list[[i]]),mm)))
  #Return dataframe
  return(newdata)
} )     

#Plot predicted values and confidence bands across all studies on top of fitted values for each individual study
lapply(1:3,function(i) {
  ggplot()+
    #Plot raw points
    geom_point(data=multifunc.avg,aes(x=richness,y=avg.fn),size=2,col="grey60",alpha=0.5)+
    #Plot curves for each study
    #geom_line(data=pred.df.list[[i]][[2]],aes(x=richness,y=avg.fn,group=Study),alpha=0.75,lwd=1,alpha=0.75)+
    #Add confidence band for mixed model predictions based on fixed effects only
    geom_ribbon(data=predint.df.list[[i]],
                aes(x=richness,y=avg.fn,ymax=avg.fn+2*fixedSE,ymin=avg.fn-2*fixedSE,group=no.fn,fill=no.fn),
                alpha=0.15)+
    #Add line for mixed model predictions
    geom_line(data=predint.df.list[[i]],aes(x=richness,y=avg.fn,group=no.fn,col=no.fn),lwd=1.5,alpha=0.9)+
    scale_color_gradient(high="red",low="blue",name="Number of\nfunctions")+
    scale_fill_gradient(high="red",low="blue",name="Number of\nfunctions")+
    scale_y_continuous(breaks=seq(0,1,0.2))+
    labs(x="\nRichness",y="Average multifunctionality\n")+
    theme_bw(base_size=18)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
} )

#Extract effect sizes and their standard errors for each level of number of functions
avginteffect.df.list=lapply(avgintmods.list,function(i) {
  ldply(2:max(multifunc.avg$no.fn),function(j) 
    data.frame(
      no.fn=j,
      Estimate=summary(i)$tTable["log(richness)","Value"]+j*summary(i)$tTable["log(richness):no.fn","Value"],
      Std.Error=sqrt( (summary(i)$tTable["log(richness)","Std.Error"]^2+summary(i)$tTable["log(richness):no.fn","Std.Error"]^2+
                         2*i$varFix["log(richness)","log(richness):no.fn"]) ) )
) } )
      
#Plot number of functions against effect size
lapply(avginteffect.df.list,function(i) {
  ggplot(i,aes(x=no.fn,y=Estimate))+
    #Add points for effect sizes
    geom_point(size=4)+
    #And error bars for the standard errors
    geom_errorbar(aes(ymax=Estimate+2*Std.Error,ymin=Estimate-2*Std.Error),width=0)+
    #Specify axis breaks
    scale_x_continuous(breaks=c(2,4,6,8,10,12))+
    labs(x="\nNumber of functions",y="Diversity effect (+/- 95% CI)\n")+
    theme_bw(base_size=18)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
} )

#######################################################################################################
#                                 MULTIPLE THRESHOLD APPROACH                                         #
#######################################################################################################

#Create vector of thresholds from 1-99%
thresholds=(1:99)/100

#Loop through thresholds and studies to calculate number of treatments > threshold
thresholds.df=ldply(thresholds,.progress="text",function(thresh) {
  ddply(multifunc,c("Reference","Study","Expt"),function(x) {
    #Melt response values and relevant metadata
    responses=melt(cbind(x[,c(2:4,6:7,9,12)],x[,Y.colnames]),id.vars=c(1:7),measure.vars=c(8:112))
    #Remove NA values
    #responses=responses[!is.na(responses$value),]
    #Create a column for richness
    responses$richness=suppressWarnings(
      ifelse(grepl("mono",responses$variable),1,as.numeric(gsub("X([0-9]+).*","\\1",responses$variable))) )
    #Determine whether each response for each functions >= some percentage (threshold) of maximum
    responses=ddply(responses,"Ydesc",function(y) cbind(y,greater.than=y$value>=thresh*max(y$value,na.rm=T)) )
    #Summarize for each treatment
    ddply(responses,c("Study","Expt","Reference","Sys1","FTG","variable","richness"),function(z) {
      data.frame(
        threshold=thresh,
        no.fn=length(z$greater.than),
        no.fn.greater=sum(z$greater.than),
        prop.fn.greater=sum(z$greater.than)/length(z$greater.than) ) } ) } ) } )

#######################################################################################################

#Use mixed models to look at trends generally across all studies by fitting raw counts to 
#quasipoisson distribution for each level of threshold
rawmods.list=dlply(thresholds.df,"threshold",.progress="text",function(i) {
  #Set lmeControl for certain thresholds
  if(i$threshold %in% c(0.02,0.99)) control=lmeControl(opt="optim",msTol=1e-6) else 
    control=lmeControl(opt="optim")
  #Function to run models
  f=function(x) glmmPQL(no.fn.greater~richness*no.fn,random=~richness|Study,
                          family=quasipoisson(link="identity"),start=c(0.15,0.05,0.1,0.001),            
                          control=control,
                          verbose=F,data=x)
  #Run model, return NA if model fails
  safef=failwith(NA,f)
  safef(i)
} ) 

#Extract predicted fit for muscle plot
#Modified from: http://glmm.wikidot.com/faq
predictraw.df=ldply(1:99,.progress="text",function(i) {
  ldply(2:12,function(j) {
    #Create dataframe for predicted values only for richness
    newdata=expand.grid(
      threshold=i/100,
      richness=1:max(subset(thresholds.df,no.fn==j)[!is.na(subset(thresholds.df,no.fn==j)$no.fn.greater),"richness"]),
      no.fn=j,
      no.fn.greater=NA)
    if(!"glmmPQL" %in% class(rawmods.list[[i]]))  {
      return(newdata)
      } else {
    #Generate predicted values from model
    newdata$no.fn.greater=predict(rawmods.list[[i]],newdata,type="response",level=0)
    #Return newdata
    return(newdata) } 
 } ) } )

#Remove predictions that exceed the number of functions measured or drop below zero
predictraw.sub.df=ddply(predictraw.df,"no.fn",function(x) {
  x=subset(x,no.fn.greater<=x$no.fn & no.fn.greater>=0)
  data.frame(
    threshold=x$threshold,
    richness=x$richness,
    no.fn=x$no.fn,
    no.fn.greater=x$no.fn.greater) 
} )
predictraw.sub.df=ddply(predictraw.df,"no.fn",function(x) subset(x,no.fn.greater<=no.fn))

#Plot predicted values against richness
ggplot(predictraw.sub.df,aes(x=richness,y=no.fn.greater,color=threshold,group=threshold))+
  geom_line(lwd=1,alpha=0.6)+
  geom_hline(aes(yintercept=no.fn),lwd=1,alpha=0.6,lty=1)+
  scale_color_gradientn(colours=rev(rainbow(5)),name="Threshold")+
  labs(x="Richness",y="Number of functions > threshold")+
  facet_wrap(~no.fn,scales="free",nrow=4)+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position=c(0.8,0.1) )

#Break out no.fn==max(no.fn) for panel b
p1=ggplot(subset(predictraw.df,no.fn==12 & no.fn.greater<=12 & no.fn.greater>=0),aes(x=richness,y=no.fn.greater,color=threshold,group=threshold))+
  #geom_hline(yintercept=1,lwd=0.8,lty=1,col="grey30")+
  geom_line(lwd=1,alpha=0.6)+
  geom_hline(yintercept=12,lwd=1,alpha=0.6,lty=1)+
  scale_color_gradientn(colours=rev(rainbow(5)),name="Threshold")+
  labs(x="Richness",y="Number of functions > threshold")+
  scale_x_continuous(breaks=c(1,20,40,60))+
  scale_y_continuous(limits=c(-0.1,14.5),breaks=c(0,4,8,12))+
  annotate("text",x=1,y=Inf,label="a",vjust=1.5,col="black",fontface="bold",size=7)+
  theme_bw(base_size=18)+
  guides(colour=guide_colourbar(title.position="top"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.background=element_blank(),
        legend.position=c(0.7,0.9),legend.direction="horizontal",legend.box="horizontal"); p1

#Extract coefficients and standard errors
rawcoefs.df=ldply(2:max(thresholds.df$no.fn),.progress="text",function(i) {
  ldply(rawmods.list,function(j) {
    if("glmmPQL" %in% class(j)) {
    data.frame(
      no.fn=i,
      Intercept=summary(j)$tTable["(Intercept)","Value"]+i*summary(j)$tTable["no.fn","Value"],
      Intercept.Std.Error=sqrt( (summary(j)$tTable["(Intercept)","Std.Error"]^2+summary(j)$tTable["no.fn","Std.Error"]^2+
                                   2*j$varFix["(Intercept)","no.fn"]) ),
      Estimate=summary(j)$tTable["richness","Value"]+i*summary(j)$tTable["richness:no.fn","Value"],
       #Std.Error=summary(j)$tTable["richness","Std.Error"] )
      Estimate.Std.Error=sqrt( (summary(j)$tTable["richness","Std.Error"]^2+summary(j)$tTable["richness:no.fn","Std.Error"]^2+
                          2*j$varFix["richness","richness:no.fn"]) ) )
    } else {
      data.frame(no.fn=i,Intercept=NA,Intercept.Std.Error=NA,Estimate=NA,Estimate.Std.Error=NA) }    
} ) } )
#Determine which thresholds are significantly different from zero
rawcoefs.df$sig=ifelse(rawcoefs.df$Estimate>2*rawcoefs.df$Estimate.Std.Error,"sig","not.sig")

#Extract threshold of max diversity effect & max threshold at which diversity has a positive effect
maxeffect.df=ddply(rawcoefs.df,"no.fn",function(x) {
  rbind(
    #Find threshold of maximum diversity effect
    cbind(type="max.div.effect",x[which.max(x$Estimate),]),
    #Find maximum threshold at which diversity has a significant positive effect
    cbind(type="max.threshold",x[rev(which(x$sig=="sig"))[1],]) ) } )

#Rename levels for plotting
#levels(maxeffect.df$type)=c("Threshold of maximum diversity effect","Maximum threshold where diversity effect > 0")
maxeffect.df$type=factor(maxeffect.df$type,levels=c("max.div.effect","max.threshold"))

#Plot threshold against effect size with 95% confidence intervals
p2=ggplot(rawcoefs.df,aes(x=threshold,y=Estimate,group=no.fn))+
  geom_hline(xintercept=0,lwd=0.8,lty=1,col="grey70")+
  geom_ribbon(aes(fill=no.fn,ymax=Estimate+2*Estimate.Std.Error,ymin=Estimate-2*Estimate.Std.Error),alpha=0.05)+
  scale_fill_gradient(high="red",low="blue",name="Number of\nfunctions")+
  geom_line(size=1,aes(col=no.fn))+
  scale_color_gradient(high="red",low="blue",name="Number of\nfunctions")+
  labs(x="Threshold",y="Diversity effect (Linear coefficients)")+
  geom_point(data=subset(maxeffect.df,no.fn %in% c(2,12)),
             aes(x=threshold,y=Estimate,shape=as.factor(threshold),fill=no.fn),col="white",size=6.5,width=1.5)+
  scale_shape_manual(values=c(21:24),guide=F)+
  annotate("text",x=0.01,y=Inf,label="b",vjust=1.5,col="black",fontface="bold",size=7)+
  theme_bw(base_size=18)+
  #guides(colour=guide_colourbar(title.position="top"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.direction="horizontal",legend.box="horizontal",legend.position=c(0.275,0.1)); p2

#Generate Figure 2
grid.arrange(arrangeGrob(p1,p2,ncol=2,widths=c(1.5,2)))

#Look at how intercept changes as a function of threshold
ggplot(rawcoefs.df,aes(x=threshold,y=Intercept))+
  geom_ribbon(aes(fill=no.fn,ymax=Intercept+2*Intercept.Std.Error,ymin=Intercept-2*Intercept.Std.Error),
              alpha=0.2)+
  scale_fill_gradient(high="red",low="blue",name="Number of\nfunctions")+
  geom_line(size=1,aes(col=no.fn))+
  scale_color_gradient(high="red",low="blue",name="Number of\nfunctions")+
  facet_wrap(~no.fn,scales="free",nrow=2)+
  labs(x="\nThreshold",y="Intercept\n")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")

#Plot against number of functions
ggplot(maxeffect.df,aes(x=no.fn,y=threshold,group=type,fill=no.fn))+
  #Add trend line
  geom_line()+
  #Scale points by the size of the diversity effect
  geom_point(aes(size=Estimate),shape=21)+
  #geom_text(data=data.frame(type=c("max.threshold","max.div.effect"),no.fn=c(2.3,2.3),threshold=c(0.935,0.795),lab=c("d","c")),
  #          aes(label=lab),size=8,fontface="bold")+
  #Break out panels by response
  facet_wrap(~type,ncol=1,scales="free_y")+
  scale_fill_gradient(high="red",low="blue",name="Number of\nfunctions")+
  scale_size(range=c(4,8),name="Effect size")+
  labs(x="\nNumber of functions",y="Threshold\n")+
  #Specify axis breaks
  scale_x_continuous(breaks=c(2,4,6,8,10,12))+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background=element_blank(),strip.text=element_blank())

#Break out into separate plots
ggplot(subset(maxeffect.df,type=="max.div.effect"),aes(x=no.fn,y=threshold,fill=no.fn))+
  #Add trend line
  geom_line()+
  #Scale points by the size of the diversity effect
  geom_point(aes(size=Estimate),shape=21)+
  scale_fill_gradient(high="red",low="blue",name="Number of\nfunctions",guide=F)+
  scale_size(range=c(4,8),name="Effect size")+
  labs(x="\nNumber of functions",y="Threshold of maximum diversity effect\n")+
  #Specify axis breaks
  scale_x_continuous(breaks=c(2,4,6,8,10,12))+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background=element_blank(),strip.text=element_blank(),
        legend.position=c(0.85,0.24))

ggplot(subset(maxeffect.df,type=="max.threshold"),aes(x=no.fn,y=threshold,fill=no.fn))+
  #Add trend line
  geom_line()+
  #Scale points by the size of the diversity effect #Scale points by the size of the diversity effect
  geom_point(aes(size=Estimate),shape=21)+
  scale_fill_gradient(high="red",low="blue",name="Number of\nfunctions",guide=F)+
  scale_size(range=c(4,8),name="Effect size")+
  labs(x="\nNumber of functions",y="Threshold at which\ndiversity effect is not significant")+
  #Specify axis breaks
  scale_x_continuous(breaks=c(2,4,6,8,10,12))+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background=element_blank(),strip.text=element_blank(),
        legend.position=c(0.85,0.28))

#######################################################################################################

#Extract slopes for each individual experiment from mixed model

#Create dataframe for predicted values for each study
predictindividual.df=ldply(1:99,.progress="text",function(i) {
  ldply(unique(thresholds.df$Study),function(j) {
    #Generate dataframe for predicted values
    newdata=expand.grid(
      Study=j,
      threshold=i/100,
      richness=1:max(subset(thresholds.df,Study==j)[!is.na(subset(thresholds.df,Study==j)$no.fn.greater),"richness"]),
      no.fn=max(subset(thresholds.df,Study==j)$no.fn))
    #Add predicted values with varying slope of richness for each study
    newdata$no.fn.greater=predict(rawmods.list[[i]],newdata,type="response")
    #Return dataframe
    return(newdata)   
} ) } )

#Add column for reference (for panel headers)
predictindividual.df$Reference=multifunc[match(predictindividual.df$Study,multifunc$Study),"Reference"]
#Generate muscle plots
ggplot(predictindividual.df,aes(x=richness,y=no.fn.greater,col=threshold,group=threshold))+
  geom_line(lwd=1,alpha=0.6)+
  scale_color_gradientn(colours=rev(rainbow(5)),name="Threshold")+
  labs(x="\nRichness",y="Number of functions > than threshold\n")+
  facet_wrap(~Reference,scales="free")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),strip.text.x=element_text(size=12),
        legend.position="right" )

#Extract coefficients for muscle plots
individualcoefs.df=ldply(1:99,.progress="text",function(i) {
  ldply(unique(thresholds.df$Study),function(j) {
    ldply(unique(thresholds.df[thresholds.df$Study==j,"no.fn"]),function(k) {
      data.frame(
        Study=rownames(coef(rawmods.list[[i]]))[j],
        no.fn=k,
        threshold=i/100,
        Intercept=coef(rawmods.list[[i]])[j,1]+k*coef(rawmods.list[[i]])[j,3],
        Estimate=coef(rawmods.list[[i]])[j,2]+k*coef(rawmods.list[[i]])[j,4] )
} ) } ) } )

#Add column for reference (for panel headers)
individualcoefs.df$Reference=multifunc[match(individualcoefs.df$Study,multifunc$Study),"Reference"]
#Plot coefs against threshold
ggplot(individualcoefs.df,aes(x=threshold,y=Estimate,group=no.fn,col=no.fn))+
  geom_hline(xintercept=0,lwd=0.8,lty=1,col="grey70")+
  geom_line(lwd=1.2)+
  facet_wrap(~Reference,scales="free",ncol=5)+
  scale_color_gradient(high="red",low="blue",name="Number of\nfunctions")+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","0.5","1"))+
  labs(x="Threshold",y="Diversity Effect")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text.x=element_text(size=12),strip.background=element_blank(),
        legend.direction="horizontal",legend.position=c(0.55,0.025))  
#Plot intercept against threshold
ggplot(individualcoefs.df,aes(x=threshold,y=Intercept,group=no.fn,col=no.fn))+
  geom_hline(xintercept=0,lwd=0.8,lty=1,col="grey70")+
  geom_line(lwd=1.2)+
  facet_wrap(~Reference,scales="free")+
  scale_color_gradient(high="red",low="blue",name="Number of\nfunctions")+
  labs(x="\nThreshold",y="Intercept\n")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),strip.text.x=element_text(size=12))  

#######################################################################################################

#Explore whether this relationship changes as a function of system or trophic level
thresholds.df$FTG=relevel(thresholds.df$FTG,"Primary Producer")
globalmods.list=dlply(thresholds.df,"threshold",.progress="text",function(i) {
  #Set lmeControl for certain thresholds
  if(i$threshold %in% c(0.46)) control=lmeControl() else 
    control=lmeControl(opt="optim",msTol=1e-6)
  #Function to run models
  f=function(x) glmmPQL(no.fn.greater~richness*no.fn+richness*Sys1+richness*FTG,random=~richness|Study,
                        family=quasipoisson(link="identity"),
                        start=c(1,0.1,0.5,rep(0,11)),
                        control=control,
                        verbose=F,data=x)
  safef=failwith(NA,f)
  safef(i)
} ) 

#Create dataframe for predicted values for each system
predictglobal.df=ldply(1:99,.progress="text",function(i) {
  ldply(unique(thresholds.df$Sys1),function(j) { 
    ldply(unique(thresholds.df$FTG),function(k) { 
      #Subset dataframe
      data=subset(thresholds.df,Sys1==j & FTG==k)
      data=data[!is.na(data$no.fn.greater),]
      if("glmmPQL" %in% class(globalmods.list[[i]]))  {
        if(nrow(data)==0) {
          data.frame()
        } else {
          #Create new dataframe to store predictions
          newdata=expand.grid(
            Sys1=j,
            FTG=k,
            threshold=i/100,
            richness=1:max(data$richness),
            no.fn=max(data$no.fn))
          #Generate predicted values
          newdata$no.fn.greater=predict(globalmods.list[[i]],newdata,type="response",level=0)
          return(newdata) }
      } else { 
        data.frame()
} } ) } ) } )
#Plot predicted results by trophic group and system
ggplot(predictglobal.df,aes(x=richness,y=no.fn.greater,col=threshold,group=threshold))+
  geom_line(lwd=1,alpha=0.6)+
  scale_color_gradientn(colours=rev(rainbow(5)),name="Threshold")+
  labs(x="\nRichness",y="Number of functions > than threshold\n")+
  facet_wrap(Sys1~FTG,scales="free",nrow=2)+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position="right" )

#Extract coefficients for muscle plots
globalcoefs.df=ldply(1:99,.progress="text",function(i) {
    ldply(unique(thresholds.df$Sys1),function(k) {
     ldply(unique(thresholds.df$FTG),function(l) {
       data=subset(thresholds.df,Sys1==k & FTG==l)
       if(nrow(data)==0) data.frame() else {
         ldply(2:max(data$no.fn),function(j) {
           if("glmmPQL" %in% class(globalmods.list[[i]])) {
            mod=globalmods.list[[i]]
            data.frame(
              threshold=i/100,
              Sys1=k,
              FTG=l,
              no.fn=j,
              Intercept=
                summary(mod)$tTable["(Intercept)","Value"]+
                ifelse(k=="Aquatic",0,summary(mod)$tTable[paste("Sys1",k,sep=""),"Value"])+
                ifelse(l=="Primary Producer",0,summary(mod)$tTable[paste("FTG",l,sep=""),"Value"])+
                j*summary(mod)$tTable["no.fn","Value"],
              Estimate=
                summary(mod)$tTable["richness","Value"]+
                ifelse(k=="Aquatic",0,summary(mod)$tTable[paste("richness:Sys1",k,sep=""),"Value"])+
                ifelse(l=="Primary Producer",0,summary(mod)$tTable[paste("richness:FTG",l,sep=""),"Value"])+
                j*summary(mod)$tTable["richness:no.fn","Value"] ) 
            } else {
             data.frame() }   
} ) } } ) } ) } )

globalcoefs.df$FTG=factor(globalcoefs.df$FTG,levels=c("Dead organic matter","Detritivore","Primary Producer",
                                                      "Herbivore","Carnivore"))
levels(globalcoefs.df$FTG)=c("Dead\norganic matter","Detritivore","Primary\nproducer","Herbivore","Carnivore")

#Plot coefs against threshold (subset out max number of functions): 12" x 4"
ggplot(ddply(subset(globalcoefs.df,FTG!="Carnivore"),c("Sys1","FTG"),function(x) subset(x,no.fn==max(no.fn))),
       aes(x=threshold,y=Estimate,col=Sys1,group=Sys1))+
  geom_rect(data=data.frame(
    FTG=levels(globalcoefs.df$FTG)[-5],
    Sys1=rep(levels(globalcoefs.df$Sys1),each=4),
    threshold=0,Estimate=0),
    aes(fill=FTG),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.15,show_guide=F)+
  scale_fill_manual(values=c("darkorange4","grey30","deepskyblue3","forestgreen"),guide="none")+
  geom_hline(xintercept=0,lwd=0.8,lty=1,col="grey30")+
  geom_line(lwd=1)+#aes(lty=Sys1),lwd=1)+
  #scale_linetype_manual(values=c(1,6))+
  scale_color_manual(values=c("blue2","darkgreen"),guide=guide_legend(ncol=1),name="")+
  facet_grid(~FTG,scales="free",space="free")+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","0.5","1"))+
  geom_text(data=data.frame(
    FTG=levels(globalcoefs.df$FTG)[-5],
    Sys1=rep(levels(globalcoefs.df$Sys1),each=4),
    labels=letters[1:4]),
    aes(x=0.1,y=Inf,label=labels),vjust=1.5,col="black",fontface="bold",size=6)+
  geom_text(data=data.frame(
    FTG=levels(globalcoefs.df$FTG)[-5],
    Sys1=rep(levels(globalcoefs.df$Sys1),each=4)),
    aes(x=0.15,y=-0.2,label=FTG),vjust=0,hjust=0,col="black",size=5)+
  labs(x="Threshold",y="Diversity Effect")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text.x=element_blank(),strip.background=element_blank(),
        legend.background=element_blank(),legend.key=element_blank(),
        legend.direction="horizontal",legend.box="horizontal",legend.position=c(0.08,0.77))

#Repeat but just for carnivores: 7" x 6"
ggplot(ddply(subset(globalcoefs.df,FTG=="Carnivore"),c("Sys1","FTG"),function(x) subset(x,no.fn==max(no.fn))),
       aes(x=threshold,y=Estimate,col=Sys1,group=Sys1))+
  geom_rect(data=data.frame(
    FTG=levels(globalcoefs.df$FTG)[-5],
    Sys1=rep(levels(globalcoefs.df$Sys1),each=4),
    threshold=0,Estimate=0),
    fill="firebrick1",xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.025,show_guide=F)+
  geom_hline(xintercept=0,lwd=0.8,lty=1,col="grey30")+
  geom_line(lwd=1)+#aes(lty=Sys1),lwd=1)+
  #scale_linetype_manual(values=c(1,6))+
  scale_color_manual(values=c("blue2","darkgreen"),guide=guide_legend(ncol=1),name="")+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","0.5","1"))+
#   geom_text(data=data.frame(
#     FTG=levels(globalcoefs.df$FTG)[5],
#     Sys1=rep(levels(globalcoefs.df$Sys1),each=1)),
#     aes(x=0.15,y=-0.2,label=FTG),vjust=0,hjust=0,col="black",size=5)+
  labs(x="Threshold",y="Diversity Effect")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text.x=element_blank(),strip.background=element_blank(),
        legend.background=element_blank(),legend.key=element_blank(),
        legend.direction="horizontal",legend.box="horizontal",legend.position=c(0.2,0.77))

#Plot intercepts against threshold
ggplot(globalcoefs.df,aes(x=threshold,y=Intercept,col=no.fn,group=no.fn))+
  geom_hline(xintercept=0,lwd=0.8,lty=1,col="grey70")+
  geom_line(lwd=1.2)+
  scale_color_gradient(high="red",low="blue",name="Number of\nfunctions")+
  facet_grid(Sys1~FTG,scales="free")+
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1))+
  labs(x="\nThreshold",y="Intercept\n")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text.x=element_text(size=12),legend.position=c(0.075,0.25) )

#######################################################################################################

#Look at sensitivity of trophic levels results to Smax
#Identify studies with Smax > 6
ref.remove=unique(
  c(as.character(multifunc[multifunc$Smax > 6,"Reference"]),
    as.character(multifunc[!multifunc$FTG %in% c("Primary Producer","Herbivore"),"Reference"])) )
#Subset thresholds.sub.df to remove all experiments that manipulated Smax > 6
thresholds.sub.df=subset(thresholds.df,!Reference %in% ref.remove)
thresholds.sub.df$FTG=relevel(thresholds.sub.df$FTG,"Primary Producer")
#Re-run GLMMs
globalmods.sub.list=dlply(thresholds.sub.df,"threshold",.progress="text",function(i) {
  #Function to run models
  f=function(x) glmmPQL(no.fn.greater~richness*no.fn+richness*Sys1+richness*FTG,random=~richness|Study,
                        family=quasipoisson(link="identity"),
                        start=c(1,0.1,0.5,rep(0,5)),
                        control=lmeControl(opt="optim",msTol=1e-4),
                        verbose=F,data=x)
  safef=failwith(NA,f)
  safef(i)
} ) 

#Extract coefficients for muscle plots
globalcoefs.sub.df=ldply(1:99,.progress="text",function(i) {
  ldply(unique(thresholds.sub.df$Sys1),function(k) {
    ldply(unique(thresholds.sub.df$FTG),function(l) {
      data=subset(thresholds.sub.df,Sys1==k & FTG==l)
      if(nrow(data)==0) data.frame() else {
        ldply(2:max(data$no.fn),function(j) {
          if("glmmPQL" %in% class(globalmods.sub.list[[i]])) {
            mod=globalmods.sub.list[[i]]
            data.frame(
              threshold=i/100,
              Sys1=k,
              FTG=l,
              no.fn=j,
              Intercept=
                summary(mod)$tTable["(Intercept)","Value"]+
                ifelse(k=="Aquatic",0,summary(mod)$tTable[paste("Sys1",k,sep=""),"Value"])+
                ifelse(l=="Primary Producer",0,summary(mod)$tTable[paste("FTG",l,sep=""),"Value"])+
                j*summary(mod)$tTable["no.fn","Value"],
              Estimate=
                summary(mod)$tTable["richness","Value"]+
                ifelse(k=="Aquatic",0,summary(mod)$tTable[paste("richness:Sys1",k,sep=""),"Value"])+
                ifelse(l=="Primary Producer",0,summary(mod)$tTable[paste("richness:FTG",l,sep=""),"Value"])+
                j*summary(mod)$tTable["richness:no.fn","Value"] ) 
          } else {
            data.frame() } }   
        ) }
    } ) } ) } )

#Plot coefs against threshold (subset out max number of functions): 9" x 5"
ggplot(
  data=ddply(subset(globalcoefs.sub.df,FTG %in% c("Primary Producer","Herbivore")),c("Sys1","FTG"),function(x) subset(x,no.fn==max(no.fn))),
  aes(x=threshold,y=Estimate,col=Sys1,group=Sys1))+
  geom_rect(data=data.frame(
    FTG=levels(globalcoefs.sub.df$FTG)[c(1,5)],
    Sys1=rep(levels(globalcoefs.sub.df$Sys1),each=4),
    threshold=0,Estimate=0),
    aes(fill=FTG),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.15,show_guide=F)+
  scale_fill_manual(values=c("deepskyblue3","forestgreen"),guide="none")+
  geom_hline(xintercept=0,lwd=0.8,lty=1,col="grey30")+
  geom_line(lwd=1)+#aes(lty=Sys1),lwd=1)+
  #scale_linetype_manual(values=c(1,6))+
  scale_color_manual(values=c("blue2","darkgreen"),guide=guide_legend(ncol=1),name="")+
  facet_grid(~FTG,scales="free",space="free")+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","0.5","1"))+
  geom_text(data=data.frame(
    FTG=levels(globalcoefs.sub.df$FTG)[c(1,5)],
    Sys1=rep(levels(globalcoefs.sub.df$Sys1),each=4),
    labels=letters[1:2]),
    aes(x=0.1,y=Inf,label=labels),vjust=1.5,col="black",fontface="bold",size=6)+
  geom_text(data=data.frame(
    FTG=levels(globalcoefs.sub.df$FTG)[c(1,5)],
    Sys1=rep(levels(globalcoefs.sub.df$Sys1),each=4)),
    aes(x=0.15,y=-0.2,label=FTG),vjust=0,hjust=0,col="black",size=5)+
  labs(x="Threshold",y="Diversity Effect")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text.x=element_blank(),strip.background=element_blank(),
        legend.background=element_blank(),legend.key=element_blank(),
        legend.direction="horizontal",legend.box="horizontal",legend.position=c(0.12,0.77))

#######################################################################################################

#Simulations to explore the dip in diversity effect around 50% threshold for high numbers of functions
#Extract data for studies that exhibit the dip for further simulations
cedarcreek.df=subset(multifunc.long,Reference=="Cedar Creek (E120)" & value!="NA")
wardle.df=subset(multifunc.long,Reference=="Wardle et al. 2003" & value!="NA")

sim.df=ldply(1:100,.progress="text",function(rep) {
  ldply(seq(0,1,0.1),function(i) { #Percent of negative functions
    ldply(list(cedarcreek.df,wardle.df),function(j) {
      #Specify number of functions
      no.fn=length(unique(j$Ydesc))
      #Specify diversity levels
      divlevels=j[j$Ydesc %in% unique(j$Ydesc)[1],"richness"]
      #Construct data.frame with appropriate dims
      newdf=cbind(
        data.frame(diversity=divlevels),
        matrix(rep(NA,length(divlevels)*no.fn),nrow=length(divlevels)) )
      #Populate with values
      newdf[,2:ncol(newdf)]=colwise(function(x) rnorm(nrow(newdf),newdf$diversity,1))(newdf[,2:ncol(newdf)])
      #Define % of functions that have a negative relationship with diversity
      negcol=round(no.fn*i)
      #Convert that number in newdf to negative
      if(negcol>0) newdf[,2:(negcol+1)]=-newdf[,2:(negcol+1)]+max(newdf[,2:(negcol+1)]) else newdf=newdf
      #Scale response
      #newdf[,2:ncol(newdf)]=colwise(function(x) x/max(x))(newdf[,2:ncol(newdf)])
      #Generate threshold data
      thresh.df=ldply(1:99/100,function(k) 
        cbind(
          diversity=newdf[,1],
          thresh=k,
          no.fn.greater=rowSums(colwise(function(x) x>=max(x)*k)(newdf[,2:ncol(newdf)])) ) )
      #Get linear slopes
      ddply(thresh.df,"thresh",function(x) {
        mod=try(glm(no.fn.greater~diversity,data=x,family=quasipoisson(link="identity"),start=c(0,0.1)))
        if(class(mod) == "try-error") 
          cbind(
            rep=rep,
            Reference=as.character(unique(j$Reference)),
            pneg=i,
            no.fn=no.fn,
            thresh=unique(x$thresh),
            Estimate=NA,
            Std.Error=NA,
            N=NA) else
          cbind(
            rep=rep,
            Reference=as.character(unique(j$Reference)),
            pneg=i,
            no.fn=no.fn,
            thresh=unique(x$thresh),
            Estimate=summary(mod)$coefficients["diversity",1],
            Std.Error=summary(mod)$coefficients["diversity",2],
            n=mod$df.residuals) } )
    } )
  } )
} )

sim.df[,3:7]=apply(sim.df[,3:7],2,function(x) as.numeric(x))

#Summarize for plotting
sim.df.summary=ddply(sim.df,c("Reference","pneg","no.fn","thresh"),function(x) 
  data.frame(
    Reference=unique(x$Reference),
    pneg=unique(x$pneg),
    no.fn=unique(x$no.fn),
    thresh=unique(x$thresh),
    Estimate=mean(x$Estimate),
    Std.Error=std.error(x$Estimate)) )
#    Std.Error=sqrt(sum(x$Std.Error/x$n)) ) )

#Graph results
ggplot(sim.df,aes(x=thresh,y=Estimate,group=rep))+
  geom_hline(yintercept=0,col="grey50",lwd=0.5)+
#  geom_ribbon(aes(ymin=Estimate-Std.Error,ymax=Estimate+Std.Error),col="blue")+
  geom_line(lwd=0.3,alpha=0.7)+
  scale_x_continuous(breaks=c(0,0.5,1),labels=c("0","0.5","1"))+
  #scale_color_gradient(high="green",low="firebrick1",name="Number of\nfunctions")+
  facet_grid(Reference~pneg,scales="free")+
  labs(x="Threshold",y="Diversity effect")+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#######################################################################################################

#Conduct sensitivity analysis by removing each study individually, and re-generating coef plot
rawmods.sensitivity.df=ldply(c("0.2","0.4","0.6","0.8"),function(i) {
  ldply(unique(thresholds.df$Study),.progress="text",function(j) {
    #Subset data
    data=subset(thresholds.df,Study!=j & threshold==i)
    #Function to run models
    f=function(x) glmmPQL(no.fn.greater~richness*no.fn,random=~richness|Study,
                          family=quasipoisson(link="identity"),start=c(0.15,0.05,0.1,0.001),            
                          control=lmeControl(opt="optim",msTol=1e-4),
                          verbose=F,data=x)
    #Run model, return NA if model fails
    safef=failwith(NA,f)
    mod=safef(data)
    #Extract coefficient and standard error and return in data.frame
    if(is.na(mod)) data.frame(Study.removed=j,threshold=i,coef=NA,coef.se=NA) else
      data.frame(
        Study.removed=j,
        threshold=i,
        coef=summary(mod)$tTable[2,1],
        coef.se=summary(mod)$tTable[2,2] )
  } )
} ) 

rawmods.sensitivity.df$Ref.removed=thresholds.df[match(rawmods.sensitivity.df$Study,thresholds.df$Study),"Reference"]
rawmods.sensitivity.df$threshold=factor(rawmods.sensitivity.df$threshold)
levels(rawmods.sensitivity.df$threshold)=c("20%","40%","60%","80%")

#Plot results: 12" x 8"
ggplot(rawmods.sensitivity.df,aes(y=coef,x=Ref.removed))+
  geom_hline(data=data.frame(
    threshold=c("20%","40%","60%","80%"),
    coef=c(summary(rawmods.list[[20]])$tTable[2,1],
           summary(rawmods.list[[40]])$tTable[2,1],
           summary(rawmods.list[[60]])$tTable[2,1],
           summary(rawmods.list[[80]])$tTable[2,1]) ),
    aes(yintercept=coef),lwd=1,alpha=0.6,lty=1)+  
  geom_point(size=3)+
  geom_errorbar(aes(ymax=coef+2*coef.se,ymin=coef-2*coef.se))+
  coord_flip()+
  facet_wrap(~threshold,nrow=1)+
  labs(y="Regression coefficient",x="Study removed")+
  theme_bw(base_size=12)+
  theme(legend.position="none",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#######################################################################################################
#                                          TURNOVER APPROACH                                          #
#######################################################################################################

#First, extract only mono data for each experiment and store in a list
mono.list=dlply(multifunc,c("Reference","Study","Expt","Smax","Sys1","FTG"),function(i) {
  #Get monocultures and metadata
  monos=i[,c(2:4,9,grep("mono",colnames(i)))]
  #Melt monos dataframe
  monos.melt=melt(monos,id.vars=c(1:4),measures.vars=c(5:ncol(monos)))
  #Split columns based on response (Y, N, SD, or ID)
  monos.melt=cbind(monos.melt[1:4],
                   t(matrix(unlist(strsplit(gsub("([0-9]+)","\\1~",monos.melt$variable),"~" )),nrow=2)),
                   value=monos.melt[,"value"])
  names(monos.melt)[5:6]=c("treatment","variable")
  #Cast variables
  monos.cast=dcast(monos.melt,Study+Expt+Reference+Ydesc+treatment~variable,value.var="value")
  #Remove rows where Y==NA
  monos.cast=monos.cast[!is.na(monos.cast$Y),]
  #Convert responses to numeric
  monos.cast[,c("Y","SD","N")]=apply(monos.cast[,c("Y","SD","N")],2,as.numeric)
  #Return dataframe
  return(monos.cast)
} )

#Determine the most extreme species for each function, then sum the number of unique most extreme species 
#across all functions within an experiment
extremesp.df=ldply(mono.list,function(i) {
  data.frame(no.sp=length(unique(i$treatment)),no.fn=length(unique(i$Ydesc)),
             no.max.sp=length(unique(ddply(i,"Ydesc",function(x) x[which.max(x$Y),"ID"])$V1)),
             no.min.sp=length(unique(ddply(i,"Ydesc",function(x) x[which.min(x$Y),"ID"])$V1)) ) } )

#Mean turnover in extreme species across all functions measured within an experiment
mean(extremesp.df$no.max.sp/extremesp.df$no.fn); std.error(extremesp.df$no.max.sp/extremesp.df$no.fn); range(extremesp.df$no.max.sp/extremesp.df$no.fn)
mean(extremesp.df$no.min.sp/extremesp.df$no.fn); std.error(extremesp.df$no.min.sp/extremesp.df$no.fn); range(extremesp.df$no.min.sp/extremesp.df$no.fn)

#Parse by system and trophic level
ddply(extremesp.df,"Sys1",function(x) data.frame(
  max.effect.size=mean(x$no.max.sp/x$no.fn),
  max.std.error=std.error(x$no.max.sp/x$no.fn),
  min.effect.size=mean(x$no.min.sp/x$no.fn),
  min.std.error=std.error(x$no.min.sp/x$no.fn) ) )

ddply(extremesp.df,"FTG",function(x) data.frame(
  max.effect.size=mean(x$no.max.sp/x$no.fn),
  max.std.error=std.error(x$no.max.sp/x$no.fn),
  min.effect.size=mean(x$no.min.sp/x$no.fn),
  min.std.error=std.error(x$no.min.sp/x$no.fn)  ) )

#######################################################################################################
#                                      MULTIPLICATIVE APPROACH                                        #
#######################################################################################################

#Calculate multiplicative level of functioning across all functions for each treatment, for each experiment
multifunc.mult=ddply(multifunc.scaled,c("Reference","Study","Expt","FTG","Sys1","Sys2"),function(x) {
  z=data.frame(
    richness=colnames(x[,Y.colnames]),
    no.fn=length(unique(x$Ydesc)),
    mult.fn=apply(x[,Y.colnames],2,prod),
    mult.fn.N=colSums(x[,N.colnames]) ) 
  #Remove rows where there is no response (i.e., mult.fn.scaled==NA)
  z=z[!is.na(z$mult.fn),]
  #Scale by nth root 
  z$mult.fn.scaled=z$mult.fn^(1/nrow(x))
  #Set richness by splitting column names
  z$richness=suppressWarnings(ifelse(grepl("mono",z$richness),1,as.numeric(gsub("X([0-9]+).*","\\1",z$richness))))
  return(z) } )

#Investigate proper functional form to use
#Group data for random effects
multifunc.mult.grouped=groupedData(mult.fn.scaled~richness|Study,data=multifunc.mult)
#Fit different functional forms using non-linear mixed models
Null=nlme(mult.fn.scaled~a,fixed=a~1,random=~a~1,start=c(a=0.2),data=multifunc.mult.grouped)
Linear=nlme(mult.fn.scaled~a+b*richness,fixed=a+b~1,random=~a+b~1,start=c(a=1.5,b=1),data=multifunc.mult.grouped)
Logarithmic=nlme(mult.fn.scaled~a+b*log(richness),fixed=a+b~1,random=~a+b~1,start=c(a=1.5,b=1),data=multifunc.mult.grouped)
Power=nlme(mult.fn.scaled~a*richness^b,fixed=a+b~1,random=~a+b~1,start=c(a=0.2,b=2),data=multifunc.mult.grouped)
Saturating=nlme(mult.fn.scaled~richness/(k+richness),fixed=k~1,random=k~1,start=c(k=1),data=multifunc.mult.grouped)
#Compare models using AIC
AIC(Null,Linear,Logarithmic,Power,Saturating)

#Fit log relationship using linear mixed effects model, allowing slopes and intercepts to vary by Study
multmods.list=lapply(c("unweighted"),function(i) {
  #Subset dataset to include only non-NA data points for each type of analysis
  if(i=="variance") { multifunc.mult=multifunc.mult[!is.na(multifunc.mult$mult.fn.scaled.SD),] 
  } else if(i=="sample.size") { multifunc.mult=multifunc.mult[!is.na(multifunc.mult$mult.fn.scaled.N),]
  } else { multifunc.mult }
  #Fit linear mixed effects model for each weighting scheme
  if(i=="unweighted") { 
    mod=glmmPQL(mult.fn.scaled~log(richness),random=~richness|Study,family=quasibinomial(link="identity"),
                start=c(0.5,0),data=multifunc.mult,verbose=F)
  } else if(i=="variance") { 
    mod=glmmPQL(mult.fn.scaled~log(richness),random=~richness|Study,weights=1/((multifunc.mult$mult.fn.scaled.SD^2)+0.01),
                family=quasibinomial(link="identity"),start=c(0.5,0),data=multifunc.mult,verbose=F)
  } else { 
    mod=glmmPQL(mult.fn.scaled~log(richness),random=~richness|Study,weights=sqrt(multifunc.mult$mult.fn.scaled.N),
                family=quasibinomial(link="identity"),start=c(0.5,0),data=multifunc.mult,verbose=F) }
  #Return model
  return(mod)
} )
#Append reduced model (S <= 16)
multmods.list=append(multmods.list,list(update(multmods.list[[1]],data=subset(multifunc.mult,richness<=16))))
#Look at output and diagnostic plots
lapply(multmods.list,summary); lapply(multmods.list,plot)

#Extract predicted fits for plot
pred.df.list=lapply(seq_along(multmods.list),function(i) {
  if(i==1) multifunc.mult=multifunc.mult else multifunc.mult=subset(multifunc.mult,!paste(Study,Expt) %in%    
                                                                      unique(paste(subset(multifunc.mult,richness>16)$Study,subset(multifunc.mult,richness>16)$Expt)))
  #Modified from: http://glmm.wikidot.com/faq 
  #Create dataframe for predicted values for overall fit
  newdata=expand.grid(richness=1:max(multifunc.mult$richness),no.fn=2:max(multifunc.mult$no.fn),mult.fn.scaled=0)
  #Generate predicted values for overall trend
  newdata$mult.fn.scaled=predict(multmods.list[[i]],newdata,type="response",level=0)
  #Obtain model matrix
  mm=model.matrix(terms(multmods.list[[i]]),newdata)
  #Obtain estimate of SE based on fixed effects only
  newdata$fixedSE=sqrt(diag(mm %*% tcrossprod(vcov(multmods.list[[i]]),mm)))
  
  #Create dataframe for predicted values for each study
  newdata2=ldply(unique(multifunc.mult$Study),function(j)
    data.frame(Study=j,
               richness=1:max(subset(multifunc.mult,Study==j)$richness),
               no.fn=max(subset(multifunc.mult,Study==j)$no.fn)) )
  #Generate predicted values for each study
  newdata2$mult.fn.scaled=predict(multmods.list[[i]],newdata2,type="response",level=1)
  
  #Return dataframes in a list
  list(newdata,newdata2)  
} )     

#Plot predicted values and confidence bands across all studies on top of fitted values for each individual study
avgplots.list=lapply(1,function(i) {
  ggplot()+
    #Plot raw points
    geom_point(data=multifunc.mult,aes(x=richness,y=mult.fn.scaled),size=2,col="grey60",alpha=0.5)+
    #Plot curves for each study
    geom_line(data=pred.df.list[[i]][[2]],aes(x=richness,y=mult.fn.scaled,group=Study),col="black",lwd=0.8,alpha=0.7)+
    #Add confidence band for mixed model predictions based on fixed effects only
    geom_ribbon(data=pred.df.list[[i]][[1]],
                aes(x=richness,y=mult.fn.scaled,ymax=mult.fn.scaled+2*fixedSE,ymin=mult.fn.scaled-2*fixedSE),
                fill="red",alpha=0.4)+
    #Add line for mixed model predictions
    geom_line(data=pred.df.list[[i]][[1]],aes(x=richness,y=mult.fn.scaled),col="red",lwd=1.5,alpha=0.9)+
    scale_x_continuous(breaks=c(1,20,40,60))+
    scale_y_continuous(limits=c(0,1.1),breaks=c(0,0.5,1))+
    labs(x="Richness",y="Average multifunctionality")+
    theme_bw(base_size=18)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
} )#; avgplots.list

#Plot figure 1 with inset from list above (5 x 4.5")
ggplot(data=subset(multifunc.mult,!paste(Study,Expt) %in%    
                     unique(paste(subset(multifunc.mult,richness>16)$Study,subset(multifunc.mult,richness>16)$Expt))),
       aes(x=richness,y=mult.fn.scaled))+
  #Plot raw points
  geom_point(size=2.5,col="grey60",alpha=0.3)+
  #Plot curves for each study
  geom_line(data=pred.df.list[[2]][[2]],aes(x=richness,y=mult.fn.scaled,group=Study),col="grey20",lwd=0.75,alpha=0.8)+
  #Add confidence band for mixed model predictions based on fixed effects only
  geom_ribbon(data=pred.df.list[[2]][[1]],
              aes(x=richness,y=mult.fn.scaled,ymax=mult.fn.scaled+2*fixedSE,ymin=mult.fn.scaled-2*fixedSE),
              fill="red",alpha=0.4)+
  #Add line for mixed model predictions
  geom_line(data=pred.df.list[[2]][[1]],aes(x=richness,y=mult.fn.scaled),col="red",lwd=2.5)+
  coord_cartesian(ylim=c(-0.05,1.1))+
  scale_x_continuous(breaks=c(1,4,8,12,16),labels=c("1","4","8","12","16"))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  labs(x="Richness",y="Multiplicative multifunctionality")+ 
  #   geom_text(data=data.frame(
  #     labels=letters[1]),
  #     aes(x=-Inf,y=Inf,label=labels),vjust=1.5,hjust=-1.5,col="black",fontface="bold",size=9)+
  theme_bw(base_size=18)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  annotation_custom(grob=ggplotGrob(avgplots.list[[1]]+labs(x="",y="")+theme(plot.margin=unit(c(0,0,0,0),"cm"))),
                    xmin=7,xmax=16,ymin=-0.075,ymax=0.385)
