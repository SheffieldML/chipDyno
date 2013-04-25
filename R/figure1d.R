rm(list=ls())
load("ResultsSpellman_200Ita.RData")

#load("ResultsSpellman_4000Ita_temp.RData")

#%% 	chipDynoActTransFact.R 
 nTrans=nrow(TransNames);
 lst=list();
 newX=array(0, dim <-c(dim(X)));
 newXVals=array(0, dim <-c(dim(X)));
 
 source("chipDynoTransFact.R")
 source("chipDynoMaxDiff.R")
 i=2

#################
#%% 	chipDynoTransFact.R
 
 name=TransNames[i,]
 transNames = TransNames
 annotations = annotation
 
 index=which(name == transNames);
 genesIn=which(X[,index]!=0);
 anno=annotations[which(X[,index]!=0)];
 nTargets=length(anno);
 npts=ncol(data);
 TF=array(0, dim=c(nTargets,npts));
 TFError=array(0, dim=c(nTargets,npts));
 TFErrorDiff=array(0, dim=c(npts,npts,nTargets));
 
 source("chipDynoExpectationsFast.R")
 i=28
 expectations = chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, transNames, annotations, name, genesIn[i]);
 
 TF[i,] = expectations[[1]];
 TFError[i,] = expectations[[2]];
#########################


#plot Errorbar Function
 plotErrorBar <- function(y,SE){
 
 add.error.bars <- function(x,y,SE,w,col=1){
 x0 = x; y0 = (y-SE); x1 =x; y1 = (y+SE);
 arrows(x0, y0, x1, y1, code=3,angle= 90,length=w,col=col);
 }
 
 x <- c(1:length(y));
 plot(x,y, type = 'l', col='green4', las=1);
 add.error.bars(x,y,SE,0.05,col='red');
 
 }

 plotErrorBar(TF[i,],TFError[i,])
