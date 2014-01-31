rm(list=ls())
load("ResultsSpellman_200Ita.RData")

#load("ResultsTu_New_10000Ita_Lovelace.RData")
#load("/home/muhammad/Dropbox/Rwork_HD/ResultsTu_New_10000Ita_Lovelace.RData")
source("/home/muhammad/Dropbox/Git/mlprojects/chipDyno/R/chipDynoLoadModules.R")
#load("ResultsSpellman_4000Ita_temp.RData")

#%% 	chipDynoActTransFact.R 
 nTrans=nrow(TransNames);
 lst=list();
 newX=array(0, dim <-c(dim(X)));
 newXVals=array(0, dim <-c(dim(X)));
 
 #source("chipDynoTransFact.R")
 #source("chipDynoMaxDiff.R")
 i=3 # ACE2 for demSpellman
 #i= 4 #for demTu

 # i=46 # for LEU3

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
 
 #i=1
 #source("chipDynoExpectationsFast.R")
 #i=28 # for SCW11
 #i=11 # for YER124C
 #i=49 # for MCR1
 #expectations = chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, transNames, annotations, name, genesIn[i]);
 #expectations = chipDynoExpectationsFastNoise(data,X,Sigma,beta, precs, gamma,mu, transNames, annotations, name, genesIn[i]);
 
 #TF[i,] = expectations[[1]];
 #TFError[i,] = expectations[[2]];
#########################


#plot Errorbar Function
 plotErrorBar <- function(y,SE, g_name){
 
 add.error.bars <- function(x,y,SE,w,col=1){
 x0 = x; y0 = (y-SE); x1 =x; y1 = (y+SE);
 arrows(x0, y0, x1, y1, code=3,angle= 90,length=w,col=col);
 }
 
 x <- c(1:length(y));
 plot(x,y, type = 'l', col='green4', las=1, xlab="time", ylab="TFA", main=bquote("Gene : " ~ .(g_name)));
 add.error.bars(x,y,SE,0.05,col='red');
 
 }

 # plotErrorBar(TF[i,],TFError[i,])

##################


M <- matrix(c(rep(1:4)), byrow=TRUE, nrow=2) # Choose the position by matrix setting!
layout(M) 

#jpeg('rplot.jpg')
#bmp("testplot.bmp")
#png(filename="TFA_data_sample1.png")
#k <- c(17,19, 28, 50)
#k <- c(46,47,50,52)
for (i in 1:4){
#for (i in k){

expectations = chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, transNames, annotations, name, genesIn[i]);
#expectations = chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs, gamma,mu, transNames, annotations, name, genesIn[i]);
 
 TF[i,] = expectations[[1]];
 TFError[i,] = expectations[[2]];
 
g_name=annotation[genesIn[i]]

plotErrorBar(TF[i,],TFError[i,], g_name)
#plotErrorBar(S1_TF[i,],S1_TFError[i,],S2_TF[i,],S2_TFError[i,],S3_TF[i,],S3_TFError[i,], actGeneNames[i])
}