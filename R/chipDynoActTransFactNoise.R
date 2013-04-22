#function [list,newX, newXVals]=chipDynoActTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#                                         TransNames, annotation);
#
#% CHIPDYNOACTTRANSFACTNOISE identifies significantly varying TFs.
#%
#%	Description:
#%	[list,newX, newXVals]=chipDynoActTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#%                                         TransNames, annotation);
#%% 	chipDynoActTransFactNoise.m version 1.4

chipDynoActTransFactNoise=function (data,X,Sigma,beta, precs, gamma,mu, TransNames, annotation) {

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));
#i=168 #For test purpose

#nTrans=size(TransNames,1);
#list=[];
#newX=zeros(size(X));
#newXVals=zeros(size(X));

source("chipDynoTransFactNoise.R")
source("chipDynoMaxDiff.R")

for (i in 1: nTrans) {
	expectations =chipDynoTransFactNoise(data,X,Sigma,beta,precs, gamma,mu, TransNames, annotation, TransNames[i,]);
	TF = expectations[[1]]
	TFError = expectations [[2]]
	TFErrorDiff = expectations [[3]]

	maxVars=chipDynoMaxDiff(TF,TFErrorDiff);

	pvals=pnorm(-maxVars) #?????????? Need to chack again!!!!!!!!!!
	sigVars = pvals[which(pvals<0.02)];
       	lst=cbind(lst, length(sigVars));
	index=which(X[,i]!=0);
	newX[index[which(pvals<0.02)],i]=1
	newXVals[index[which(pvals<0.02)],i]=pvals[which(pvals<0.02)];
}

#for i=1:nTrans
#    [TF,TFError,TFErrorDiff]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#                                         TransNames, annotation, ...
#                                        TransNames(i));
#    maxVars=chipDynoMaxDiff(TF,TFErrorDiff);
#    pvals=cumGaussian(-maxVars));
#    sigVars=pvals(find(pvals<0.02));
#    list=[list, size(sigVars,2)];
#    index=find(X(:,i));
#    newX(index(find(pvals<0.02)),i)=1;
#     newXVals(index(find(pvals<0.02)),i)=pvals((find(pvals<0.02)));
#end


f=list(lst,newX, newXVals)
return(f)
}
