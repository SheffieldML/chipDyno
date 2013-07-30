#function [TF,TFError,TFErrorDiff]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#                                         transNames, annotations, ...
#                                        name);#

#% CHIPDYNOTRANSFACTNOISE given a transcription factor, provides TFAs.
#%
#%	Description:
#%	[TF,TFError,TFErrorDiff]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#%                                         transNames, annotations, ...
#%                                        name);
#%% 	chipDynoTransFactNoise.m version 0.1.0


chipDynoTransFactNoise = function(data, X, Sigma, beta, precs, gamma, mu, transNames, annotations, name) {


###
# Only for test purpose
#name= "ZAP1"
#transNames = TransNames
#annotations = annotation
#i=8
# name= TransNames[i]
##

index=which(name == transNames);
genesIn=which(X[,index]!=0);
anno=annotations[which(X[,index]!=0)];
nTargets=length(anno);
npts=ncol(data);
TF=array(0, dim=c(nTargets,npts));
TFError=array(0, dim=c(nTargets,npts));
TFErrorDiff=array(0, dim=c(npts,npts,nTargets));

#index=find(strcmp(name,transNames));
#genesIn=find(X(:,index));
#anno=annotations(find(X(:,index)));
#nTargets=size(anno,1);
#npts=size(data,2);
#TF=zeros(nTargets,npts);
#TFError=zeros(nTargets,npts);
#TFErrorDiff=zeros(npts,npts,nTargets);

source("chipDynoExpectationsFastNoise.R")

for (i in 1 : nTargets) {
	expectations = chipDynoExpectationsFastNoise(data,X,Sigma,beta, precs, gamma,mu, transNames, annotations, name, genesIn[i]);

	TF[i,] = expectations[[1]];
	TFError[i,] = expectations[[2]];
	TFErrorDiff[ , ,i] = expectations[[3]] ; 									 
}


#for i=1:nTargets
#  [TF(i,:),TFError(i,:),TFErrorDiff(:,:,i)]=chipDynoExpectationsFastNoise#(data,X,Sigma,beta,precs,gamma,mu, ...
#                                         transNames, annotations, ...
#                                         name,genesIn(i));
#end


expectations = list(TF,TFError,TFErrorDiff)
return(expectations)

}
