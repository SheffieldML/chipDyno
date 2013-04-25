#function [TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, transNames, annotations, name);
#
#% CHIPDYNOTRANSFACT provides gene-specific TFAs with errorbars.
#%
#%	Description:
#%	[TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, ...
#%                                         transNames, annotations, ...
#%                                        name);
#%% 	chipDynoTransFact.R version 0.1.0

####
# For Figure 1
# anno[48] is "YER124C"
# anno[34] is "SCW11"
###


chipDynoTransFact = function(data, X, Sigma, beta, gamma, mu, transNames, annotations, name) {


###
# Only for test purpose
# name= "ACE2" 
# name = TransNames[i,]
# transNames = TransNames
# annotations = annotation
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

source("chipDynoExpectationsFast.R")

for (i in 1 : nTargets) {
	expectations = chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, transNames, annotations, name, genesIn[i]);

	TF[i,] = expectations[[1]];
	TFError[i,] = expectations[[2]];
	TFErrorDiff[ , ,i] = expectations[[3]] ; 									 
}


#for i=1:nTargets
#  [TF(i,:),TFError(i,:),TFErrorDiff(:,:,i)]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
#                                         transNames, annotations, ...
#                                         name,genesIn(i));
#end

#

expectations = list(TF,TFError,TFErrorDiff)
return(expectations)

}
