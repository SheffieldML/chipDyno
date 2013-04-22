#function [list,newX, newXVals]=chipDynoActTransFact(data,X,Sigma,beta,gamma,mu, TransNames, annotation,sigLev);
#
#% CHIPDYNOACTTRANSFACT identifies significantly varying TFs.
#%
#%	Description:
#%	[list,newX, newXVals]=chipDynoActTransFact(data,X,Sigma,beta,gamma,mu, ...
#%                                         TransNames, annotation,sigLev);
#%% 	chipDynoActTransFact.R version 0.1.0

chipDynoActTransFact=function (data,X,Sigma,beta,gamma,mu, TransNames, annotation, sigLev) {

# sigLev= 10; # For Tu data Set! unknown!! just for development!!!
# sigLev= 10; # Spellman data Set unknown! just for development!!!

#load("ResultsTu_New_500Ita.RData")
#load("ResultsSpellman_200Ita.RData")
#load("/home/muhammad/H-drive/CElegans/Results_cElegans_100Ita.RData")

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));

#nTrans=size(TransNames,1);
#list=list[];
#newX=zeros(size(X));
#newXVals=zeros(size(X));

source("chipDynoTransFact.R")
source("chipDynoMaxDiff.R")

for (i in 1: nTrans) {
	expectations =chipDynoTransFact(data,X,Sigma,beta,gamma,mu, TransNames, annotation, TransNames[i,]);
	TF = expectations[[1]]
	TFError = expectations [[2]]
	TFErrorDiff = expectations [[3]]
	#    %vars=max(abs((TF-mu(i)*ones(size(TF)))'./TFError'));

	maxVars=chipDynoMaxDiff(TF,TFErrorDiff);

	sigVars = maxVars[which(maxVars>sigLev)];
       	lst=cbind(lst, length(sigVars));
	index=which(X[,i]!=0);
	newX[index[which(maxVars>sigLev)],i]=1
	newXVals[index[which(maxVars>sigLev)],i]=maxVars[which(maxVars>sigLev)];
}

#### Plot thw ErrorBar ###
# source("plotErrorBar.R")
# plotErrorBar(TF[1,],TFError[1,]);
#####

#for i=1:nTrans
#    [TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, ...
#                                         TransNames, annotation, ...
#                                        TransNames(i));
#    %vars=max(abs((TF-mu(i)*ones(size(TF)))'./TFError'));
#    maxVars=chipDynoMaxDiff(TF,TFErrorDiff);
#    sigVars=maxVars(find(maxVars>sigLev));
#    list=[list, size(sigVars,2)];
#    index=find(X(:,i));
#    newX(index(find(maxVars>sigLev)),i)=1;
#    newXVals(index(find(maxVars>sigLev)),i)=maxVars((find(maxVars>sigLev)));
#end
f=list(lst,newX, newXVals)
return(f)
}
