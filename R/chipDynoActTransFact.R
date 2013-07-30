# CHIPDYNOACTTRANSFACT identifies significantly varying TFs.
# CHIPDYNO toolbox
# chipDynoActTransFact.m version 1.1
# FORMAT [list,newX, newXVals]=chipDynoActTransFact(data,X,Sigma,beta,gamma,mu, ...
#                                         TransNames, annotation,sigLev);
# DESC identifies significantly varying TFs.
# ARG data : point estimate of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG Sigma : prior covariance matrix
# ARG beta :
# ARG gamma : degree of temporal continuity
# ARG mu : mean value of the transcription factor activity
# ARG TransNames : Transcription factors
# ARG annotation : Gene names
# ARG sigLev : threshold value
# RETURN f[[1]] : (list) list of regulators for a specific gene
# RETURN f[[2]] : (newX) 
# RETURN f[[1]] : (newXVals)
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoTransFact, chipDynoTransFactNoise, chipDynoActTransFactNoise

chipDynoActTransFact=function (data,X,Sigma,beta,gamma,mu, TransNames, annotation, sigLev) {

# sigLev= 10; # For Tu data Set! unknown!! just for development!!!
# sigLev= 10; # Spellman data Set unknown! just for development!!!

#load("/home/muhammad/H-drive/CElegans/Results_cElegans_Optim_Sample1.RData")

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
