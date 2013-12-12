# CHIPDYNOACTTRANSFACT identifies significantly varying TFs.
# CHIPDYNO toolbox
# chipDynoActTransFact.R version 1.1
# FORMAT chipDynoActTransFact <- function (data,X,Sigma,beta,gamma,mu, 
#                                         TransNames, annotation,sigLev)
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
# RETURN f : concatenated dataframe of list of regulators for a specific gene,
# its index and values
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoTransFact, chipDynoTransFactNoise, chipDynoActTransFactNoise

chipDynoActTransFact <- function (data,X,Sigma,beta,gamma,mu, TransNames, annotation, sigLev) {

# sigLev= 10; # For Tu data Set! unknown!! just for development!!!
# sigLev= 10; # Spellman data Set unknown! just for development!!!

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));

#source("chipDynoTransFact.R")
#source("chipDynoMaxDiff.R")

for (i in 1: nTrans) {
	expectations =chipDynoTransFact(data,X,Sigma,beta,gamma,mu, TransNames, annotation, TransNames[i,]);
	TF = expectations[[1]]
	TFError = expectations [[2]]
	TFErrorDiff = expectations [[3]]

	maxVars=chipDynoMaxDiff(TF,TFErrorDiff);

	sigVars = maxVars[which(maxVars>sigLev)];
       	lst=cbind(lst, length(sigVars));
	index=which(X[,i]!=0);
	newX[index[which(maxVars>sigLev)],i]=1
	newXVals[index[which(maxVars>sigLev)],i]=maxVars[which(maxVars>sigLev)];
}

f=list(lst,newX, newXVals)
return(f)
}
