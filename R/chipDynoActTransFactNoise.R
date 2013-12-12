# CHIPDYNOACTTRANSFACTNOISE identifies significantly varying TFs with uncertainty of
# expression level.
# CHIPDYNO toolbox
# chipDynoActTransFactNoise.R version 1.1
# FORMAT chipDynoActTransFactNoise <- function (data,X,Sigma,beta, precs, gamma,mu,
#                                         TransNames, annotation,sigLev)
# DESC identifies significantly varying TFs.
# ARG data : point estimate of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG Sigma : prior covariance matrix of TFA
# ARG beta :
# ARG precs : uncertainty of the expression level 
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
# SEEALSO : chipDynoTransFact, chipDynoTransFactNoise, chipDynoActTransFact

chipDynoActTransFactNoise <- function (data,X,Sigma,beta, precs, gamma,mu, 
										TransNames, annotation) {

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));

#source("chipDynoTransFactNoise.R")
#source("chipDynoMaxDiff.R")

for (i in 1: nTrans) {
	expectations =chipDynoTransFactNoise(data,X,Sigma,beta,precs, gamma, mu, 
						TransNames, annotation, TransNames[i,]);
	TF = expectations[[1]]
	TFError = expectations [[2]]
	TFErrorDiff = expectations [[3]]

	maxVars=chipDynoMaxDiff(TF,TFErrorDiff);

	pvals=ecdf(-maxVars) #?????????? Need to chack again!!!!!!!!!!
	sigVars = pvals[which(pvals<0.02)];
       	lst=cbind(lst, length(sigVars));
	index=which(X[,i]!=0);
	newX[index[which(pvals<0.02)],i]=1
	newXVals[index[which(pvals<0.02)],i]=pvals[which(pvals<0.02)];
}

f=list(lst,newX, newXVals)
return(f)
}