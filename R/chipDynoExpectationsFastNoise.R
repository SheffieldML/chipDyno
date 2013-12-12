# CHIPDYNOEXPECTATIONSFASTNOISE computes posterior expectations TFA
# considering uncertainty of the expression level
# CHIPDYNO toolbox
# chipDynoExpectationsFastNoise.m version 1.4
# FORMAT chipDynoExpectationsFastNoise <- function(data, X, Sigma, beta, 
#			precs, gamma,mu, transNames, annotations, transName, geneName)
# DESC computes posterior expectations of TFA considering uncertainty of the expression level.
# ARG data : point estimate of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG Sigma : prior covariance matrix of TFA
# ARG beta :
# ARG precs : uncertainty of the expression level
# ARG gamma : degree of temporal continuity
# ARG mu : mean value of the transcription factor activity
# ARG TransNames : Transcription factors
# ARG annotations : Gene names
# ARG transName : specific transcription factor
# ARG geneName : specific gene name
# RETURN expectations : concateneted dataframe of gene specific transcription 
# factor activity, error in gene specific transcription factor activity 
# and corresponding difference in error
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoExpectationsFast

chipDynoExpectationsFastNoise <- function(data,X,Sigma,beta, precs, gamma,mu, transNames, 
											annotations, transName, geneName){

## Only for test purpose
# annotations = annotation
# transNames=TransNames
# transName= "ACE2"
# transName=activeNames[i]
# geneName="YHR143W"
####

npts=ncol(data);
nTrans=ncol(X);
c = class(geneName)

v= mat.or.vec(length(annotations),1)

if (c == 'character'){
	for (i in 1: nrow(annotations)) { 
		v[i] <- geneName==annotations[i,1]
		}
	x=data.matrix(X[which(v==1),])
	data=t(data[which(v==1),])
	precs = precs[which(v==1),]

} else if (c=='integer'){
	x=data.matrix(X[geneName,])
	data=data[geneName,]
	precs=precs[geneName,] ### t() ???

} else {
	print('Error: Genes can be identified either by number or name')
}

#source("chipDynoStatPostEstNoise.R")
expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);

expectations.b=expectations[[1]]
expectations.tfError=expectations[[2]]
expectations.tfErrorDiffs=expectations[[3]]

index=which(transName==transNames);
if (x[index,]== 0) {
 print('Error: The gene selected is not a target of the transcription factor')
}

tf=expectations.b[,index];
ind=which(transName == transNames[which(x!=0)]);
tfErrors=expectations.tfError[,ind];
tfErrorsDiffs=expectations.tfErrorDiffs[,,ind];

expectations = list(tf,tfErrors,tfErrorsDiffs);

return(expectations)
}
