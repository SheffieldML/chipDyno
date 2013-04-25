# CHIPDYNOEXPECTATIONSFAST computes posterior expectations of TFA.
# CHIPDYNO toolbox
# chipDynoExpectationsFast.R version 1.1
# FORMAT chipDynoExpectationsFast<-function(data,X,Sigma,beta,gamma,mu,
#                                      transNames, annotations, transName,geneName);
# DESC computes posterior expectations of TFA.
# ARG data : point estimate of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG Sigma : prior covariance matrix of TFA
# ARG beta :
# ARG gamma : degree of temporal continuity
# ARG mu : mean value of the transcription factor activity
# ARG transNames : Transcription factors
# ARG annotations : Gene names
# ARG transName : specific transcription factor
# ARG geneName : specific gene name
# RETURN expectations : concateneted dataframe of gene specific transcription 
# factor activity, error in gene specific transcription factor activity 
# and corresponding difference in error
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoExpectationsFastNoise


chipDynoExpectationsFast<-function(data,X,Sigma,beta,gamma,mu, transNames, annotations, transName, geneName){

## Only for test purpose
#annotations = annotation
#transNames=TransNames
#transName= name
#geneName= genesIn[i]
#transName= "ACE2" "FKH1" "FKH2" "MBP1"
#geneName="YHR143W"
####

npts=ncol(data);
nTrans=ncol(X);
c = class(geneName)

v= mat.or.vec(length(annotations),1)

if (c == 'character'){
	for (i in 1: length(annotations)) { 
		v[i] <- geneName==annotations[i]
		}
	x=data.matrix(X[which(v==1),])
	data=t(data[which(v==1),])

} else if (c=='integer'){
	x=data.matrix(X[geneName,])
	data=data[geneName,]	

} else {
	print('Error: Genes can be identified either by number or name')
}

source("chipDynoStatPostEst.R")
expectations= chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);

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

return (expectations)

}