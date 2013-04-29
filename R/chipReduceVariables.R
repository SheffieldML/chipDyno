# CHIPREDUCEVARIABLES reduce number of variables in chipDyno model
# CHIPDYNO toolbox
# chipReduceVariables.R version 1.0.1
# FORMAT chipReduceVariables <- function(X)
# DESC reduce  number of variables in chipDyno model
# ARG X : connectivity measurement between genes and transcription factors
# RETURN val : concatenated datafrane of same length integer vectors 
# specifying the row indices, column indices, the values of the non-zero 
# entries of the sparce matrix, and the number of effective genes
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : 

chipReduceVariables <- function(X){

nGenes=nrow(X)
preSigma1=X[1,]
preSigma1=matrix(preSigma1,1,)
preSigma2=preSigma1

for (i in 2: nGenes){
	preSigma2 = rbind(preSigma2,X[i,])
	matrix1=t(preSigma1)%*%preSigma1
	matrix2=t(preSigma2)%*%preSigma2
	decider=min(matrix1[which(matrix2 !=0)])
	if (decider==0)
		preSigma1=rbind(preSigma1,X[i,])
}

nEffectGenes=nrow(preSigma1)
R = row(preSigma1)[which(preSigma1 != 0)]
C = col(preSigma1)[which(preSigma1 != 0)]
V = preSigma1[which(preSigma1 != 0)]

val= list(R,C,V,nEffectGenes)
return (val)
}
