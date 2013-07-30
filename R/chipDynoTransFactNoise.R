# CHIPDYNOTRANSFACTNOISE given a transcription factor, provides TFAs with errorbars .
# CHIPDYNO toolbox
# chipDynoTransFactNoise.R version 1.0.1
# FORMAT chipDynoTransFactNoise <- function(data, X, Sigma, beta, precs, 
#									gamma, mu, transNames, annotations, name)
# DESC given a transcription factor, provides gene-specific TFAs with errorbars.
# ARG data : point estimate of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG Sigma : prior covariance matrix of TFA
# ARG beta :
# ARG precs : uncertainty of the expression level
# ARG gamma : degree of temporal continuity
# ARG mu : mean value of the transcription factor activity
# ARG transNames : Transcription factors
# ARG annotations : Gene names
# ARG name : given transcription factor name
# RETURN expectations : concatenated dataframe of transcription factor 
# activity and its error
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoTransFact

chipDynoTransFactNoise <- function(data, X, Sigma, beta, precs, 
								gamma, mu, transNames, annotations, name) {

###
# Only for test purpose
# name= "ZAP1"
# transNames = TransNames
# annotations = annotation
# i=1
##

index=which(name == transNames)
genesIn=which(X[,index]!=0)
anno=annotations[which(X[,index]!=0)]
nTargets=length(anno)
npts=ncol(data)
TF=array(0, dim=c(nTargets,npts))
TFError=array(0, dim=c(nTargets,npts))
TFErrorDiff=array(0, dim=c(npts,npts,nTargets))

source("chipDynoExpectationsFastNoise.R")

for (i in 1 : nTargets) {
	expectations = chipDynoExpectationsFastNoise(data,X,Sigma,beta, precs,
						gamma,mu, transNames, annotations, name, genesIn[i])

	TF[i,] = expectations[[1]]
	TFError[i,] = expectations[[2]]
	TFErrorDiff[ , ,i] = expectations[[3]]
}

expectations = list(TF,TFError,TFErrorDiff)
return(expectations)
}