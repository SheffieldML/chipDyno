# CHIPDYNOGENEACTNOISE given a gene, lists activators in decreasing order
# CHIPDYNO toolbox
# chipDynoGeneActNoise.R version 1.0.1
# FORMAT chipDynoGeneAct <- function(data, X, Sigma,beta,precs, gamma,mu, 
#					transNames, annotation, geneName);
# DESC given a gene, lists activators in decreasing order
# ARG data : point estimate of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG Sigma : prior covariance matrix of TFA
# ARG beta :
# ARG precs :
# ARG gamma : degree of temporal continuity
# ARG mu : mean value of the transcription factor activity
# ARG transNames : Transcription factors
# ARG annotation : Gene names
# ARG geneName : specific gene name
# RETURN list : concatenated data frame of lists activators for a given gene 
# in decreasing order, maximum activity and maximum activity error
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoGeneAct, chipDynoActTransFact, chipDynoActTransFactNoise

chipDynoGeneActNoise = function(data, X, Sigma, beta, precs, gamma,mu, transNames, 
	annotation, geneName) {

## Only for test purpose
# source("test_chipDynoGeneActNoise.R")
# load("/home/muhammad/H-drive/CElegans/Results_cElegans_Optim_Sample1.RData")
# annotations = annotation
# transNames=TransNames
# geneName="YHR143W"
# geneName="AGA1"
# geneName = "W02D7.4" # from Wormnet
# geneName = annotation[i]
####

require("Matrix")
I=which(geneName==annotation);
activeNames=transNames[which(X[I,]!=0)];
nTransFact=sum(X[I,]);
maxActivity=list();
maxActivityError=list();

#source("chipDynoExpectationsFastNoise.R")

for (i in 1: nTransFact) {
	expectations=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, 
		transNames, annotation,  activeNames[i], geneName);
	tf = expectations[[1]]
	tfError = expectations[[2]]
	tfErrorDiffs = expectations[[3]]

	ind=which(activeNames[i]==transNames);
	tf = tf- mu[ind]*matrix(1, length(tf), 1);
	act=max(tf);
	index=which.max(tf)
	maxActivity=cbind(maxActivity,act);
	maxActivityError=cbind(maxActivityError,tfError[index]);
}

temp_ind = sort(unlist(maxActivity), decreasing = TRUE, index.return = TRUE)
maxActivity = temp_ind$x # '$x' return the sorted value 
index = temp_ind$ix	# '$ix' return the sorted value's index
maxActivityError=unlist(maxActivityError[index]);
list=activeNames[index];

name_value_error = list(list,maxActivity,maxActivityError)
return(name_value_error)
}
