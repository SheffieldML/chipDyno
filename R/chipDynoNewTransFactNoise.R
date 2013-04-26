# CHIPDYNONEWTRANSFACTNOISE transcription factors active for us and not for Tu et al.
# CHIPDYNO toolbox
# chipDynoNewTransFactNoise.R version 1.0.1
# FORMAT list=chipDynoNewTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#                                         TransNames, annotation);
# DESC determins tfs active for us and not for Tu et al
# ARG data : point estimate of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG Sigma : prior covariance matrix of TFA
# ARG beta :
# ARG precs : uncertainty of the expression level
# ARG gamma : degree of temporal continuity
# ARG mu : mean value of the transcription factor activity
# ARG TransNames : Transcription factors
# ARG annotations : Gene names
# RETURN list: transcription factors active for us and not for Tu et al
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoTransFactNoise, chipDynoNewTransFactNoise

chipDynoNewTransFactNoise<-function(data, X, Sigma, beta, precs,  
									gamma, mu, TransNames, annotation){

file= "./data/MetabolData/PerTransFact.txt"

TuTransFact <- read.table(file, header = TRUE, sep = "\t", quote = "", dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

source("chipDynoActTransFactNoise.R")
list1_newX_newXVals=chipDynoActTransFactNoise(data,X,Sigma,beta,precs,gamma,mu,TransNames, annotation);

list1 = list1_newX_newXVals[[1]]
newX = list1_newX_newXVals[[2]]
newXVals = list1_newX_newXVals[[3]]

list2= TransNames[which(list1>4)]
index=list();

for (i in 1: nrow(list2)){
#### TODO
}

#for i=1:size(list2,1)
#  index=[index,1-size(find(strcmp(list2(i),TuTransFact)),1)];
#end
#list=list2(find(index));
#  

return(list)
}