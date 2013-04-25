#function list=chipDynoNewTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#                                         TransNames, annotation);
#
#% CHIPDYNONEWTRANSFACTNOISE tfs active for us and not for Tu et al.
#%
#%	Description:
#%	list=chipDynoNewTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#%                                         TransNames, annotation);
#%% 	chipDynoNewTransFactNoise.R version 0.1.0

chipDynoNewTransFactNoise = function(data, X, Sigma, beta, precs, gamma, mu, TransNames, annotation){


file= "./data/MetabolData/PerTransFact.txt"

TuTransFact <- read.table(file, header = TRUE, sep = "\t", quote = "", dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

#TuTransFact=textread('./data/MetabolData/PerTransFact.txt', ...
#                     '%q');

source("chipDynoActTransFactNoise.R")
list1_newX_newXVals=chipDynoActTransFactNoise(data,X,Sigma,beta,precs,gamma,mu,TransNames, annotation);

list1 = list1_newX_newXVals[[1]]
newX = list1_newX_newXVals[[2]]
newXVals = list1_newX_newXVals[[3]]

list2= TransNames[which(list1>4)]
index=list();
#list2=TransNames(find(list1>4));
#index=[];


for (i in 1: nrow(list2)){

}

#for i=1:size(list2,1)
#  index=[index,1-size(find(strcmp(list2(i),TuTransFact)),1)];
#end
#list=list2(find(index));
#  

return(list)
}
