#function [list,maxActivity,maxActivityError]=chipDynoGeneAct(data, ...
#                                                  X,Sigma,beta,gamma,mu, ...
#                                                  transNames, ...
#                                                  annotation,geneName);
#
#% CHIPDYNOGENEACT given a gene, lists activators in decreasing order
#%
#%	Description:
#%	[list,maxActivity,maxActivityError]=chipDynoGeneAct(data, ...
#%                                                  X,Sigma,beta,gamma,mu, ...
#%                                                  transNames, ...
#%                                                  annotation,geneName);
#%% 	chipDynoGeneAct.R version 0.1.0
#


chipDynoGeneAct = function(data, X, Sigma,beta,gamma,mu, transNames, annotation, geneName) {


## Only for test purpose
# annotations = annotation
# transNames=TransNames
# geneName="YHR143W"
# geneName="AGA1"
# geneName = "T20B12.8" # from Wormnet
####


I=which(geneName==annotation);
activeNames=transNames[which(X[I,]!=0)];
nTransFact=sum(X[I,]);
maxActivity=list();
maxActivityError=list();

#I=find(strcmp(geneName,annotation));
#activeNames=transNames(find(X(I,:)));
#nTransFact=sum(X(I,:));
#maxActivity=[];
#maxActivityError=[];

source("chipDynoExpectationsFast.R")

for (i in 1: nTransFact) {
	expectations=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, transNames, annotation,  activeNames[i], geneName);
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



#for i=1:nTransFact
#  [tf,tfError]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
#                                         transNames, annotation, ...
#                                         activeNames(i),geneName);
#  ind=find(strcmp(activeNames(i),transNames));
#  tf=tf-mu(ind)*ones(size(tf));
#  [act,index]=max(tf);
#  maxActivity=[maxActivity,act];
#  maxActivityError=[maxActivityError,tfError(index)];
#end


temp_ind = sort(unlist(maxActivity), decreasing = TRUE, index.return = TRUE)
maxActivity = temp_ind$x # '$x' return the sorted value 
index = temp_ind$ix	# '$ix' return the sorted value's index
maxActivityError=unlist(maxActivityError[index]);
list=activeNames[index];

#[maxActivity,index]=sort(maxActivity,'descend');
#maxActivityError=maxActivityError(index);
#list=activeNames(index);

name_value_error = list(list,maxActivity,maxActivityError)
return(name_value_error)

}
