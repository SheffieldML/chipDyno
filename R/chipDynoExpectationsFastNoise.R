#function [tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#                                         TransNames, annotation, ...
#                                         transName,geneName);
#
#% CHIPDYNOEXPECTATIONSFASTNOISE computes posterior expectations TFA.
#%
#%	Description:
#%	[tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
#%                                         TransNames, annotation, ...
#%                                         transName,geneName);
#%% 	chipDynoExpectationsFastNoise.m version 0.1.0

chipDynoExpectationsFastNoise=function(data,X,Sigma,beta, precs, gamma,mu, transNames, annotations, transName, geneName){


## Only for test purpose
#annotations = annotation
#transNames=TransNames
#transName= "ACE2"
#geneName="YHR143W"
####


npts=ncol(data);
nTrans=ncol(X);
c = class(geneName)

#npts=size(data,2);
#nTrans=size(X,2);
#c=class(geneName);

#v= mat.or.vec(nrow(annotations),1)
v= mat.or.vec(length(annotations),1)

if (c == 'character'){
	for (i in 1: nrow(annotations)) { 
		v[i] <- geneName==annotations[i,1]
		}
	x=data.matrix(X[which(v==1),])
	#x=t(data.matrix(X[which(v==1),]))
	data=t(data[which(v==1),])
	precs = precs[which(v==1),]

} else if (c=='integer'){
	x=data.matrix(X[geneName,])
	#x=t(data.matrix(X[geneName,]))
	#data=t(data[geneName,])
	data=data[geneName,]
	precs=precs[geneName,] ### t() ???

} else {
	print('Error: Genes can be identified either by number or name')
}

#if c(1)=='c'
#    x=X(find(strcmp(geneName,annotation)),:)';
#    data=data(find(strcmp(geneName,annotation)),:);
#    precs=precs(find(strcmp(geneName,annotation)),:);
#elseif c(1)=='d'
#    x=X(geneName,:)';
#    data=data(geneName,:);
#    precs=precs(geneName,:);
#else
#    error('Genes can be identified either by number or name\n')
#end


source("chipDynoStatPostEstNoise.R")
expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);

expectations.b=expectations[[1]] # First segment of the returned list
expectations.tfError=expectations[[2]] # Second segment of the returned list
expectations.tfErrorDiffs=expectations[[3]] # Third segment of the returned list

index=which(transName==transNames);
if (x[index,]== 0) {
 print('Error: The gene selected is not a target of the transcription factor')
}


#expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);
#index=find(strcmp(transName,TransNames));
#if x(index)==0
#  error('The gene selected is not a target of the transcription factor \n')
#end

tf=expectations.b[,index];
ind=which(transName == transNames[which(x!=0)]);
tfErrors=expectations.tfError[,ind];
tfErrorsDiffs=expectations.tfErrorDiffs[,,ind];

#tf=expectations.b(:,index);
#ind=find(strcmp(transName,TransNames(find(x))));
#tfErrors=expectations.tfError(:,ind);
#tfErrorsDiffs=[expectations.tfErrorDiffs(:,:,ind)];

expectations = list(tf,tfErrors,tfErrorsDiffs);

return(expectations)

}
