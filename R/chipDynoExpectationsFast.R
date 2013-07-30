#% CHIPDYNOEXPECTATIONSFAST computes posterior expectations of TFA.
#%
#%	Description:
#%	[tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
#%                                         transNames, annotations, ...
#%                                         transName,geneName);
#%% 	chipDynoExpectationsFast.R version 0.1.0


chipDynoExpectationsFast=function(data,X,Sigma,beta,gamma,mu, transNames, annotations, transName, geneName){

#function [tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
#                                         transNames, annotations, ...
#                                         transName,geneName);

#% CHIPDYNOEXPECTATIONSFAST computes posterior expectations of TFA.
#%
#%	Description:
#%	[tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
#%                                         transNames, annotations, ...
#%                                         transName,geneName);
#%% 	chipDynoExpectationsFast.m version 1.5


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

#npts=size(data,2);
#nTrans=size(X,2);
#c=class(geneName);


#v= mat.or.vec(nrow(annotations),1)
v= mat.or.vec(length(annotations),1)

if (c == 'character'){
	for (i in 1: length(annotations)) { 
		v[i] <- geneName==annotations[i]
		}
	x=data.matrix(X[which(v==1),])
	#x=t(data.matrix(X[which(v==1),]))
	data=t(data[which(v==1),])

} else if (c=='integer'){
	x=data.matrix(X[geneName,])
	#x=t(data.matrix(X[geneName,]))
	#data=t(data[geneName,])
	data=data[geneName,]	

} else {
	print('Error: Genes can be identified either by number or name')
}


#if c(1)=='c'
#    x=X(find(strcmp(geneName,annotations)),:)';
#    data=data(find(strcmp(geneName,annotations)),:);
#elseif c(1)=='d'
#    x=X(geneName,:)';
#    data=data(geneName,:);
#else
#    error('Genes can be identified either by number or name\n')
#end


source("chipDynoStatPostEst.R")
expectations= chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);

expectations.b=expectations[[1]] # First segment of the returned list
expectations.tfError=expectations[[2]] # Second segment of the returned list
expectations.tfErrorDiffs=expectations[[3]] # Third segment of the returned list

index=which(transName==transNames);

if (x[index,]== 0) {
 print('Error: The gene selected is not a target of the transcription factor')
}


tf=expectations.b[,index];
ind=which(transName == transNames[which(x!=0)]);
tfErrors=expectations.tfError[,ind];
tfErrorsDiffs=expectations.tfErrorDiffs[,,ind];


#expectations=chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);
#index=find(strcmp(transName,transNames));
#if x(index)==0
#  error('The gene selected is not a target of the transcription factor \n')
#end

#tf=expectations.b(:,index);
#ind=find(strcmp(transName,transNames(find(x))));
#tfErrors=expectations.tfError(:,ind);
#tfErrorsDiffs=[expectations.tfErrorDiffs(:,:,ind)];

expectations = list(tf,tfErrors,tfErrorsDiffs);

return (expectations)

}
