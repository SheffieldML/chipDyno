# from chipReduceVariables.m
# function [R,C,V,nEffectGenes]=chipReduceVariables(X);

# CHIPREDUCEVARIABLES reduce  number of variables in chipDyno model
#
#	Description:
#	[R,C,V,nEffectGenes]=chipReduceVariables(X);
## 	chipReduceVariables.R version 0.1.0

chipReduceVariables = function(X){

#nGenes=size(X,1);
#preSigma1=X(1,:);
#preSigma2=preSigma1;

nGenes=nrow(X)
preSigma1=X[1,]
preSigma1=matrix(preSigma1,1,)
preSigma2=preSigma1

#for i=2:nGenes
#
#  preSigma2=[preSigma2;X(i,:)];
#  matrix1=preSigma1'*preSigma1;
#  matrix2=preSigma2'*preSigma2;         
#  decider=min(matrix1(find(matrix2)));
#  if decider==0
#    preSigma1=[preSigma1;X(i,:)];
#  end
#end

for (i in 2: nGenes){
	preSigma2 = rbind(preSigma2,X[i,])
	matrix1=t(preSigma1)%*%preSigma1
	matrix2=t(preSigma2)%*%preSigma2
	decider=min(matrix1[which(matrix2 !=0)])
	if (decider==0)
		preSigma1=rbind(preSigma1,X[i,])
}


#nEffectGenes=size(preSigma1,1);
#[R,C,V]=find(preSigma1);

nEffectGenes=nrow(preSigma1)
R = row(preSigma1)[which(preSigma1 != 0)]
C = col(preSigma1)[which(preSigma1 != 0)]
V = preSigma1[which(preSigma1 != 0)]

val= list(R,C,V,nEffectGenes)
return (val)
}
