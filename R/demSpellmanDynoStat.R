#% DEMSPELLMANDYNOSTAT demonstrates dynamical chipCHIP on Spellman data.
#%
#%	Description:
#%	% 	demSpellmanDynoStat.R version 0.1.0
#rm(list=ls())
#clear all
#randn('seed',39)

source("chipDynoLoadData.R")
data_X_annotation_TransNames=chipDynoLoadData()
data = data_X_annotation_TransNames[[1]]
X = data_X_annotation_TransNames [[2]]
annotation = data_X_annotation_TransNames[[3]]
TransNames = data_X_annotation_TransNames[[4]]


#[data,X,annotation,TransNames]=chipDynoLoadData();

nGenes=nrow(data);
npts=ncol(data);
nTrans=ncol(X);

#options=foptions;
#options(1)=1;
#options(14)=10000;
g=cov(t(data))
muIn=array(0, dim <- c(nTrans,1));

#
options = array(0, dim=c(1,18))
options[1]=1;
options[2]=0.0001
options[3]=0.0001
options[14]=5 # No of iteration
options[17]=0.1
#

#nGenes=size(data,1);
#npts=size(data,2);
#nTrans=size(X,2);
#options=foptions;
#options(1)=1;
#options(14)=10000;
#g=cov(data');
#muIn=zeros(nTrans,1);

source('chipReduceVariables.R') 
R_C_V_nEffectGenes = chipReduceVariables(X);
R = R_C_V_nEffectGenes[[1]]
C = R_C_V_nEffectGenes[[2]]
V = R_C_V_nEffectGenes[[3]]
nEffectGenes = R_C_V_nEffectGenes[[4]]

#[R,C,V,nEffectGenes]=chipReduceVariables(X);

#### remove TODO only for test purpose
#diagonal=matrix(mat.or.vec(1,nTrans),1,) # create nX1 matrix of 0's
#diagonal[,]=0.5
######

diagonal=(t(rnorm(nTrans)))^2;
beta=3;
gamma=pi/4;


#diagonal=(randn(1,nTrans)).^2;
#beta=3;
#gamma=pi/4;

##remove TODO
#params=matrix(c(beta,gamma,t(muIn),t(V), diagonal),1,);
#####

params=matrix(c(beta,gamma,t(muIn),rnorm(length(V))*V, diagonal),1,);

source("chipDynoLikeStat.R")
source("chipDynoLikeStatGrad.R")
source("SCGoptimSpellDS.R")

params = SCGoptimSpellDS(params, options, data, X, nEffectGenes, R, C)

#params=[beta,gamma,muIn',randn(1,size(V,1)).*V', diagonal];
#params=scg('chipDynoLikeStat',params,options,'chipDynoLikeStatGrad',data, ...
#           X, nEffectGenes,R, C);

#save.image("ResultsSpellman_5000_temp.RData")

library(Matrix)
V=params[(3+nTrans):(length(params)-nTrans)]
preSigma <- as.matrix(sparseMatrix(R, C, x=V, dims = c(nEffectGenes,nTrans)))
#preSigma <- sparseMatrix(R, C, x=V, dims = c(nEffectGenes,nTrans))
diagonal = params[(length(params)-nTrans+1) :(length(params))];
Sigma = t(preSigma)%*%preSigma + diag(diagonal*diagonal)
beta = params[1]
gamma = params[2]
mu = params[3:(2+nTrans)]

# plot(mu,type="l",col="#22EEC6")
# plot(mu,type="l",col="red")

save.image("ResultsSpellman_5000Ita.RData")



#V=params(nTrans+3:end-nTrans)';
#preSigma=sparse(R,C,V,nEffectGenes,nTrans);
#diagonal=params(end-nTrans+1:end);
#Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
#beta=params(1);
#gamma=params(2);
#mu=params(3:2+nTrans);
#save ResultsSpellmanNewInit params
