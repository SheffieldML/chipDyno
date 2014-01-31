# DEMSPELLMANDYNOSTAT demonstrates dynamical chipCHIP on Spellman data.
# CHIPDYNO toolbox
# demSpellmanDynoStat.R version 1.0.1
# FORMAT 
# DESC demonstrates dynamical chipCHIP on Spellman data.
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : demTu

rm(list=ls())

cat("demSpellmanDynoStat.R demonstrates dynamical chipCHIP on Spellman data. \n");
cat("Loading and pre-prosessing data files ...\n ");

source("chipDynoLoadData.R")
data_X_annotation_TransNames=chipDynoLoadData()
data = data_X_annotation_TransNames[[1]]
X = data_X_annotation_TransNames [[2]]
annotation = data_X_annotation_TransNames[[3]]
TransNames = data_X_annotation_TransNames[[4]]
annotations = annotation # Both of the variavle contain the same data
transNames = TransNames # Both of the variavle contain the same data

nGenes=nrow(data);
npts=ncol(data);
nTrans=ncol(X);

g=cov(t(data))
muIn=array(0, dim <- c(nTrans,1));

options = array(0, dim=c(1,18))
options[1]=1;
options[2]=0.0001
options[3]=0.0001
options[14]=10000 # No of iteration
options[17]=0.1

cat("Creating a sparse matrix for gene vs TF connectivity  ...\n");

source('chipReduceVariables.R') 
R_C_V_nEffectGenes = chipReduceVariables(X);
R = R_C_V_nEffectGenes[[1]]
C = R_C_V_nEffectGenes[[2]]
V = R_C_V_nEffectGenes[[3]]
nEffectGenes = R_C_V_nEffectGenes[[4]]

diagonal=(t(rnorm(nTrans)))^2;
beta=3;
gamma=pi/4;

cat("Optimizing parameters... \n");

params=matrix(c(beta,gamma,t(muIn),rnorm(length(V))*V, diagonal),1,);

source("chipDynoLikeStat.R")
source("chipDynoLikeStatGrad.R")
source("SCGoptim.R")

params = SCGoptim(params, options, data, X, nEffectGenes, R, C)

library(Matrix)
V=params[(3+nTrans):(length(params)-nTrans)]
preSigma <- as.matrix(sparseMatrix(R, C, x=V, dims = c(nEffectGenes,nTrans)))
#preSigma <- sparseMatrix(R, C, x=V, dims = c(nEffectGenes,nTrans))
diagonal = params[(length(params)-nTrans+1) :(length(params))];
Sigma = t(preSigma)%*%preSigma + diag(diagonal*diagonal)
beta = params[1]
gamma = params[2]
mu = params[3:(2+nTrans)]

save.image("ResultsSpellman_ita_10000.RData")