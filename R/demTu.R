# DEMTU demonstrates dynamical chipCHIP on Tu data.
#
#	Description:
#	
## 	rdemTu.R version 0.01
#	Written on 14.02.2011


#####
#save.image("filename.RData")
# load("filename.RData")
####

rm(list = ls()) # Clear the workspace

# read the value of [data,vars,X,annotation,TransNames]

library(Matrix);

source("chipDynoTuLoadData.R");
data_vars_X_annotation_TransNames = chipDynoTuLoadData();

data = data_vars_X_annotation_TransNames[[1]] 		# Experimental data
vars = data_vars_X_annotation_TransNames[[2]]		# Standard Error
X = data_vars_X_annotation_TransNames[[3]]		# Binary connectivity matrix
annotation = data_vars_X_annotation_TransNames[[4]]	# Gene Names
TransNames = data_vars_X_annotation_TransNames[[5]]	# Transcripton factors

nGenes= nrow(data)
npts= ncol(data)
nTrans = ncol(X)
#options <- c(0,1e-4,1e-4,1e-6,0,0,0,0,0,0,0,0,0,0,0,1e-8,0.1,0) # Will be obsolute
## Source : http://www.engr.colostate.edu/~echong/ece520/matlab_demos/foptions.m

#options(1)=1;
#options(14)=50000;
#muIn=zeros(nTrans,1);

#options[1]=1
#options[14]=50000

muIn = array(0, dim <-c(nTrans,1)); 

#[R,C,V,nEffectGenes]=chipReduceVariables(X);
source('chipReduceVariables.R') 
R_C_V_nEffectGenes = chipReduceVariables(X);
R = R_C_V_nEffectGenes[[1]]
C = R_C_V_nEffectGenes[[2]]
V = R_C_V_nEffectGenes[[3]]
nEffectGenes = R_C_V_nEffectGenes[[4]]

diagonal = array(0.5, dim <- c(1,nTrans))
		
precs_mat = array(1, dim <- c(nrow(vars),ncol(vars)))
precs = precs_mat/(vars^2)

beta=3;
gamma=pi/4;
params=matrix(c(beta,gamma,t(muIn),0.1*t(V), diagonal),1,)


###
options = array(0, dim <- c(1,18))
options[1]=1;
options[2]=0.0001
options[3]=0.0001
options[14]=10000 # No of iteration
options[17]=0.1
##

source("chipDynoLikeStatNoise.R")
source("chipDynoLikeStatNoiseGrad.R")
source("SCGoptimTU.R")

params = SCGoptimTU(params, options, data, precs, X, nEffectGenes, R, C)


V=params[(3+nTrans):(length(params)-nTrans)]
#preSigma <- sparseMatrix(R, C, x=V, dims = c(nEffectGenes,nTrans))
preSigma <- as.matrix(sparseMatrix(R, C, x=V, dims = c(nEffectGenes,nTrans)))
diagonal = params[(length(params)-nTrans+1) :(length(params))];
Sigma = t(preSigma)%*%preSigma + diag(diagonal*diagonal)

beta = params[1]
gamma = params[2]
mu = params[3:(2+nTrans)]

# plot(mu,type="l",col="#22AAC6")
# plot(mu,type="l",col="red")

save.image("ResultsTu_New_10000Ita_Lovelace.RData")


# V=params(nTrans+3:end-nTrans)';
# preSigma=sparse(R,C,V,nEffectGenes,nTrans);
# diagonal=params(end-nTrans+1:end);
# Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
