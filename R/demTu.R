# DEMTU demonstrates dynamical chipCHIP on Tu data.
# CHIPDYNO toolbox
# demTu.R version 1.0.1
# FORMAT 
# DESC demonstrates dynamical chipCHIP on Tu data.
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : demSpellmanDynoStat

rm(list = ls()) 

library(Matrix);


cat("demTu.R demonstrates dynamical chipCHIP on Tu data. \n");
cat("Loading and pre-prosessing data files ...\n ");

file_dictionary <- "./data/MetabolData/dictionary.txt";
file_probeIDTu <- "./data/MetabolData/probeIDTu.txt";
file_data <- "./data/MetabolData/YeastMetabolism_exprs.txt";
file_vars <- "./data/MetabolData/YeastMetabolism_se.txt";

file_dataChip <- "./data/Connectivity2.txt";
file_annotation <- "./data/annotations2.txt"
file_transNames <- "./data/Trans_Names2.txt"

source("chipDynoTuLoadData.R");
data_vars_X_annotation_TransNames = chipDynoTuLoadData(file_dictionary, 
	file_probeIDTu, file_data, file_vars, file_dataChip, 
	file_annotation, file_transNames);

data = data_vars_X_annotation_TransNames[[1]]
vars = data_vars_X_annotation_TransNames[[2]]
X = data_vars_X_annotation_TransNames[[3]]
annotation = data_vars_X_annotation_TransNames[[4]]
TransNames = data_vars_X_annotation_TransNames[[5]]
annotations = annotation # Both of the variavle contain the same data
transNames = TransNames # Both of the variavle contain the same data

nGenes= nrow(data)
npts= ncol(data)
nTrans = ncol(X)
muIn = array(0, dim <-c(nTrans,1)); 

cat("Creating a sparse matrix for gene vs TF connectivity  ...\n");

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

options = array(0, dim <- c(1,18))
options[1]=1;
options[2]=0.0001
options[3]=0.0001
options[14]= 5  # No of iteration
options[17]=0.1

source("chipDynoLikeStatNoise.R")
source("chipDynoLikeStatNoiseGrad.R")
source("SCGoptimTU.R")

cat("Optimizing parameters... \n");

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

save.image("ResultsTu.RData")