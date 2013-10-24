require("Matrix")

load("/home/muhammad/Dropbox/CElegans/cluster/3exp_15dp_After_meeting/Results_cElegans_cluster3_Optim_CE_CX_Sample115.RData")
#load("/home/muhammad/Dropbox/CElegans/cluster/3exp_15dp_After_meeting/Results_cElegans_cluster4_Optim_CE_CX_Sample115.RData")

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));

noRegulatedGenes = mat.or.vec(1, length(TransNames))

for (i in 1 : length(TransNames)){

expectations =chipDynoTransFactNoise(data,X,Sigma,beta, precs, gamma,mu, TransNames, annotation, TransNames[i,]);
S1_TF = expectations[[1]]
S1_TFError = expectations [[2]]
S1_TFErrorDiff = expectations [[3]]
noRegulatedGenes[i] = nrow(S1_TF)
}

TFnoRegulatedGenes = cbind(TransNames,t(noRegulatedGenes))