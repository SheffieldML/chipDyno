# CHIPDYNOLIKESTAT marginal likelihood for chipChip dynamical model
# CHIPDYNO toolbox
# chipDynoLikeStat.R version 1.0.1
# FORMAT chipDynoLikeStat <- function(params,data, X,nEffectGenes,R,C)
# DESC compute the marginal likelihood for chipChip dynamical model
# ARG params: concatenated vector of multiple parameters(beta, gamma, 
# initial mean of the transcription factors, and 
# a vector to create diagonal matrix used to reduce the sparsity of covariance)
# ARG data : point estimate of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG nEffectGenes : effectice gene name
# ARG R, C : same length integer vectors specifying the row and column 
# indices of the non-zero entries of the sparce matrix
# RETURN likelihood : marginal likelihood
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoLikeStatNoise

chipDynoLikeStat <- function(params,data, X,nEffectGenes,R,C){

nGenes=nrow(data)
nTrans=ncol(X)
npts=ncol(data)
beta=params[1]
gamma=params[2]

mu=t(matrix(params[,3:(2+nTrans)],1,))
V=t(matrix(params[,(3+nTrans):(ncol(params)-nTrans)],1,))
V=as.vector(V)
library(Matrix)
preSigma <- as.matrix(sparseMatrix(R, C, x=V, dims = c(nEffectGenes,nTrans)))

diagonal = as.vector(params[(ncol(params)-nTrans+1):ncol(params)])

Sigma=t(preSigma)%*%preSigma+diag(diagonal*diagonal); 
factor=cos(gamma)^2
preLike= mat.or.vec(1, nGenes)
k=mat.or.vec(1, npts)
alpha=mat.or.vec(1, npts)

for (i in 1: nGenes){
	coeff=(1-factor^2)*X[i,] %*% Sigma %*% X[i,];
	number=X[i,]%*%mu;
	alpha[1]=beta^2;
	k[1]=(data[i,npts]-(1-factor)*number)/factor;
	sigma=beta^-2+(alpha[1]^-1+coeff)/factor^2;
	#update=0.5*log(sigma)+ 0.5*((k[1]-data[i,npts-1])^2)/(sigma); 
	preLike[i]=preLike[i]+ (0.5*log(sigma)+ 0.5*((k[1]-data[i,npts-1])^2)/(sigma))	

	for (j in 2: (npts-1)){
    		alpha[j]= beta^2+factor^2*(alpha[j-1]^-1+coeff)^-1;
		k[j]= alpha[j]^-1*(data[i,npts-j+1]*
			beta^2+factor^2*k[j-1]*(alpha[j-1]^-1+coeff)^-1);
		sigma=beta^-2+(alpha[j]^-1+coeff)/factor^2;
		preLike[i]= preLike[i]+ 0.5*log(sigma)+
			0.5*((k[j]-data[i,npts-j])^2)/(sigma);
	}
	
	alpha[length(alpha)]= beta^2+factor^2*(alpha[length(alpha)-1]^-1+coeff)^-1;
	k[length(k)]= alpha[length(alpha)]^-1*(data[i,1]*beta^2+
		factor^2*k[length(k)-1]*(alpha[length(alpha)-1]^-1+coeff)^-1);
	preLike[i]=preLike[i]+ (0.5*log(alpha[length(alpha)]^-1+X[i,]%*%Sigma%*%X[i,])+
		0.5*((k[length(k)]-number)^2)/(alpha[length(alpha)]^-1+X[i,]%*%Sigma%*%X[i,]));
}

likelihood=rowSums(preLike);

return(likelihood)
}