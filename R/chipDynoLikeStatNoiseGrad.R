# CHIPDYNOLIKESTATNOISEGRAD gradient of CHIPDYNOLIKESTATNOISE.
# CHIPDYNO toolbox
# chipDynoLikeStatNoiseGrad.R version 1.0.1
# FORMAT chipDynoLikeStatNoiseGrad <- function(params, data, precs, X, nEffectGenes, R, C)
# DESC compute the gradient of chipDynoLikeStatNoise for chipChip dynamical model
# ARG params: concatenated vector of multiple parameters(beta, gamma, 
# initial mean of the transcription factors, and 
# a vector to create diagonal matrix used to reduce the sparsity of covariance)
# ARG data : point estimate of the expression level
# ARG precs : uncertainty of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG nEffectGenes : effectice gene name
# ARG R, C : same length integer vectors specifying the row and column 
# indices of the non-zero entries of the sparce matrix
# RETURN f : a list of beta,gamma, mu, Sigma, diagonal
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoLikeStatGrad

chipDynoLikeStatNoiseGrad <- function(params, data, precs, X, nEffectGenes, R, C) {

nGenes <- nrow(data)
nTrans <- ncol(X)
npts <- ncol(data)
beta <- params[1]
gamma <- params[2]

factor <- cos(gamma)^2
mu <- t(matrix(params[,3:(2+nTrans)],1,))
V <- t(matrix(params[,(3+nTrans):(ncol(params)-nTrans)],1,))
V <- as.vector(V)

library(Matrix)
preSigma <- sparseMatrix(R, C, x=V, dims = c(nEffectGenes,nTrans))
preSigma <- as.matrix(preSigma)
diagonal <- params[(ncol(params)-nTrans+1):ncol(params)]

Sigma <- t(preSigma)%*%preSigma+diag(diagonal*diagonal);
preGradBeta <- mat.or.vec(1, nGenes);
preGradSigmaScalar <- mat.or.vec(1, nGenes);
preGradMuScalar <- mat.or.vec(1, nGenes);
preGradGamma <- 0;
k <- mat.or.vec(1,npts);
alpha <- mat.or.vec(1,npts);

gradAlphaBeta <- mat.or.vec(1,npts);
gradAlphaGamma <- mat.or.vec(1,npts);
gradAlphaSigmaScalar <- mat.or.vec(1,npts);
gradKBeta <- mat.or.vec(1,npts);
gradKGamma <- mat.or.vec(1,npts);
gradKSigmaScalar <- mat.or.vec(1,npts);
gradKMuScalar <- mat.or.vec(1,npts);

for(i in 1:nGenes){
	number <- X[i,] %*% mu;
	coeff <- (1-factor^2)*(t(X[i,])%*% Sigma %*% X[i,])
	alpha[1] <- beta^2+precs[i,npts];
	k[1] <- (data[i,npts]-(1-factor)*number)/factor;
	sigma <- (beta^2+precs[i,(npts-1)])^-1+(alpha[1]^-1+coeff)/factor^2;
	gradientSigma <- 0.5*(sigma^-1-(k[1]-data[i,(npts-1)])^2/sigma^2);
	gradK <- (k[1]-data[i,(npts-1)])/sigma;

	gradAlphaBeta[1] <- 2*beta;
  	gradAlphaSigmaScalar[1] <- 0;
	gradAlphaGamma[1] <- 0;
	gradKBeta[1] <- 0;
	gradKSigmaScalar[1] <- 0;
	gradKGamma[1] <- -k[1]/factor+number/factor;
	gradKMuScalar[1] <- -(1-factor)/factor;
	preGradBeta[i] <- gradientSigma*(-2*(beta^2+precs[i,npts-1])^-2*beta-
		alpha[1]^-2*gradAlphaBeta[1]/factor^2);
	preGradGamma <- preGradGamma+gradK*gradKGamma[1]+gradientSigma*(-2*(alpha[1]^-1+coeff)/
		factor^3-2*coeff/(factor*(1-factor^2)));
	preGradSigmaScalar[i]<- gradientSigma*(1-factor^2)/factor^2;
	preGradMuScalar[i] <- gradK*gradKMuScalar[1];
	
	for(j in 2:(npts-1)){
		gradAlphaPrev <- factor^2*(alpha[j-1]^-1+coeff)^-2*alpha[j-1]^-2;  
    		gradAlphaBeta[j] <- 2*beta+gradAlphaPrev*gradAlphaBeta[j-1];
		gradAlphaSigmaScalar[j] <- -factor^2*(1-factor^2)*
			(alpha[j-1]^-1+coeff)^-2+gradAlphaPrev*gradAlphaSigmaScalar[j-1];
		gradAlphaGamma[j] <- 2*factor^3*coeff*(1-factor^2)^-1*
			(alpha[j-1]^-1+coeff)^-2+gradAlphaPrev*gradAlphaGamma[j-1]+
			2*factor*(alpha[j-1]^-1+coeff)^-1;
    		alpha[j] <- (beta^2+precs[i,(npts-j+1)]+factor^2*
			(alpha[j-1]^-1+coeff)^-1);
		k[j] <- alpha[j]^-1*(data[i,(npts-j+1)]*beta^2+factor^2*
			k[j-1]*(alpha[j-1]^-1+coeff)^-1);
		sigma <- (beta^2+precs[i,npts-j])^-1+(alpha[j]^-1+coeff)/factor^2; 

		gradientSigma <- 0.5*(sigma^-1-(k[j]-data[i,npts-j])^2/sigma^2);
    		gradK <- (k[j]-data[i,npts-j])/sigma;
    		gradKPrev <- alpha[j]^-1*factor^2*(alpha[j-1]^-1+coeff)^-1;
    		gradKAlpha <- -k[j]/alpha[j];
    		gradKAlphaPrev <- factor^2*alpha[j]^-1*k[j-1]*alpha[j-1]^-2*(alpha[j-1]^-1+coeff)^-2;

		gradKBeta[j] <- 2*alpha[j]^-1*data[i,(npts-j+1)]*beta+
			gradKPrev*gradKBeta[j-1]+gradKAlpha*gradAlphaBeta[j]+
			gradKAlphaPrev*gradAlphaBeta[j-1];

		gradKGamma[j] <- 2*factor*alpha[j]^-1*(alpha[j-1]^-1+coeff)^-1*k[j-1]*
			(1+factor^2*coeff*(alpha[j-1]^-1+coeff)^-1*(1-factor^2)^-1)+ 
			gradKAlpha*gradAlphaGamma[j]+
			gradKAlphaPrev*gradAlphaGamma[j-1]+
			gradKPrev*gradKGamma[j-1];

		gradKSigmaScalar[j] <- -factor^2*alpha[j]^-1*k[j-1]*(alpha[j-1]^-1+
			coeff)^-2*(1-factor^2)+gradKAlpha*gradAlphaSigmaScalar[j]+
			gradKAlphaPrev*gradAlphaSigmaScalar[j-1]+
			gradKPrev*gradKSigmaScalar[j-1];     

		gradKMuScalar[j] <- gradKPrev*gradKMuScalar[j-1];
		gradSigmaBeta <- -2*(beta^2+precs[i,(npts-j)])^-2*beta;
		gradSigmaAlpha <- -alpha[j]^-2*factor^-2;
		gradSigmaGamma <- -2*(alpha[j]^-1+coeff)*factor^-3-2*coeff*factor^-1*(1-factor^2)^-1;
		gradSigmaSigmaScalar <- (1-factor^2)*factor^-2;

		preGradBeta[i] <- preGradBeta[i]+gradientSigma*
			(gradSigmaBeta+gradSigmaAlpha*gradAlphaBeta[j])+gradK*gradKBeta[j];

		preGradSigmaScalar[i] <- preGradSigmaScalar[i]+gradientSigma*
			(gradSigmaSigmaScalar+gradSigmaAlpha*gradAlphaSigmaScalar[j])+
			gradK*gradKSigmaScalar[j];
		preGradGamma <- preGradGamma+gradientSigma*(gradSigmaGamma+
			gradSigmaAlpha*gradAlphaGamma[j])+gradK*gradKGamma[j];

		preGradMuScalar[i] <- preGradMuScalar[i]+gradK*gradKMuScalar[j];
	}

	gradAlphaPrev <- factor^2*(alpha[length(alpha)-1]^-1+coeff)^-2*alpha[length(alpha)-1]^-2;
  	gradAlphaBeta[length(gradAlphaBeta)] <- 2*beta+gradAlphaPrev*
		gradAlphaBeta[length(gradAlphaBeta)-1];
	gradAlphaSigmaScalar[length(gradAlphaSigmaScalar)] <- -factor^2*(1-factor^2)*
		(alpha[length(alpha)-1]^-1+coeff)^-2+gradAlphaPrev*
		gradAlphaSigmaScalar[length(gradAlphaSigmaScalar)-1];
	gradAlphaGamma[length(gradAlphaGamma)] <- 2*factor^3*coeff*
		(1-factor^2)^-1*(alpha[length(alpha)-1]^-1+coeff)^-2+
		gradAlphaPrev*gradAlphaGamma[length(gradAlphaGamma)-1]+
		2*factor*(alpha[length(alpha)-1]^-1+coeff)^-1;
	alpha[length(alpha)] <- beta^2+precs[i,1]+factor^2*
		(alpha[length(alpha)-1]^-1+coeff)^-1;	
	k[length(k)] <- alpha[length(alpha)]^-1*(data[i,1]*beta^2+factor^2*
		k[length(k)-1]*(alpha[length(alpha)-1]^-1+coeff)^-1);

	sigma <- t(X[i,])%*%Sigma%*%X[i,]+alpha[length(alpha)]^-1;
	#gradientSigma=0.5*(sigma^-1-(k[length(k)]-number)^2/sigma^2);
	#gradK=(k[length(k)]-number)/sigma;
	gradientSigma <- 0.5*(sigma^-1-(k[length(k)]-number)^2/sigma^2);
	gradK <- (k[length(k)]-number)/sigma;
	gradKPrev <- alpha[length(alpha)]^-1*factor^2*(alpha[length(alpha)-1]^-1+coeff)^-1;
	gradKAlpha <- -k[length(k)]/alpha[length(alpha)];
	gradKAlphaPrev <- factor^2*alpha[length(alpha)]^-1*k[length(k)-1]*
		alpha[length(alpha)-1]^-2*(alpha[length(alpha)-1]^-1+coeff)^-2;

	gradKBeta[length(gradKBeta)] <- 2*alpha[length(alpha)]^-1*
		data[i,1]*beta+gradKPrev*gradKBeta[length(gradKBeta)-1]+
		gradKAlpha*gradAlphaBeta[length(gradAlphaBeta)]+
		gradKAlphaPrev*gradAlphaBeta[length(gradAlphaBeta)-1];

	gradKGamma[length(gradKGamma)] <- 2*factor*alpha[length(alpha)]^-1*
		(alpha[length(alpha)-1]^-1+coeff)^-1*k[length(k)-1]*
		(1+factor^2*coeff*(alpha[length(alpha)-1]^-1+coeff)^-1*
		(1-factor^2)^-1)+gradKAlpha*gradAlphaGamma[length(gradAlphaGamma)]+
		gradKAlphaPrev*gradAlphaGamma[length(gradAlphaGamma)-1]+
		gradKPrev*gradKGamma[length(gradKGamma)-1];

	gradKSigmaScalar[length(gradKSigmaScalar)] <- -factor^2*
		alpha[length(alpha)]^-1*k[length(k)-1]*
		(alpha[length(alpha)-1]^-1+coeff)^-2*(1-factor^2)+
		gradKAlpha*gradAlphaSigmaScalar[length(gradAlphaSigmaScalar)]+
		gradKAlphaPrev*gradAlphaSigmaScalar[length(gradAlphaSigmaScalar)-1]+
		gradKPrev*gradKSigmaScalar[length(gradKSigmaScalar)-1];     

	gradKMuScalar[length(gradKMuScalar)] <- gradKPrev*
		gradKMuScalar[length(gradKMuScalar)-1];
	gradSigmaBeta <- 0;
	gradSigmaAlpha <- -alpha[length(alpha)]^-2;
	gradSigmaGamma <- 0;
	gradSigmaSigmaScalar <- 1;
	preGradBeta[i] <- preGradBeta[i]+gradientSigma*
		(gradSigmaAlpha*gradAlphaBeta[length(gradAlphaBeta)])+
		gradK*gradKBeta[length(gradKBeta)];
    
	preGradSigmaScalar[i] <- preGradSigmaScalar[i]+ gradientSigma*
		(gradSigmaSigmaScalar+gradSigmaAlpha*
		gradAlphaSigmaScalar[length(gradAlphaSigmaScalar)])+
		gradK*gradKSigmaScalar[length(gradKSigmaScalar)];

	preGradGamma <- preGradGamma+gradientSigma*(gradSigmaAlpha*
		gradAlphaGamma[length(gradAlphaGamma)])+gradK*gradKGamma[length(gradKGamma)];

	preGradMuScalar[i] <- preGradMuScalar[i]+
		gradK*gradKMuScalar[length(gradKMuScalar)]-(k[length(k)]-number)/sigma;
}


preGradSigma <- t(X)%*%(diag(as.vector(preGradSigmaScalar)))%*%(X);
postGradSigma <- 2*preSigma%*%preGradSigma;
gradSigma <- mat.or.vec(nEffectGenes,nTrans);

for (i in 1:length(R)) {
    gradSigma[R[i],C[i]] <- postGradSigma[R[i],C[i]];
}

gradSigma <- gradSigma[which(gradSigma!='0')]

gradDiagonal <- 2*diagonal*t(diag(preGradSigma))
gradMu <- (preGradMuScalar%*%X);
gradBeta <- rowSums(preGradBeta);
gradGamma <- as.vector((-1)*(preGradGamma*sin(2*gamma)));

f=cbind(gradBeta,gradGamma,gradMu, matrix(gradSigma,1,), gradDiagonal)

return(f)
}
