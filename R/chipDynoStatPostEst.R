# CHIPDYNOSTATPOSTEST computes posterior expectations
# CHIPDYNO toolbox
# chipDynoStatPostEst.R version 1.0.1
# FORMAT chipDynoStatPostEst <- function(data,x,Sigma,beta,gamma,mu)
# DESC computes posterior expectations
# ARG data : point estimate of the expression level
# ARG x : connectivity measurement between genes and transcription factors
# ARG Sigma : prior covariance matrix of TFA
# ARG beta :
# ARG gamma : degree of temporal continuity
# ARG mu : mean value of the transcription factor activity
# RETURN expectations: computes posterior expectations
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipDynoStatPostEstNoise

chipDynoStatPostEst <- function(data,x,Sigma,beta,gamma,mu) {

npts=length(data);
nTrans=length(x);
source("chipStatMatrixInverter.R")
invC=chipStatMatrixInverter(Sigma, gamma, beta, x, npts); 
invC.Sigma=invC[[1]]
invC.YYT=invC[[2]]

Y=Sigma%*%x;
lambda=as.vector(t(Y)%*%x);
factor=cos(gamma)^2;
coeff=(t(x)) %*% (mu)
YYT=Y%*%t(Y);

Mean=beta^2*sum(data*(invC.Sigma[1,]+lambda*invC.YYT[1,]))%*%t(Y)+
      (1+factor)^-1*(invC.Sigma[1,1]*mu+coeff%*%invC.YYT[1,1]%*%t(Y))+
      (1-factor)*(1+factor)^-1*((sum(invC.Sigma[1,2:(ncol(invC.Sigma)-1)]))*mu+ 
        coeff%*%sum(invC.YYT[1,2:(ncol(invC.YYT)-1)])%*%t(Y))+
	(1+factor)^-1*(invC.Sigma[1,ncol(invC.Sigma)]*mu+coeff%*%invC.YYT[1,ncol(invC.YYT)]%*%t(Y));

for (i in 2:(npts-1)){
	Mean_temp=beta^2*sum(data*(invC.Sigma[i,]+lambda*invC.YYT[i,]))%*%t(Y)+
	(1+factor)^-1*(invC.Sigma[i,1]*mu+coeff%*%invC.YYT[i,1]%*%t(Y))+
	(1-factor)*(1+factor)^-1*((sum(invC.Sigma[i,2:(ncol(invC.Sigma)-1)]))*mu+
	  coeff%*%sum(invC.YYT[i,2:(ncol(invC.Sigma)-1)])%*%t(Y))+
	  (1+factor)^-1*(invC.Sigma[i,ncol(invC.Sigma)]*mu+coeff%*%invC.YYT[i,ncol(invC.YYT)]%*%t(Y));
	
	Mean=rbind(Mean[1:nrow(Mean),],Mean_temp[1,])

Mean_end=beta^2*sum(data*(invC.Sigma[ncol(invC.Sigma),]+lambda*invC.YYT[ncol(invC.YYT),]))%*%t(Y)+
      (1+factor)^-1*(invC.Sigma[ncol(invC.Sigma),1]*mu+coeff%*%invC.YYT[ncol(invC.YYT),1]%*%t(Y))+
      (1-factor)*(1+factor)^-1*((sum(invC.Sigma[ncol(invC.Sigma),(2:(ncol(invC.Sigma)-1))]))*mu+
      coeff%*% sum(invC.YYT[ncol(invC.YYT),2:(ncol(invC.YYT)-1)])%*%t(Y))+
	(1+factor)^-1*(invC.Sigma[ncol(invC.Sigma),ncol(invC.Sigma)]*mu+coeff%*%invC.YYT[ncol(invC.YYT),ncol(invC.YYT)]%*%t(Y));

Mean=rbind(Mean[,],Mean_end[1,])

expectations.b=Mean;
gigio=which(x!=0);
expectations.tfError=array(0,dim=c(npts,sum(x)))
expectations.tfErrorDiffs=array(0,dim=c(npts,npts,sum(x)))
preDiffs=array(0,dim=c(npts,npts));

for (i in 1: sum(x)){
	postCov=invC.Sigma*Sigma[gigio[i],gigio[i]]+invC.YYT*Y[gigio[i]]^2;
    	#%[var,u,lambda]=ppca(postCov,1);
	auxMat=postCov-as.vector(matrix(1,1,npts)%*%postCov%*%matrix(1,npts,1))*matrix(1,npts,npts)/npts^2;
	expectations.tfError[,i]=sqrt(diag(postCov));
	for (j in 1: (npts-1)){
		for ( l in (j+1) : npts){
			preDiffs[j,l]=sqrt((auxMat[j,j]+auxMat[l,l]-2*auxMat[j,l])/2);
		}
	}
	expectations.tfErrorDiffs[ , ,i]=preDiffs+t(preDiffs)+diag(npts);
}

expectations = list(expectations.b, expectations.tfError, expectations.tfErrorDiffs)

return(expectations)
}