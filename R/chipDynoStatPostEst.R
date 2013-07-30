#function expectations=chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);

# CHIPDYNOSTATPOSTEST computes posterior expectations
#
#	Description:
#	expectations=chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);
# 	chipDynoStatPostEst.R version 0.1.0

chipDynoStatPostEst = function(data,x,Sigma,beta,gamma,mu) {

npts=length(data);
nTrans=length(x);
#npts=ncol(data);
#nTrans=nrow(x);

source("chipStatMatrixInverter.R")
invC=chipStatMatrixInverter(Sigma, gamma, beta, x, npts); 
## 'invC.Sigma' will be assigned by 'invC[[1]]'; invC.Sigma=invC[[1]]
## 'invC.YYT' will be assigned by 'invC[[2]]';invC.YYT=invC[[2]]

invC.Sigma=invC[[1]]
invC.YYT=invC[[2]]

## need to check the dimention of 'mu'
#mu=t(mu)
###

Y=Sigma%*%x;
lambda=as.vector(t(Y)%*%x);
factor=cos(gamma)^2;
coeff=(t(x)) %*% (mu); #?????????????? coeff=(t(x)) %*% t(mu)
YYT=Y%*%t(Y);

#npts=size(data,2);
#nTrans=size(x,1);
#invC=chipStatMatrixInverter(Sigma, gamma, beta,x,npts);
#Y=Sigma*x;
#lambda=Y'*x;
#factor=cos(gamma)^2;
#coeff=x'*mu';
#YYT=Y*Y';

#require("mnormt") #mnormt package is required to inverse the positive definite matrix
#Sigma = data.matrix(Sigma) # 'Sigma' have to be a data frame
#invSigma = pd.solve(Sigma)

#invSigma=pdinv(Sigma);
# %Z=invSigma*mu';

#sum(invC.Sigma[1,2:(ncol(invC.Sigma)-1)])
#sum(invC.YYT[1,2:(ncol(invC.YYT)-1)])

Mean=beta^2*sum(data*(invC.Sigma[1,]+lambda*invC.YYT[1,]))%*%t(Y)+
      (1+factor)^-1*(invC.Sigma[1,1]*mu+coeff%*%invC.YYT[1,1]%*%t(Y))+
      (1-factor)*(1+factor)^-1*((sum(invC.Sigma[1,2:(ncol(invC.Sigma)-1)]))*mu+ 
        coeff%*%sum(invC.YYT[1,2:(ncol(invC.YYT)-1)])%*%t(Y))+
	(1+factor)^-1*(invC.Sigma[1,ncol(invC.Sigma)]*mu+coeff%*%invC.YYT[1,ncol(invC.YYT)]%*%t(Y));


#Mean=[beta^2*sum(data.*(invC.Sigma(1,:)+lambda*invC.YYT(1,:)))*Y'+ ...
#      (1+factor)^-1*(invC.Sigma(1,1)*mu+coeff*invC.YYT(1,1)*Y')+ ...
#     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(1,2:end-1))*mu+coeff* ...
#        sum(invC.YYT(1,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(1, ...
#                                                  end)*mu+coeff* ...
#     invC.YYT(1,end)*Y')];


for (i in 2:(npts-1)){
	Mean_temp=beta^2*sum(data*(invC.Sigma[i,]+lambda*invC.YYT[i,]))%*%t(Y)+
	(1+factor)^-1*(invC.Sigma[i,1]*mu+coeff%*%invC.YYT[i,1]%*%t(Y))+
	(1-factor)*(1+factor)^-1*((sum(invC.Sigma[i,2:(ncol(invC.Sigma)-1)]))*mu+
	  coeff%*%sum(invC.YYT[i,2:(ncol(invC.Sigma)-1)])%*%t(Y))+
	  (1+factor)^-1*(invC.Sigma[i,ncol(invC.Sigma)]*mu+coeff%*%invC.YYT[i,ncol(invC.YYT)]%*%t(Y));
	
	Mean=rbind(Mean[1:nrow(Mean),],Mean_temp[1,])
#	Mean=rbind(Mean[,],Mean_temp[1,])
}

#for i=2:npts-1

#  Mean=[Mean;beta^2*sum(data.*(invC.Sigma(i,:)+lambda*invC.YYT(i,:)))*Y'+ ...
#        (1+factor)^-1*(invC.Sigma(i,1)*mu+coeff*invC.YYT(i,1)*Y')+ ...
#     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(i,2:end-1))*mu+coeff* ...
#        sum(invC.YYT(i,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(i, ...
#                                                  end)*mu+coeff* ...
#     invC.YYT(i,end)*Y')];
#end

Mean_end=beta^2*sum(data*(invC.Sigma[ncol(invC.Sigma),]+lambda*invC.YYT[ncol(invC.YYT),]))%*%t(Y)+
      (1+factor)^-1*(invC.Sigma[ncol(invC.Sigma),1]*mu+coeff%*%invC.YYT[ncol(invC.YYT),1]%*%t(Y))+
      (1-factor)*(1+factor)^-1*((sum(invC.Sigma[ncol(invC.Sigma),(2:(ncol(invC.Sigma)-1))]))*mu+
      coeff%*% sum(invC.YYT[ncol(invC.YYT),2:(ncol(invC.YYT)-1)])%*%t(Y))+
	(1+factor)^-1*(invC.Sigma[ncol(invC.Sigma),ncol(invC.Sigma)]*mu+coeff%*%invC.YYT[ncol(invC.YYT),ncol(invC.YYT)]%*%t(Y));

Mean=rbind(Mean[,],Mean_end[1,])


#Mean=[Mean;beta^2*sum(data.*(invC.Sigma(end,:)+lambda*invC.YYT(end,:)))*Y'+ ...
#      (1+factor)^-1*(invC.Sigma(end,1)*mu+coeff*invC.YYT(end,1)*Y')+ ...
#     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(end,2:end-1))*mu+coeff* ...
#        sum(invC.YYT(end,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(end, ...
#                                                  end)*mu+coeff* ...
#     invC.YYT(end,end)*Y')];
#


expectations.b=Mean;
gigio=which(x!=0);
expectations.tfError=array(0,dim=c(npts,sum(x)))
expectations.tfErrorDiffs=array(0,dim=c(npts,npts,sum(x)))
preDiffs=array(0,dim=c(npts,npts));

#expectations.b=Mean;
#gigio=find(x);
#expectations.tfError=zeros(npts,sum(x));
#expectations.tfErrorDiffs=zeros(npts,npts, sum(x));
#preDiffs=zeros(npts,npts);


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


#for i=1:sum(x)
#    postCov=invC.Sigma*Sigma(gigio(i),gigio(i))+invC.YYT* ...
#            Y(gigio(i))^2;
#    %[var,u,lambda]=ppca(postCov,1);
#    auxMat=postCov-(ones(1,npts)*postCov*ones(npts,1))* ...
#                                        ones(npts,npts)/npts^2;
#    expectations.tfError(:,i)=sqrt(diag(postCov));
#    for j=1:npts-1
#        for l=j+1:npts
#           preDiffs(j,l)=sqrt((auxMat(j,j)+auxMat(l, ...
#                                                          l)-2* ...
#                                           auxMat(j,l))/2);
#        end
#    end
#    expectations.tfErrorDiffs(:,:,i)=preDiffs+preDiffs'+eye(npts);
#end

expectations = list(expectations.b, expectations.tfError, expectations.tfErrorDiffs)

return(expectations)
}
