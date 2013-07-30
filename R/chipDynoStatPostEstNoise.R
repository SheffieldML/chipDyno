#function expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);
#
#% CHIPDYNOSTATPOSTESTNOISE computes posterior expectations.
#%
#%	Description:
#%	expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);
#%% 	chipDynoStatPostEstNoise.m version 1.4

chipDynoStatPostEstNoise = function(data,x,Sigma,beta,precs,gamma,mu) {

npts=length(data);
nTrans=nrow(x);
source("chipStatMatrixInverterNoise.R")
invC=chipStatMatrixInverterNoise(Sigma, gamma, beta, precs, x,npts);

#npts=size(data,2);
#nTrans=size(x,1);
#invC=chipStatMatrixInverterNoise(Sigma, gamma, beta,precs,x,npts);

## 'invC.Sigma' will be assigned by 'invC[[1]]'; invC.Sigma=invC[[1]]
## 'invC.YYT' will be assigned by 'invC[[2]]';invC.YYT=invC[[2]]

invC.Sigma=invC[[1]]
invC.YYT=invC[[2]]

Y=Sigma%*%x;
lambda=as.vector(t(Y)%*%x);
factor=cos(gamma)^2;
coeff=(t(x)) %*% (mu); #?????????????? coeff=(t(x)) %*% t(mu)
coeff=coeff[1,1]

#YYT=Y%*%t(Y);
#require("mnormt") #mnormt package is required to inverse the positive definite matrix
#Sigma = data.matrix(Sigma) # 'Sigma' have to be a data frame
#invSigma = pd.solve(Sigma)
#Z=invSigma%*%t(mu);


#Y=Sigma*x;
#lambda=Y'*x;
#factor=cos(gamma)^2;
#coeff=x'*mu';
#%YYT=Y*Y';
#%invSigma=pdinv(Sigma);
#%Z=invSigma*mu';

Mean = list()

for (i in 1:npts){
	tempMean=sum((beta^-2*array(1,dim<-c(1,npts))+precs^-1)^-1*data*(invC.Sigma[i,]+lambda*invC.YYT[i,]))*t(Y)+
	(1+factor)^-1*(invC.Sigma[i,1]*mu+coeff*invC.YYT[i,1]%*%t(Y))+
	(1-factor)*(1+factor)^-1*(sum(invC.Sigma[i,2:(ncol(invC.Sigma)-1)])*mu+coeff* sum(invC.YYT[i,2:(ncol(invC.YYT)-1)])%*%t(Y))+
	(1+factor)^-1*(invC.Sigma[i,ncol(invC.Sigma)]*mu+coeff*invC.YYT[i,ncol(invC.YYT)]*t(Y));
	if (i==1){
		Mean=tempMean;
	} else
		Mean=rbind(Mean[1:nrow(Mean),],tempMean[1,])
}


#Mean=[];
#for i=1:npts

#  Mean=[Mean;sum((beta^-2*ones(1,npts)+precs.^-1).^-1.*data.*(invC.Sigma(i,:)+lambda*invC.YYT(i,:)))*Y'+ ...
#        (1+factor)^-1*(invC.Sigma(i,1)*mu+coeff*invC.YYT(i,1)*Y')+ ...
#     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(i,2:end-1))*mu+coeff* ...
#        sum(invC.YYT(i,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(i, ...
#                                                  end)*mu+coeff* ...
#     invC.YYT(i,end)*Y')];
#end
       
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
