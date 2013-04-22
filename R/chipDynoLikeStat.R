#function f=chipDynoLikeStat(params,data,X,nEffectGenes,R,C);

#% CHIPDYNOLIKESTAT marginal likelihood for chipChip dynamical model
#%
#%	Description:
#%	f=chipDynoLikeStat(params,data,X,nEffectGenes,R,C);
#%% 	chipDynoLikeStat.m version 1.4
chipDynoLikeStat = function(params,data, X,nEffectGenes,R,C){
##

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


#nGenes=size(data,1);
#nTrans=size(X,2);
#npts=size(data,2);
#beta=params(1);
#gamma=params(2);
#mu=params(3:2+nTrans)';
#V=params(3+nTrans:end-nTrans)';
#preSigma=sparse(R,C,V,nEffectGenes,nTrans);
#diagonal=params(end-nTrans+1:end);
#Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
#factor=cos(gamma)^2;
#preLike=zeros(1, nGenes);
#k=zeros(1,npts);
#alpha=zeros(1,npts);

for (i in 1: nGenes){
	coeff=(1-factor^2)*X[i,] %*% Sigma %*% X[i,];
	number=X[i,]%*%mu;
	alpha[1]=beta^2;
	k[1]=(data[i,npts]-(1-factor)*number)/factor;
	sigma=beta^-2+(alpha[1]^-1+coeff)/factor^2;
	#update=0.5*log(sigma)+ 0.5*((k[1]-data[i,npts-1])^2)/(sigma); 
	preLike[i]=preLike[i]+ (0.5*log(sigma)+ 0.5*((k[1]-data[i,npts-1])^2)/(sigma))	
	# Have to check the preLike...
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




#for i=1:nGenes
#  coeff=(1-factor^2)*X(i,:)*Sigma*X(i,:)';
#  number=X(i,:)*mu;
#  alpha(1)=beta^2;
#  k(1)=(data(i,npts)-(1-factor)*number)/factor;
#  sigma=beta^-2+(alpha(1)^-1+coeff)/factor^2;
#  preLike(i)=preLike(i)+ 0.5*log(sigma)+...
#       0.5*((k(1)-data(i,npts-1))^2)/(sigma);
#  for j=2:npts-1
#    alpha(j)=(beta^2+factor^2*(alpha(j-1)^-1+coeff)^-1);
#    k(j)=alpha(j)^-1*(data(i,npts-j+1)*beta^2+factor^2*k(j-1)*(alpha(j-1)^-1+coeff)^-1);
#    sigma=beta^-2+(alpha(j)^-1+coeff)/factor^2;
#    preLike(i)=preLike(i)+ 0.5*log(sigma)+...
#      0.5*((k(j)-data(i,npts-j))^2)/(sigma);
#  end
#  alpha(end)=(beta^2+factor^2*(alpha(end-1)^-1+coeff)^-1);
#  k(end)=alpha(end)^-1*(data(i,1)*beta^2+factor^2*k(end-1)*(alpha(end-1)^-1+coeff)^-1);
#  preLike(i)=preLike(i)+0.5*log(alpha(end)^-1+X(i,:)*Sigma*X(i,:)')+...
#      0.5*((k(end)-number)^2)/(alpha(end)^-1+X(i,:)*Sigma*X(i,:)');
#end
#f=sum(preLike);

f=rowSums(preLike);

##test
return(f)
}


