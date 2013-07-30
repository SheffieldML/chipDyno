#function f=chipStatMatrixInverterNoise(Sigma, gamma, beta,precs, x, npts);
#
#% CHIPSTATMATRIXINVERTERNOISE inverts block tridiagonal matrices for chipChip
#%
#%	Description:
#%	f=chipStatMatrixInverterNoise(Sigma, gamma, beta,precs, x, npts);
#%% 	chipStatMatrixInverterNoise.R version 0.1.0

chipStatMatrixInverterNoise = function(Sigma, gamma, beta, precs, x, npts){

####
#x=X[1,] # Temporary value for development/ transpose of colled value ??  
####


lambda = t(x) %*% Sigma %*% x
lambda = as.vector(lambda) ## need to cheek about its dimention
nTrans = length(x) ## This depends on dimention of 'x'
Y = Sigma %*% x
factor=cos(gamma)^2
UcoeffInvSigma=mat.or.vec(1,npts)
UcoeffXXT=mat.or.vec(1,npts)
LcoeffId=mat.or.vec(1,npts-1) 		# computes LU dec exploiting simple block
LcoeffXYT=mat.or.vec(1,npts-1)          # structure

#lambda=x'*Sigma*x;
#nTrans=size(x,1);
#%invSigma=pdinv(Sigma);
#Y=Sigma*x;
#factor=cos(gamma)^2;
#UcoeffInvSigma=zeros(1,npts);
#UcoeffXXT=zeros(1,npts);
#LcoeffId=zeros(1,npts-1); %computes LU dec exploiting simple block
#LcoeffXYT=zeros(1,npts-1);                        %structure

UcoeffXXT[1]=(beta^-2+precs[1]^-1)^-1;
UcoeffInvSigma[1]=(1-factor^2)^-1;
LcoeffXYT[1]=factor*(1-factor^2)^-1*(beta^-2+precs[1]^-1)^-1/(UcoeffInvSigma[1]* (UcoeffInvSigma[1]+(beta^-2+precs[1]^-1)^-1*lambda));
LcoeffId[1]=-factor;

#UcoeffXXT(1)=(beta^-2+precs(1)^-1)^-1;
#UcoeffInvSigma(1)=(1-factor^2)^-1;
#LcoeffXYT(1)=factor*(1-factor^2)^-1*(beta^-2+precs(1)^-1)^-1/(UcoeffInvSigma(1)* ...
#                                             (UcoeffInvSigma(1)+(beta^-2+precs(1)^-1)^-1*lambda));
#LcoeffId(1)=-factor;




for (i in 2:(npts-1)){
	UcoeffXXT[i]=(beta^-2+precs[i]^-1)^-1+factor*(1-factor^2)^-1*LcoeffXYT[(i-1)];
	UcoeffInvSigma[i]=(1+factor^2)*(1-factor^2)^-1+factor*(1-factor^2)^-1*LcoeffId[(i-1)];
	LcoeffXYT[i]=factor*(1-factor^2)^-1*UcoeffXXT[i]/(UcoeffInvSigma[i]*(UcoeffInvSigma[i]+UcoeffXXT[i]*lambda));
	LcoeffId[i]=-factor*(1-factor^2)^-1*UcoeffInvSigma[i]^-1;
}


#for i=2:npts-1
#  UcoeffXXT(i)=(beta^-2+precs(i)^-1)^-1+factor*(1-factor^2)^-1*LcoeffXYT(i-1);
#  UcoeffInvSigma(i)=(1+factor^2)*(1-factor^2)^-1+factor*(1- ...
#                                                    factor^2)^-1*LcoeffId(i-1);
#  LcoeffXYT(i)=factor*(1-factor^2)^-1*UcoeffXXT(i)/ ...
#      (UcoeffInvSigma(i)*(UcoeffInvSigma(i)+UcoeffXXT(i)*lambda));
#  
#  
#  
#  LcoeffId(i)=-factor*(1-factor^2)^-1*UcoeffInvSigma(i)^-1;
#end

UcoeffXXT[ncol(UcoeffXXT)]=(beta^-2+precs[length(precs)]^-1)^-1+factor*(1-factor^2)^-1*LcoeffXYT[ncol(LcoeffXYT)];
UcoeffInvSigma[ncol(UcoeffInvSigma)]=(1-factor^2)^-1+factor*(1-factor^2)^-1*LcoeffId[ncol(LcoeffId)];
#%lambda=Y'*x;
invL.XYT=mat.or.vec(npts,npts); #%computes the inverse of the L bit
invL.Id=mat.or.vec(npts,npts);


#UcoeffXXT(end)=(beta^-2+precs(end)^-1)^-1+factor*(1-factor^2)^-1*LcoeffXYT(end);
#UcoeffInvSigma(end)=(1-factor^2)^-1+factor*(1-factor^2)^-1*LcoeffId(end);
#%lambda=Y'*x;
#invL.XYT=zeros(npts,npts);%computes the inverse of the L bit
#invL.Id=zeros(npts,npts);


for (i in 1:npts){
	invL.Id[i,i]=1;
}

#for i=1:npts
#    invL.Id(i,i)=1;
#end

for (i in 2:npts) {
	invL.Id[i,(i-1)]=(-1)^(2*i-1)*LcoeffId[i-1];
	invL.XYT[i,(i-1)]=(-1)^(2*i-1)*LcoeffXYT[i-1];
}

#for i=2:npts
#    invL.Id(i,i-1)=(-1)^(2*i-1)*LcoeffId(i-1);
#    invL.XYT(i,i-1)=(-1)^(2*i-1)*LcoeffXYT(i-1);
#end

for (i in 3 : npts) {
	for (j in 1:(i-2)){
		invL.Id[i,j]=invL.Id[(i-1),j]*invL.Id[i,(i-1)];
		invL.XYT[i,j]=(invL.Id[(i-1),j]*invL.XYT[i,(i-1)]+ invL.Id[i,(i-1)]*invL.XYT[(i-1),j]+ invL.XYT[i,(i-1)]*invL.XYT[(i-1),j]*lambda);

	}
}


#for i=3:npts
#    for j=1:i-2
#    invL.Id(i,j)=invL.Id(i-1,j)*invL.Id(i,i-1);
#    invL.XYT(i,j)=(invL.Id(i-1,j)*invL.XYT(i,i-1)+...
#        invL.Id(i,i-1)*invL.XYT(i-1,j)+...
#        invL.XYT(i,i-1)*invL.XYT(i-1,j)*lambda);

#    end
#end


invU.Sigma=mat.or.vec(npts,npts);
invU.YYT=mat.or.vec(npts,npts);
for (i in 1:(npts-1)){
	invU.Sigma[i,i]=-LcoeffId[i]*(1-factor^2)/factor;
	invU.YYT[i,i]=-LcoeffXYT[i]*(1-factor^2)/factor;
}


#invU.Sigma=zeros(npts,npts);
#invU.YYT=zeros(npts,npts);
#for i=1:npts-1
#    invU.Sigma(i,i)=-LcoeffId(i)*(1-factor^2)/factor;
#    invU.YYT(i,i)=-LcoeffXYT(i)*(1-factor^2)/factor;
#end


invU.Sigma[nrow(invU.Sigma),ncol(invU.Sigma)]=UcoeffInvSigma[length(UcoeffInvSigma)]^-1;
invU.YYT[nrow(invU.YYT),ncol(invU.YYT)]=-UcoeffXXT[length(UcoeffXXT)]/(UcoeffInvSigma[length(UcoeffInvSigma)]*(UcoeffInvSigma[length(UcoeffInvSigma)]+UcoeffXXT[length(UcoeffXXT)]*lambda));


#invU.Sigma(end,end)=UcoeffInvSigma(end)^-1;
#invU.YYT(end,end)=-UcoeffXXT(end)/(UcoeffInvSigma(end)* ...
#                                              (UcoeffInvSigma(end)+ ...
#                                               UcoeffXXT(end)*lambda));


for (i in 1:(npts-1)){
	for (j in 1:i){
		invU.Sigma[(npts-i),(npts-j+1)]=factor*(1-factor^2)^-1* invU.Sigma[(npts-i+1),(npts-j+1)]*invU.Sigma[(npts-i),(npts-i)];
        	invU.YYT[(npts-i),(npts-j+1)]=factor*(1-factor^2)^-1*(invU.Sigma[(npts-i+1),(npts-j+1)]*invU.YYT[(npts-i),(npts-i)]+invU.Sigma[(npts-i),(npts-i)]*invU.YYT[(npts-i+1),(npts-j+1)]+ invU.YYT[(npts-i+1),(npts-j+1)]*invU.YYT[(npts-i),(npts-i)]*lambda);
	}
}


#for i=1:npts-1
#    for j=1:i
#        invU.Sigma(npts-i,npts-j+1)=factor*(1-factor^2)^-1*...
#            invU.Sigma(npts-i+1,npts-j+1)*invU.Sigma(npts-i,npts-i);
#        invU.YYT(npts-i,npts-j+1)=factor*(1-factor^2)^-1*...
#            (invU.Sigma(npts-i+1,npts-j+1)*invU.YYT(npts-i,npts-i)+...
#        invU.Sigma(npts-i,npts-i)*invU.YYT(npts-i+1,npts-j+1)+...
#        invU.YYT(npts-i+1,npts-j+1)*invU.YYT(npts-i,npts-i)*lambda);
#    end
#end


invC.Sigma=mat.or.vec(npts,npts); #%computes the inverses of the matrix;
invC.YYT=mat.or.vec(npts,npts);

#invC.Sigma=zeros(npts,npts);%computes the inverses of the matrix;
#invC.YYT=zeros(npts,npts);


for (i in 1:npts) {
	for (j in 1:npts) {
        	invC.Sigma[i,j]=invU.Sigma[i,]%*%invL.Id[,j];
		invC.YYT[i,j]=invU.Sigma[i,]%*%invL.XYT[,j]+invU.YYT[i,]%*% invL.Id[,j]+lambda*invU.YYT[i,]%*%invL.XYT[,j];
	}
}

#for i=1:npts
#    for j=1:npts
#        invC.Sigma(i,j)=invU.Sigma(i,:)*invL.Id(:,j);
#        invC.YYT(i,j)=invU.Sigma(i,:)*invL.XYT(:,j)+invU.YYT(i,:)* ...
#            invL.Id(:,j)+lambda*invU.YYT(i,:)*invL.XYT(:,j);
#    end
#end



invC <- list(invC.Sigma, invC.YYT)
#f=invC;

## 'invC.Sigma' will be assigned by 'invC[[1]]'; invC.Sigma=invC[[1]]
## 'invC.YYT' will be assigned by 'invC[[2]]';invC.YYT=invC[[2]]

return(invC)
}

