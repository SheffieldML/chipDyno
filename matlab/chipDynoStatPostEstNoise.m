function expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);

% CHIPDYNOSTATPOSTESTNOISE computes posterior expectations.
%
%	Description:
%	expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);
%% 	chipDynoStatPostEstNoise.m version 1.4



npts=size(data,2);
nTrans=size(x,1);
invC=chipStatMatrixInverterNoise(Sigma, gamma, beta,precs,x,npts);
Y=Sigma*x;
lambda=Y'*x;
factor=cos(gamma)^2;
coeff=x'*mu';
%YYT=Y*Y';
%invSigma=pdinv(Sigma);
%Z=invSigma*mu';
Mean=[];

for i=1:npts

  Mean=[Mean;sum((beta^-2*ones(1,npts)+...
      precs.^-1).^-1.*data.*(invC.Sigma(i,:)+lambda*invC.YYT(i,:)))*Y'+ ...
        (1+factor)^-1*(invC.Sigma(i,1)*mu+coeff*invC.YYT(i,1)*Y')+ ...
     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(i,2:end-1))*mu+coeff* ...
        sum(invC.YYT(i,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(i, ...
                                                  end)*mu+coeff* ...
     invC.YYT(i,end)*Y')];
end

       
expectations.b=Mean;
gigio=find(x);
expectations.tfError=zeros(npts,sum(x));
expectations.tfErrorDiffs=zeros(npts,npts, sum(x));
preDiffs=zeros(npts,npts);
for i=1:sum(x)
    postCov=invC.Sigma*Sigma(gigio(i),gigio(i))+invC.YYT* ...
            Y(gigio(i))^2;
    %[var,u,lambda]=ppca(postCov,1);
    auxMat=postCov-(ones(1,npts)*postCov*ones(npts,1))* ...
                                        ones(npts,npts)/npts^2;
    expectations.tfError(:,i)=sqrt(diag(postCov));
    for j=1:npts-1
        for l=j+1:npts
           preDiffs(j,l)=sqrt((auxMat(j,j)+auxMat(l, ...
                                                          l)-2* ...
                                           auxMat(j,l))/2);
        end
    end
    expectations.tfErrorDiffs(:,:,i)=preDiffs+preDiffs'+eye(npts);
end
     
     