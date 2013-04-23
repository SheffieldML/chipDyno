function expectations=chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);

% CHIPDYNOSTATPOSTEST computes posterior expectations
% CHIPDYNO toolbox
% chipDynoStatPostEst.m version 1.4
% FORMAT expectations=chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);
% DESC computes posterior expectations
% ARG data : point estimate of the expression level
% ARG x : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix of TFA
% ARG beta :
% ARG gamma : degree of temporal continuity
% ARG mu : mean value of the transcription factor activity
% RETURN expectations: computes posterior expectations
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoStatPostEstNoise

npts=size(data,2);
nTrans=size(x,1);
invC=chipStatMatrixInverter(Sigma, gamma, beta,x,npts);
Y=Sigma*x;
lambda=Y'*x;
factor=cos(gamma)^2;
coeff=x'*mu';
YYT=Y*Y';
invSigma=pdinv(Sigma);
%Z=invSigma*mu';
Mean=[beta^2*sum(data.*(invC.Sigma(1,:)+lambda*invC.YYT(1,:)))*Y'+ ...
      (1+factor)^-1*(invC.Sigma(1,1)*mu+coeff*invC.YYT(1,1)*Y')+ ...
     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(1,2:end-1))*mu+coeff* ...
        sum(invC.YYT(1,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(1, ...
                                                  end)*mu+coeff* ...
     invC.YYT(1,end)*Y')];

for i=2:npts-1

  Mean=[Mean;beta^2*sum(data.*(invC.Sigma(i,:)+lambda*invC.YYT(i,:)))*Y'+ ...
        (1+factor)^-1*(invC.Sigma(i,1)*mu+coeff*invC.YYT(i,1)*Y')+ ...
     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(i,2:end-1))*mu+coeff* ...
        sum(invC.YYT(i,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(i, ...
                                                  end)*mu+coeff* ...
     invC.YYT(i,end)*Y')];
end
Mean=[Mean;beta^2*sum(data.*(invC.Sigma(end,:)+lambda*invC.YYT(end,:)))*Y'+ ...
      (1+factor)^-1*(invC.Sigma(end,1)*mu+coeff*invC.YYT(end,1)*Y')+ ...
     (1-factor)*(1+factor)^-1*(sum(invC.Sigma(end,2:end-1))*mu+coeff* ...
        sum(invC.YYT(end,2:end-1))*Y')+(1+factor)^-1*(invC.Sigma(end, ...
                                                  end)*mu+coeff* ...
     invC.YYT(end,end)*Y')];
        
       
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
