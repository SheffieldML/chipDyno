function f=chipDynoLikeStatGrad(params,data,X,nEffectGenes,R,C);

% CHIPDYNOLIKESTATGRAD gradient of CHIPDYNOLIKESTAT.
% CHIPDYNO toolbox
% chipDynoLikeStatGrad.m version 1.4
% FORMAT f=chipDynoLikeStatGrad(params,data,X,nEffectGenes,R,C);
% DESC compute the gradient of chipDynoLikeStat for chipChip dynamical model
% ARG params: concatenated vector of multiple parameters(beta, gamma, 
% initial mean of the transcription factors, and 
% a vector to create diagonal matrix used to reduce the sparsity of covariance)
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG nEffectGenes : effectice gene name
% ARG R, C : same length integer vectors specifying the row and column 
% indices of the non-zero entries of the sparce matrix
% RETURN f : gradient of concatenated parameters (beta,gamma, mu, Sigma, diagonal)
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoLikeStatGradNoise

nGenes=size(data,1);
nTrans=size(X,2);
npts=size(data,2);
beta=params(1);
gamma=params(2);
factor=cos(gamma)^2;
mu=params(3:2+nTrans)';
V=params(3+nTrans:end-nTrans)';
preSigma=sparse(R,C,V,nEffectGenes,nTrans);
diagonal=params(end-nTrans+1:end);
Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
preGradBeta=zeros(1, nGenes);
preGradSigmaScalar=zeros(1, nGenes);
preGradMuScalar=zeros(1,nGenes);
preGradGamma=0;
k=zeros(1,npts);
alpha=zeros(1,npts);
gradAlphaBeta=zeros(1,npts);
gradAlphaGamma=zeros(1,npts);
gradAlphaSigmaScalar=zeros(1,npts);
gradKBeta=zeros(1,npts);
gradKGamma=zeros(1,npts);
gradKSigmaScalar=zeros(1,npts);
gradKMuScalar=zeros(1,npts);
for i=1:nGenes
  number=X(i,:)*mu;
  coeff=(1-factor^2)*X(i,:)*Sigma*X(i,:)';
 
  alpha(1)=beta^2;
  k(1)=(data(i,npts)-(1-factor)*number)/factor;
  sigma=beta^-2+(alpha(1)^-1+coeff)/factor^2;
  gradientSigma=0.5*(sigma^-1-(k(1)-data(i,npts-1))^2/sigma^2);
  gradK=(k(1)-data(i,npts-1))/sigma;
  gradAlphaBeta(1)=2*beta;
  gradAlphaSigmaScalar(1)=0;
  gradAlphaGamma(1)=0;
  gradKBeta(1)=0;
  gradKSigmaScalar(1)=0;
  gradKGamma(1)=-k(1)/factor+number/factor;
  gradKMuScalar(1)=-(1-factor)/factor;
  preGradBeta(i)=gradientSigma*(-2*beta^-3-alpha(1)^-2*gradAlphaBeta(1)/factor^2);
  preGradGamma=preGradGamma+gradK*gradKGamma(1)+gradientSigma*(-2*(alpha(1)^-1+coeff)/factor^3-2*coeff/(factor*(1-factor^2)));
  preGradSigmaScalar(i)=gradientSigma*(1-factor^2)/factor^2;
  preGradMuScalar(i)=gradK*gradKMuScalar(1);
  for j=2:npts-1
    gradAlphaPrev=factor^2*(alpha(j-1)^-1+coeff)^-2*alpha(j-1)^-2;  
    gradAlphaBeta(j)=2*beta+gradAlphaPrev*gradAlphaBeta(j-1);
    gradAlphaSigmaScalar(j)=-factor^2*(1-factor^2)*(alpha(j-1)^-1+coeff)^-2+gradAlphaPrev*gradAlphaSigmaScalar(j-1);
    gradAlphaGamma(j)=2*factor^3*coeff*(1-factor^2)^-1*(alpha(j-1)^-1+coeff)^-2+gradAlphaPrev*gradAlphaGamma(j-1)+2*factor*(alpha(j-1)^-1+coeff)^-1;
    alpha(j)=(beta^2+factor^2*(alpha(j-1)^-1+coeff)^-1);
    k(j)=alpha(j)^-1*(data(i,npts-j+1)*beta^2+factor^2*k(j-1)*(alpha(j-1)^-1+coeff)^-1);
    sigma=beta^-2+(alpha(j)^-1+coeff)/factor^2; 
    gradientSigma=0.5*(sigma^-1-(k(j)-data(i,npts-j))^2/sigma^2);
    gradK=(k(j)-data(i,npts-j))/sigma;
    gradKPrev=alpha(j)^-1*factor^2*(alpha(j-1)^-1+coeff)^-1;
    gradKAlpha=-k(j)/alpha(j);
    gradKAlphaPrev=factor^2*alpha(j)^-1*k(j-1)*alpha(j-1)^-2*(alpha(j-1)^-1+coeff)^-2;
    gradKBeta(j)=2*alpha(j)^-1*data(i,npts-j+1)*beta+gradKPrev*gradKBeta(j-1)+gradKAlpha*gradAlphaBeta(j)+gradKAlphaPrev*gradAlphaBeta(j-1);
    gradKGamma(j)=2*factor*alpha(j)^-1*(alpha(j-1)^-1+coeff)^-1*k(j-1)*(1+factor^2*coeff*(alpha(j-1)^-1+coeff)^-1*(1-factor^2)^-1)+...
        gradKAlpha*gradAlphaGamma(j)+gradKAlphaPrev*gradAlphaGamma(j-1)+gradKPrev*gradKGamma(j-1);

    gradKSigmaScalar(j)=-factor^2*alpha(j)^-1*k(j-1)*(alpha(j-1)^-1+coeff)^-2*(1-factor^2)+...
        gradKAlpha*gradAlphaSigmaScalar(j)+gradKAlphaPrev*gradAlphaSigmaScalar(j-1)+gradKPrev*gradKSigmaScalar(j-1);     
    gradKMuScalar(j)=gradKPrev*gradKMuScalar(j-1);
    gradSigmaBeta=-2*beta^-3;
    gradSigmaAlpha=-alpha(j)^-2*factor^-2;
    gradSigmaGamma=-2*(alpha(j)^-1+coeff)*factor^-3-2*coeff*factor^-1*(1-factor^2)^-1;
    gradSigmaSigmaScalar=(1-factor^2)*factor^-2;
    preGradBeta(i)=preGradBeta(i)+gradientSigma*(gradSigmaBeta+gradSigmaAlpha*gradAlphaBeta(j))+gradK*gradKBeta(j);
    
    preGradSigmaScalar(i)=preGradSigmaScalar(i)+...
        gradientSigma*(gradSigmaSigmaScalar+gradSigmaAlpha*gradAlphaSigmaScalar(j))+gradK*gradKSigmaScalar(j);
    preGradGamma=preGradGamma+gradientSigma*(gradSigmaGamma+gradSigmaAlpha*gradAlphaGamma(j))+gradK*gradKGamma(j);
    preGradMuScalar(i)=preGradMuScalar(i)+gradK*gradKMuScalar(j);
  end
  gradAlphaPrev=factor^2*(alpha(end-1)^-1+coeff)^-2*alpha(end-1)^-2;  
  gradAlphaBeta(end)=2*beta+gradAlphaPrev*gradAlphaBeta(end-1);
  gradAlphaSigmaScalar(end)=-factor^2*(1-factor^2)*(alpha(end-1)^-1+coeff)^-2+gradAlphaPrev*gradAlphaSigmaScalar(end-1);
  gradAlphaGamma(end)=2*factor^3*coeff*(1-factor^2)^-1*(alpha(end-1)^-1+coeff)^-2+gradAlphaPrev*gradAlphaGamma(end-1)+2*factor*(alpha(end-1)^-1+coeff)^-1;
  alpha(end)=(beta^2+factor^2*(alpha(end-1)^-1+coeff)^-1);
  k(end)=alpha(end)^-1*(data(i,1)*beta^2+factor^2*k(end-1)*(alpha(end-1)^-1+coeff)^-1);
  sigma=X(i,:)*Sigma*X(i,:)'+alpha(end)^-1;  
  gradientSigma=0.5*(sigma^-1-(k(end)-number)^2/sigma^2);
  gradK=(k(end)-number)/sigma;
  gradKPrev=alpha(end)^-1*factor^2*(alpha(end-1)^-1+coeff)^-1;
  gradKAlpha=-k(end)/alpha(end);
  gradKAlphaPrev=factor^2*alpha(end)^-1*k(end-1)*alpha(end-1)^-2*(alpha(end-1)^-1+coeff)^-2;
  gradKBeta(end)=2*alpha(end)^-1*data(i,1)*beta+gradKPrev*gradKBeta(end-1)+gradKAlpha*gradAlphaBeta(end)+gradKAlphaPrev*gradAlphaBeta(end-1);
  gradKGamma(end)=2*factor*alpha(end)^-1*(alpha(end-1)^-1+coeff)^-1*k(end-1)*(1+factor^2*coeff*(alpha(end-1)^-1+coeff)^-1*(1-factor^2)^-1)+...
        gradKAlpha*gradAlphaGamma(end)+gradKAlphaPrev*gradAlphaGamma(end-1)+gradKPrev*gradKGamma(end-1);

  gradKSigmaScalar(end)=-factor^2*alpha(end)^-1*k(end-1)*(alpha(end-1)^-1+coeff)^-2*(1-factor^2)+...
        gradKAlpha*gradAlphaSigmaScalar(end)+gradKAlphaPrev*gradAlphaSigmaScalar(end-1)+gradKPrev*gradKSigmaScalar(end-1);     
  gradKMuScalar(end)=gradKPrev*gradKMuScalar(end-1);
  gradSigmaBeta=0;
  gradSigmaAlpha=-alpha(end)^-2;
  gradSigmaGamma=0;
  gradSigmaSigmaScalar=1;
  preGradBeta(i)=preGradBeta(i)+gradientSigma*(gradSigmaAlpha*gradAlphaBeta(end))+gradK*gradKBeta(end);
    
  preGradSigmaScalar(i)=preGradSigmaScalar(i)+...
        gradientSigma*(gradSigmaSigmaScalar+gradSigmaAlpha*gradAlphaSigmaScalar(end))+gradK*gradKSigmaScalar(end);
  preGradGamma=preGradGamma+gradientSigma*(gradSigmaAlpha*gradAlphaGamma(end))+gradK*gradKGamma(end);
  preGradMuScalar(i)=preGradMuScalar(i)+gradK*gradKMuScalar(end)-(k(end)-number)/sigma;
end
preGradSigma=X'*diag(preGradSigmaScalar)*X;
postGradSigma=2*preSigma*preGradSigma;
gradSigma=zeros(nEffectGenes,nTrans);
for i=1:size(R,1)
    gradSigma(R(i),C(i))=postGradSigma(R(i),C(i));
end
[row,col,gradSigma]=find(gradSigma);

gradDiagonal=2*diagonal.*diag(preGradSigma)'; 
gradMu=(preGradMuScalar*X);
gradBeta=sum(preGradBeta);
gradGamma=-preGradGamma*sin(2*gamma);

f=[gradBeta,gradGamma,gradMu,gradSigma', gradDiagonal];
  
