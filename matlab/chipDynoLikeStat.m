function likelihood=chipDynoLikeStat(params,data,X,nEffectGenes,R,C);

% CHIPDYNOLIKESTAT marginal likelihood for chipChip dynamical model
% CHIPDYNO toolbox
% chipDynoLikeStat.m version 1.4
% FORMAT likelihood=chipDynoLikeStat(params,data,X,nEffectGenes,R,C);
% DESC compute the marginal likelihood for chipChip dynamical model
% ARG params: concatenated vector of multiple parameters(beta, gamma, 
% initial mean of the transcription factors, and 
% a vector to create diagonal matrix used to reduce the sparsity of covariance)
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG nEffectGenes : effectice gene name
% ARG R, C : same length integer vectors specifying the row and column 
% indices of the non-zero entries of the sparce matrix
% RETURN likelihood : marginal likelihood
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoLikeStatNoise

nGenes=size(data,1);
nTrans=size(X,2);
npts=size(data,2);
beta=params(1);
gamma=params(2);
mu=params(3:2+nTrans)';
V=params(3+nTrans:end-nTrans)';
preSigma=sparse(R,C,V,nEffectGenes,nTrans);
diagonal=params(end-nTrans+1:end);
Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
factor=cos(gamma)^2;
preLike=zeros(1, nGenes);
k=zeros(1,npts);
alpha=zeros(1,npts);
for i=1:nGenes
  coeff=(1-factor^2)*X(i,:)*Sigma*X(i,:)';
  number=X(i,:)*mu;
  alpha(1)=beta^2;
  k(1)=(data(i,npts)-(1-factor)*number)/factor;
  sigma=beta^-2+(alpha(1)^-1+coeff)/factor^2;
  preLike(i)=preLike(i)+ 0.5*log(sigma)+...
       0.5*((k(1)-data(i,npts-1))^2)/(sigma);
  for j=2:npts-1
    alpha(j)=(beta^2+factor^2*(alpha(j-1)^-1+coeff)^-1);
    k(j)=alpha(j)^-1*(data(i,npts-j+1)*beta^2+factor^2*k(j-1)*(alpha(j-1)^-1+coeff)^-1);
    sigma=beta^-2+(alpha(j)^-1+coeff)/factor^2;
    preLike(i)=preLike(i)+ 0.5*log(sigma)+...
      0.5*((k(j)-data(i,npts-j))^2)/(sigma);
  end
  alpha(end)=(beta^2+factor^2*(alpha(end-1)^-1+coeff)^-1);
  k(end)=alpha(end)^-1*(data(i,1)*beta^2+factor^2*k(end-1)*(alpha(end-1)^-1+coeff)^-1);
  preLike(i)=preLike(i)+0.5*log(alpha(end)^-1+X(i,:)*Sigma*X(i,:)')+...
      0.5*((k(end)-number)^2)/(alpha(end)^-1+X(i,:)*Sigma*X(i,:)');
end
likelihood=sum(preLike);
