% DEMTU demonstrates dynamical chipCHIP on Tu data.
% CHIPDYNO toolbox
% demTu.m version 1.4
% FORMAT 
% DESC demonstrates dynamical chipCHIP on Tu data.
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : demSpellmanDynoStat

clear all
#[data,vars,X,annotation,TransNames]=chipDynoTuLoadData();
[data,vars,X,annotations,transNames]=chipDynoTuLoadData();
nGenes=size(data,1);
npts=size(data,2);
nTrans=size(X,2);
options=foptions;

options(1)=1;
options(14)=50000;
muIn=zeros(nTrans,1);
[R,C,V,nEffectGenes]=chipReduceVariables(X);
diagonal=0.5*ones(1,nTrans);
precs=ones(size(vars,1),size(vars,2))./(vars.^2);
beta=3;
gamma=pi/4;
params=[beta,gamma,muIn',0.1*V', diagonal];
params=scg('chipDynoLikeStatNoise',params,options,'chipDynoLikeStatNoiseGrad',data, precs, ...
           X, nEffectGenes,R, C);
V=params(nTrans+3:end-nTrans)';
preSigma=sparse(R,C,V,nEffectGenes,nTrans);
diagonal=params(end-nTrans+1:end);
Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
beta=params(1);
gamma=params(2);
mu=params(3:2+nTrans);
save results/ResultsTu params
