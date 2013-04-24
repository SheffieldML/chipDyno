% DEMSPELLMANP1P2 demonstrates dynamical chipCHIP on Spellman data.
% CHIPDYNO toolbox
% demSpellmanp1p2.m version 1.2
% FORMAT 
% DESC demonstrates dynamical chipCHIP on Spellman data.
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : demTu

clear all
randn('seed',39)
[data,X,probeNames,annotation,TransNames]=chipDynoLoadDatap2();
[data,X,probeNames,annotation,TransNames]=chipDynoLoadDatap1p2(TransNames,probeNames);
nGenes=size(data,1);
npts=size(data,2);
nTrans=size(X,2);
options=foptions;
options(1)=1;
options(14)=10000;
g=cov(data');
muIn=zeros(nTrans,1);

[R,C,V,nEffectGenes]=chipReduceVariables(X);
diagonal=(randn(1,nTrans)).^2;
beta=3;
gamma=pi/4;
params=[beta,gamma,muIn',randn(1,size(V,1)).*V', diagonal];
params=scg('chipDynoLikeStat',params,options,'chipDynoLikeStatGrad',data, ...
           X, nEffectGenes,R, C);
V=params(nTrans+3:end-nTrans)';
preSigma=sparse(R,C,V,nEffectGenes,nTrans);
diagonal=params(end-nTrans+1:end);
Sigma=preSigma'*preSigma+diag(diagonal.*diagonal);
beta=params(1);
gamma=params(2);
mu=params(3:2+nTrans);
save ResultsSpellmanp1p2 params
