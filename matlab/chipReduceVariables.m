function [R,C,V,nEffectGenes]=chipReduceVariables(X);

% CHIPREDUCEVARIABLES reduce  number of variables in chipDyno model
%
%	Description:
%	[R,C,V,nEffectGenes]=chipReduceVariables(X);
%% 	chipReduceVariables.m version 1.4


nGenes=size(X,1);
preSigma1=X(1,:);
preSigma2=preSigma1;

for i=2:nGenes

  preSigma2=[preSigma2;X(i,:)];
  matrix1=preSigma1'*preSigma1;
  matrix2=preSigma2'*preSigma2;         
  decider=min(matrix1(find(matrix2)));
  if decider==0
    preSigma1=[preSigma1;X(i,:)];
  end
end
nEffectGenes=size(preSigma1,1);
[R,C,V]=find(preSigma1);
