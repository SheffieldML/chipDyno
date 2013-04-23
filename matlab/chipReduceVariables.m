function [R,C,V,nEffectGenes]=chipReduceVariables(X);

% CHIPREDUCEVARIABLES reduce  number of variables in chipDyno model
% CHIPDYNO toolbox
% chipReduceVariables.m version 1.4
% FORMAT [R,C,V,nEffectGenes]=chipReduceVariables(X);
% DESC reduce  number of variables in chipDyno model
% ARG X : connectivity measurement between genes and transcription factors
% RETURN R, C : same length integer vectors specifying the row and column 
% indices of the non-zero entries of the sparce matrix
% RETURN V : same length integer vectors with R and C specifying the
% values of the non-zero entries of the sparce matrix
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : 

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
