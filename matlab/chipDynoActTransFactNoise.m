function [list,newX, newXVals]=chipDynoActTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation);

% CHIPDYNOACTTRANSFACTNOISE identifies significantly varying TFs with uncertainty of..
% expression level.
% CHIPDYNO toolbox
% chipDynoActTransFactNoise.m version 1.4
% FORMAT [list,newX, newXVals]=chipDynoActTransFact(data,X,Sigma,beta, precs, gamma,mu, ...
%                                         TransNames, annotation,sigLev);
% DESC identifies significantly varying TFs.
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix of TFA
% ARG beta :
% ARG precs : uncertainty of the expression level 
% ARG gamma : degree of temporal continuity
% ARG mu : mean value of the transcription factor activity
% ARG TransNames : Transcription factors
% ARG annotation : Gene names
% ARG sigLev : threshold value
% RETURN list: list of regulators for a specific gene
% RETURN newX: 
% RETURN newXVals : 
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoTransFact, chipDynoTransFactNoise, chipDynoActTransFact

nTrans=size(TransNames,1);
list=[];
newX=zeros(size(X));
newXVals=zeros(size(X));
for i=1:nTrans
    [TF,TFError,TFErrorDiff]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation, ...
                                        TransNames(i));
    maxVars=chipDynoMaxDiff(TF,TFErrorDiff);
    pvals=cumGaussian(-maxVars));
    sigVars=pvals(find(pvals<0.02));
    list=[list, size(sigVars,2)];
    index=find(X(:,i));
    newX(index(find(pvals<0.02)),i)=1;
     newXVals(index(find(pvals<0.02)),i)=pvals((find(pvals<0.02)));
end
