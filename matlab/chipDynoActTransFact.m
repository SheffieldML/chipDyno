function [list,newX, newXVals]=chipDynoActTransFact(data,X,Sigma,beta,gamma,mu, ...
                                         TransNames, annotation,sigLev);

% CHIPDYNOACTTRANSFACT identifies significantly varying TFs.
% CHIPDYNO toolbox
% chipDynoActTransFact.m version 1.1
% FORMAT [list,newX, newXVals]=chipDynoActTransFact(data,X,Sigma,beta,gamma,mu, ...
%                                         TransNames, annotation,sigLev);
% DESC identifies significantly varying TFs.
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix
% ARG beta :
% ARG gamma :
% ARG mu : mean value of the transcription factor activity
% ARG TransNames : Transcription factors
% ARG annotation : Gene names
% ARG sigLev : threshold value
% RETURN list: list of regulators for a specific gene
% RETURN newX: 
% RETURN newXVals : 
% COPYRIGHT : Neil D. Lawrence, 2005
% COPYRIGHT : Guido Sanguinetti, 2005
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoTransFact, chipDynoTransFactNoise, chipDynoActTransFactNoise


nTrans=size(TransNames,1);
list=[];
newX=zeros(size(X));
newXVals=zeros(size(X));
for i=1:nTrans
    [TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, ...
                                         TransNames, annotation, ...
                                        TransNames(i));
    %vars=max(abs((TF-mu(i)*ones(size(TF)))'./TFError'));
    maxVars=chipDynoMaxDiff(TF,TFErrorDiff);
    sigVars=maxVars(find(maxVars>sigLev));
    list=[list, size(sigVars,2)];
    index=find(X(:,i));
    newX(index(find(maxVars>sigLev)),i)=1;
    newXVals(index(find(maxVars>sigLev)),i)=maxVars((find(maxVars>sigLev)));
end
