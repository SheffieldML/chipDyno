function [tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFastNoise(data, X, Sigma, ...
					beta, precs, gamma, mu, TransNames, ...
					annotation, transName,geneName);

% CHIPDYNOEXPECTATIONSFASTNOISE computes posterior expectations TFA
% considering uncertainty of the expression level
% CHIPDYNO toolbox
% chipDynoExpectationsFastNoise.m version 1.4
% FORMAT [tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFastNoise(data,X,Sigma, ...
%					beta, precs, gamma, mu, TransNames, ...
%					annotations, transName,geneName);
% DESC computes posterior expectations of TFA considering uncertainty of the expression level.
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix of TFA
% ARG beta :
% ARG precs : uncertainty of the expression level
% ARG gamma : degree of temporal continuity
% ARG mu : mean value of the transcription factor activity
% ARG TransNames : Transcription factors
% ARG annotations : Gene names
% ARG transName : specific transcription factor
% ARG geneName : specific gene name
% RETURN tf: gene specific transcription factor activity
% RETURN tfErrors: error in gene specific transcription factor activity
% RETURN tfErrorsDiffs : 
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoExpectationsFast

npts=size(data,2);
nTrans=size(X,2);

c=class(geneName);
if c(1)=='c'
    x=X(find(strcmp(geneName,annotation)),:)';
    data=data(find(strcmp(geneName,annotation)),:);
    precs=precs(find(strcmp(geneName,annotation)),:);
elseif c(1)=='d'
    x=X(geneName,:)';
    data=data(geneName,:);
    precs=precs(geneName,:);
else
    error('Genes can be identified either by number or name\n')
end
expectations=chipDynoStatPostEstNoise(data,x,Sigma,beta,precs,gamma,mu);
index=find(strcmp(transName,TransNames));
if x(index)==0
  error('The gene selected is not a target of the transcription factor \n')
end
tf=expectations.b(:,index);
ind=find(strcmp(transName,TransNames(find(x))));
tfErrors=expectations.tfError(:,ind);
tfErrorsDiffs=[expectations.tfErrorDiffs(:,:,ind)];
