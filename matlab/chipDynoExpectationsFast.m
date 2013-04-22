function [tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, transName,geneName);

% CHIPDYNOEXPECTATIONSFAST computes posterior expectations of TFA.
% CHIPDYNO toolbox
% chipDynoExpectationsFast.m version 1.5
% FORMAT [tf,tfErrors,tfErrorsDiffs]=chipDynoExpectationsFast(data,X,Sigma,beta, gamma, ...
%                                      mu, transNames, annotations, transName,geneName);
% DESC computes posterior expectations of TFA.
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix of TFA
% ARG beta :
% ARG gamma : degree of temporal continuity
% ARG mu : mean value of the transcription factor activity
% ARG transNames : Transcription factors
% ARG annotations : Gene names
% ARG transName : specific transcription factor
% ARG geneName : specific gene name
% RETURN tf: gene specific transcription factor activity
% RETURN tfErrors: error in gene specific transcription factor activity
% RETURN tfErrorsDiffs : 
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoExpectationsFastNoise


npts=size(data,2);
nTrans=size(X,2);

c=class(geneName);
if c(1)=='c'
    x=X(find(strcmp(geneName,annotations)),:)';
    data=data(find(strcmp(geneName,annotations)),:);
elseif c(1)=='d'
    x=X(geneName,:)';
    data=data(geneName,:);
else
    error('Genes can be identified either by number or name\n')
end
expectations=chipDynoStatPostEst(data,x,Sigma,beta,gamma,mu);
index=find(strcmp(transName,transNames));
if x(index)==0
  error('The gene selected is not a target of the transcription factor \n')
end
tf=expectations.b(:,index);
ind=find(strcmp(transName,transNames(find(x))));
tfErrors=expectations.tfError(:,ind);
tfErrorsDiffs=[expectations.tfErrorDiffs(:,:,ind)];
