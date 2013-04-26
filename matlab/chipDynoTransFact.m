function [TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, name);

% CHIPDYNOTRANSFACT provides gene-specific TFAs with errorbars.
% CHIPDYNO toolbox
% chipDynoTransFact.m version 1.4
% FORMAT [TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma, ...
%                                         mu, transNames, annotations, name);
% DESC given a transcription factor, provides gene-specific TFAs with errorbars.
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix of TFA
% ARG beta :
% ARG gamma : degree of temporal continuity
% ARG mu : mean value of the transcription factor activity
% ARG transNames : Transcription factors
% ARG annotations : Gene names
% ARG name : given transcription factor name
% RETURN TF : transcription factor activity
% RETURN TFError : error in transcription factor activity
% RETURN TFErrorDiff :
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoTransFactNoise

index=find(strcmp(name,transNames));
genesIn=find(X(:,index));
anno=annotations(find(X(:,index)));
nTargets=size(anno,1);
npts=size(data,2);
TF=zeros(nTargets,npts);
TFError=zeros(nTargets,npts);
TFErrorDiff=zeros(npts,npts,nTargets);
for i=1:nTargets
  [TF(i,:),TFError(i,:),TFErrorDiff(:,:,i)]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                         name,genesIn(i));
end
