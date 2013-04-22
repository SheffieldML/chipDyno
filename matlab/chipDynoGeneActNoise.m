function [list,maxActivity,maxActivityError]=chipDynoGeneActNoise(data, ...
                                                  X,Sigma,beta,precs,gamma,mu, ...
                                                  TransNames, ...
                                                  annotation,geneName);

% CHIPDYNOGENEACTNOISE considering uncertainty of the expression level,  
% lists activators in decreasing order for a given gene
% CHIPDYNO toolbox
% chipDynoGeneActNoise.m version 1.4
% FORMAT [list,maxActivity,maxActivityError]=chipDynoExpectationsFastNoise(data,X,Sigma, ...
%					beta, precs, gamma, mu, TransNames, ...
%					annotation, geneName);
% DESC considering uncertainty of the expression level,  
% lists activators in decreasing order for a given gene
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix of TFA
% ARG beta :
% ARG precs : uncertainty of the expression level
% ARG gamma : degree of temporal continuity
% ARG mu : mean value of the transcription factor activity
% ARG TransNames : Transcription factors
% ARG annotation : Gene names
% ARG geneName : specific gene name
% RETURN list : lists activators for a given gene in decreasing order
% RETURN maxActivity : 
% RETURN maxActivityError : 
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoGeneAct, chipDynoActTransFact, chipDynoActTransFactNoise


I=find(strcmp(geneName,annotation));
activeNames=TransNames(find(X(I,:)));
nTransFact=sum(X(I,:));
maxActivity=[];
maxActivityError=[];
for i=1:nTransFact
  [tf,tfError]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation, ...
                                         activeNames(i),geneName);
  [act,index]=max(tf-mu(find(strcmp(activeNames(i),TransNames)))*ones(size(tf)));
  maxActivity=[maxActivity,act];
  maxActivityError=[maxActivityError,tfError(index)];
end
[maxActivity,index]=sort(maxActivity,'descend');
maxActivityError=maxActivityError(index);
list=activeNames(index);
