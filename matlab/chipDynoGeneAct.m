function [list,maxActivity,maxActivityError]=chipDynoGeneAct(data, ...
                                                  X,Sigma,beta,gamma,mu, ...
                                                  transNames, ...
                                                  annotation,geneName);

% CHIPDYNOGENEACT given a gene, lists activators in decreasing order
% CHIPDYNO toolbox
% chipDynoGeneAct.m version 1.4
% FORMAT [list,maxActivity,maxActivityError]=chipDynoExpectationsFastNoise(data,X,Sigma, ...
%					beta, gamma, mu, transNames, ...
%					annotation, geneName);
% DESC given a gene, lists activators in decreasing order
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix of TFA
% ARG beta :
% ARG gamma : degree of temporal continuity
% ARG mu : mean value of the transcription factor activity
% ARG transNames : Transcription factors
% ARG annotation : Gene names
% ARG geneName : specific gene name
% RETURN list : lists activators for a given gene in decreasing order
% RETURN maxActivity : 
% RETURN maxActivityError : 
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoGeneActNoise, chipDynoActTransFact, chipDynoActTransFactNoise

I=find(strcmp(geneName,annotation));
activeNames=transNames(find(X(I,:)));
nTransFact=sum(X(I,:));
maxActivity=[];
maxActivityError=[];
for i=1:nTransFact
  [tf,tfError]=chipDynoExpectationsFast(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotation, ...
                                         activeNames(i),geneName);
  ind=find(strcmp(activeNames(i),transNames));
  tf=tf-mu(ind)*ones(size(tf));
  [act,index]=max(tf);
  maxActivity=[maxActivity,act];
  maxActivityError=[maxActivityError,tfError(index)];
end
[maxActivity,index]=sort(maxActivity,'descend');
maxActivityError=maxActivityError(index);
list=activeNames(index);
