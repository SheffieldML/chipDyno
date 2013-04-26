function list=chipDynoNewTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation);

% CHIPDYNONEWTRANSFACTNOISE transcription factors active for us and not for Tu et al.
% CHIPDYNO toolbox
% chipDynoNewTransFactNoise.m version 1.4
% FORMAT list=chipDynoNewTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
%                                         TransNames, annotation);
% DESC determins tfs active for us and not for Tu et al
% ARG data : point estimate of the expression level
% ARG X : connectivity measurement between genes and transcription factors
% ARG Sigma : prior covariance matrix of TFA
% ARG beta :
% ARG precs : uncertainty of the expression level
% ARG gamma : degree of temporal continuity
% ARG mu : mean value of the transcription factor activity
% ARG TransNames : Transcription factors
% ARG annotations : Gene names
% RETURN list: transcription factors active for us and not for Tu et al
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoTransFactNoise, chipDynoNewTransFactNoise

TuTransFact=textread('./data/MetabolData/PerTransFact.txt', ...
                     '%q');
[list1,newX,newXVals]=chipDynoActTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation);
list2=TransNames(find(list1>4));
index=[];
for i=1:size(list2,1)
  index=[index,1-size(find(strcmp(list2(i),TuTransFact)),1)];
end
list=list2(find(index));
  
