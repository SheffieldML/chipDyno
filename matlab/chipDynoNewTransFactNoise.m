function list=chipDynoNewTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         TransNames, annotation);

% CHIPDYNONEWTRANSFACTNOISE tfs active for us and not for Tu et al.
%
%	Description:
%	list=chipDynoNewTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
%                                         TransNames, annotation);
%% 	chipDynoNewTransFactNoise.m version 1.4


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
  