function [TF,TFError,TFErrorDiff]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         transNames, annotations, ...
                                        name);

% CHIPDYNOTRANSFACTNOISE given a transcription factor, provides TFAs.
%
%	Description:
%	[TF,TFError,TFErrorDiff]=chipDynoTransFactNoise(data,X,Sigma,beta,precs,gamma,mu, ...
%                                         transNames, annotations, ...
%                                        name);
%% 	chipDynoTransFactNoise.m version 1.4


index=find(strcmp(name,transNames));
genesIn=find(X(:,index));
anno=annotations(find(X(:,index)));
nTargets=size(anno,1);
npts=size(data,2);
TF=zeros(nTargets,npts);
TFError=zeros(nTargets,npts);
TFErrorDiff=zeros(npts,npts,nTargets);
for i=1:nTargets
  [TF(i,:),TFError(i,:),TFErrorDiff(:,:,i)]=chipDynoExpectationsFastNoise(data,X,Sigma,beta,precs,gamma,mu, ...
                                         transNames, annotations, ...
                                         name,genesIn(i));
end
