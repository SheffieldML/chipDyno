function [TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, ...
                                         transNames, annotations, ...
                                        name);

% CHIPDYNOTRANSFACT provides gene-specific TFAs with errorbars.
%
%	Description:
%	[TF,TFError,TFErrorDiff]=chipDynoTransFact(data,X,Sigma,beta,gamma,mu, ...
%                                         transNames, annotations, ...
%                                        name);
%% 	chipDynoTransFact.m version 1.4


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
