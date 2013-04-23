function [data,vars,X,annotation,TransNames]=chipDynoTuLoadData();

% CHIPDYNOTULOADDATA loads Tu Data with Lee et al ChIP data.
% CHIPDYNO toolbox
% chipDynoTuLoadData.m version 1.5
% FORMAT [data,vars,X,annotation,TransNames]=chipDynoTuLoadData();
% DESC loads Tu Data with Lee et al ChIP data.
% RETURN data : point estimate of the expression level
% RETURN vars : uncertainty of the expression level
% RETURN X : connectivity measurement between genes and transcription factors
% RETURN annotation : Gene names
% RETURN TransNames : transcription factors
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : 

[probeName, data, vars] = chipTuTextRead;
data=data(find(sum(vars,2)),:);
probeName=probeName(find(sum(vars,2)));
vars=vars(find(sum(vars,2)),:);
dataChip=load('./data/Connectivity2.txt');
[probeName2, annota] = textread('./data/annotations2.txt','%q%q');
TransNames=textread('./data/Trans_Names2.txt','%q');
%TransNames=TransNames(2:end);
redundancy=ones(size(probeName));
index=zeros(size(dataChip,1),1);
for i=1:size(dataChip,1)
   vec=strcmp(probeName2(i),probeName);
   index(i)=sum(vec);
   if index(i)>1
     pippo=find(vec);
     redundancy(pippo(2:end))=0;
   end
   
   
end
dataChip=dataChip(find(index),:);
annota=annota(find(index));
probeName2=probeName2(find(index));
probeName=probeName(find(redundancy));
data=data(find(redundancy),:);
vars=vars(find(redundancy),:);
index=zeros(size(data,1),1);
preX=[];
annotation=[];
for i=1:size(data,1)
    index(i)=sum(strcmp(probeName(i),probeName2));
    if index(i)
      preX=[preX;dataChip(find(strcmp(probeName(i), probeName2)),:)];
      annotation=[annotation; annota(find(strcmp(probeName(i),probeName2)))];
    end
end
data=data(find(index),:);
vars=vars(find(index),:);
probeName=probeName(find(index));
X=zeros(size(preX,1),size(preX,2));
I=find(preX<1e-3);
X(I)=1;

fakeX=sum(X,2);
X=X(find(fakeX),:);
annotation=annotation(find(fakeX));
effectX=sum(X,1);
TransNames=TransNames(find(effectX));
X=X(:,find(effectX));
data=data(find(fakeX),:);
vars=vars(find(fakeX),:);
