function [data,X,probeNames,annotation,TransNames]=chipDynoLoadDatap5();

% CHIPDYNOLOADDATAP5 loads Spellman Data with Lee et al ChIP data.
%
%	Description:
%	[data,X,probeNames,annotation,TransNames]=chipDynoLoadDatap5();
%% 	chipDynoLoadDatap5.m version 1.3

  
[probeName, data] = chipTextRead(['./data/' ...
                    'SpellmanMicro.txt']);
[row,col,how]=find(data==0);
fakeData=sparse(row,col,how,size(data,1),size(data,2));
data=data(find(sum(fakeData,2)<5),:);
probeName=probeName(find(sum(fakeData,2)<5));
[probeName2, annotation, dataChip] = chipChipTextRead(['./data/' ...
                    'Yeast_Connectivity.txt'], './data/Connectivity_Matrix.txt');
TransNames=textread('./data/Trans_Names.txt','%q', ...
                    'headerlines',1,'whitespace','','delimiter','\t');
TransNames=TransNames(4:end-2);
index=zeros(size(dataChip,1),1);
for i=1:size(dataChip,1)
    index(i)=sum(strcmp(probeName2(i),probeName));
end
dataChip=dataChip(find(index),:);
annotation=annotation(find(index));
index=zeros(size(data,1),1);
for i=1:size(data,1)
    index(i)=sum(strcmp(probeName(i),probeName2));
end
data=data(find(index),:);
probeName=probeName(find(index));
X=zeros(size(dataChip,1),size(dataChip,2));
I=find(dataChip<5e-3);
X(I)=1;
%X=X(:,7:25);


fakeX=sum(X,2);
X=X(find(fakeX),:);
annotation=annotation(find(fakeX));
effectX=sum(X,1);
TransNames=TransNames(find(effectX));
X=X(:,find(effectX));
data=data(find(fakeX),:);
probeNames=probeName(find(fakeX));
