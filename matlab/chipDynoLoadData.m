function [data,X,annotation,TransNames]=chipDynoLoadData();

% CHIPDYNOLOADDATA loads Spellman Data with Lee et al ChIP data.
% CHIPDYNO toolbox
% chipDynoLoadData.m version 1.5
% FORMAT [data,X,annotation,TransNames]=chipDynoLoadData();
% DESC loads Spellman Data with Lee et al ChIP data.
% RETURN data : point estimate of the expression level
% RETURN X : connectivity measurement between genes and transcription factors
% RETURN annotation : gene names
% RETURN TransNames : transcription factors
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : 

[probeName, data] = chipTextRead('./data/SpellmanMicro.txt');
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
I=find(dataChip<1e-3);
X(I)=1;

fakeX=sum(X,2);
X=X(find(fakeX),:);
annotation=annotation(find(fakeX));
effectX=sum(X,1);
TransNames=TransNames(find(effectX));
X=X(:,find(effectX));
data=data(find(fakeX),:);
