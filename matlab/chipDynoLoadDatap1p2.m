function [data,X,probeNames,annotation,TransNames] = chipDynoLoadDatap1p2(TransNamesp2, probeNamesp2);

% CHIPDYNOLOADDATAP1P2 loads Spellman Data with Lee et al ChIP data.
%
%	Description:
%	[data,X,probeNames,annotation,TransNames] = chipDynoLoadDatap1p2(TransNamesp2, probeNamesp2);
%% 	chipDynoLoadDatap1p2.m version 1.3


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
I=find(dataChip<1e-3);
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
newX=zeros(size(probeNamesp2,1),size(TransNamesp2,1));
newData=zeros(size(probeNamesp2,1),size(data,2));
newAnnotation=[];
for i=1:size(newX,1)
  for j=1:size(newX,2)
    newX(i,j)=X(find(strcmp(probeNamesp2(i),probeNames)), ...
                find(strcmp(TransNamesp2(j),TransNames)));
  end
  newData(i,:)=data(find(strcmp(probeNamesp2(i),probeNames)),:);
  newAnnotation=[newAnnotation;annotation(find(strcmp(probeNamesp2(i),probeNames)))];
end
X=newX;
probeNames=probeNamesp2;
TransNames=TransNamesp2;
data=newData;
annotation=newAnnotation;
