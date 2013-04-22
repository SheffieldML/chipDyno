function list=chipDynoPerGenes(anno)

% CHIPDYNOPERGENES genes periodic according to Tu et al
%
%	Description:
%	list=chipDynoPerGenes(anno)
%% 	chipDynoPerGenes.m version 1.4


TuPerGenes=textread('./data/MetabolData/PerGenes.txt', ...
                     '%q');
index=[];
for i=1:size(anno,1)
  index=[index,size(find(strcmp(anno(i),TuPerGenes)),1)];
end
list=anno(find(index));