#function [data,X,annotation,TransNames]=chipDynoLoadData();

#% CHIPDYNOLOADDATA loads Spellman Data with Lee et al ChIP data.
#%
#%	Description:
#%	[data,X,annotation,TransNames]=chipDynoLoadData();
#%% 	chipDynoLoadData.m version 1.5

chipDynoLoadData = function() {

data_file="./data/SpellmanMicro.txt"
source("chipTextRead.R")
probeName_data = chipTextRead(data_file);

probeName = probeName_data[[1]] 	# Gene Names
data = probeName_data [[2]]		# Time series data

#[probeName, data] = chipTextRead('./data/SpellmanMicro.txt');

geneName_annotation_file = "./data/Yeast_Connectivity.txt"
connectivity_file="./data/Connectivity_Matrix.txt"
source("chipChipTextRead.R")
probeName2_annotation_dataChip = chipChipTextRead(geneName_annotation_file, connectivity_file);
probeName2 = probeName2_annotation_dataChip[[1]]
annotation = probeName2_annotation_dataChip[[2]]
dataChip = probeName2_annotation_dataChip[[3]]

#[probeName2, annotation, dataChip] = chipChipTextRead(['./data/Yeast_Connectivity.txt'], './data/Connectivity_Matrix.txt');

transcription_file="./data/Trans_Names.txt"

TransNames <- read.table(transcription_file, header = TRUE, sep = "\t", quote = "", dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

########??? Start  
#TransNames=textread('./data/Trans_Names.txt','%q', ...
#                    'headerlines',1,'whitespace','','delimiter','\t');

TransNames=TransNames[4:(length(TransNames)-2)]
TransNames = matrix(t(TransNames),,1)
index = array(0,dim <- c(nrow(dataChip),1))
#TransNames=TransNames(4:end-2);
#index=zeros(size(dataChip,1),1);

#probeName=matrix(probeName,,1)
#probeName2=matrix(probeName2[2:length(probeName2)],,1)
#annotation =matrix (annotation,,1)
 
for (i in 1: nrow(dataChip)){
	vec <- probeName2[i]==probeName
	index[i]=colSums(vec)
}

#for i=1:size(dataChip,1)
#    index(i)=sum(strcmp(probeName2(i),probeName));
#end

dataChip=dataChip[which(index!=0),]
annotation=annotation[which(index!=0),]

#dataChip=dataChip(find(index),:);
#annotation=annotation(find(index));

index=array(0, dim<-c(nrow(data),1))

for (i in 1: nrow(data)){
	vec <- probeName[i]==probeName2
	index[i]=colSums(vec)
}
data=data[which(index!=0),];
probeName=probeName[which(index!=0),];

#index=zeros(size(data,1),1);
#for i=1:size(data,1)
#    index(i)=sum(strcmp(probeName(i),probeName2));
#end

#data=data(find(index),:);
#probeName=probeName(find(index));

X=array(0, dim <-c(nrow(dataChip),ncol(dataChip)));
I <- c(which(dataChip<1e-3))
X[I]=1;

#I=find(dataChip<1e-3);
#X(I)=1;
#X=zeros(size(dataChip,1),size(dataChip,2));

fakeX = rowSums(X)
X=X[which(fakeX!=0),]
annotation = matrix(annotation,,1)
annotation=annotation[which(fakeX!=0),]

effectX=colSums(X)
TransNames=TransNames[which(effectX!=0),]
TransNames=matrix(TransNames,,1)
X=X[,which(effectX!=0)]
data= data[which(fakeX!=0),]

#fakeX=sum(X,2);
#X=X(find(fakeX),:);
#annotation=annotation(find(fakeX));
#effectX=sum(X,1);
#TransNames=TransNames(find(effectX));
#X=X(:,find(effectX));
#data=data(find(fakeX),:);
val=list(data,X,annotation,TransNames)

return(val)
}
