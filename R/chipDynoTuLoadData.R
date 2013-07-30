# CHIPDYNOTULOADDATA loads Tu Data with Lee et al ChIP data.
#
#	Description:
#	[data,vars,X,annotation,TransNames]=chipDynoTuLoadData();
## 	rChipDynoTuLoadData.R version 0.01
#	Written on 09.11.2011

chipDynoTuLoadData <- function() {

file_dictionary <- "./data/MetabolData/dictionary.txt";
file_probeIDTu <- "./data/MetabolData/probeIDTu.txt";
file_data <- "./data/MetabolData/YeastMetabolism_exprs.txt";
file_vars <- "./data/MetabolData/YeastMetabolism_se.txt";

file_dataChip <- "./data/Connectivity2.txt";
file_annotation <- "./data/annotations2.txt"
file_transNames <- "./data/Trans_Names2.txt"

#file_dictionary, file_probeIDTu, file_data, file_vars

source("chipTuTextRead.R") # This will collect the values of variable "ORF", "data", "vars" form "rChipTuTextRead.R" file.
ORF_data_vars = chipTuTextRead(file_dictionary, file_probeIDTu, file_data, file_vars);

ORF = ORF_data_vars[[1]]
data = ORF_data_vars[[2]]
vars = ORF_data_vars[[3]]

# data = as.matrix(read.table("./data/MetabolData/YeastMetabolism_exprs.txt")) #9335X36 matrix
# vars=as.matrix(read.table("./data/MetabolData/YeastMetabolism_se.txt")) #9335X36 matrix

zeroValueRow = which(rowSums(vars)==0)

data = data[-zeroValueRow, ] #9307X36 matrix
vars = vars[-zeroValueRow, ] #9307X36 matrix

probeName = ORF # Rename ORF
probeName = matrix(probeName[-zeroValueRow, ],,1) #9307X1 matrix


dataChip=as.matrix(read.table(file_dataChip)) #6229X204 matrix

probe_anno <- read.table(file_annotation, header = FALSE, sep = "\t", quote = "", dec = ".", col.names=c("prob", "anno"), na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = FALSE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

probeName2=matrix(probe_anno$prob,,1) 

annota=matrix(probe_anno$anno,,1)

TransNames_tab <- read.table(file_transNames, header = FALSE, sep = "\t", quote = "", dec = ".", col.names=c("TN"), na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = FALSE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

TransNames = matrix(TransNames_tab$TN,,1)

noOfprobe=nrow(probeName)
redundancy=matrix(mat.or.vec(noOfprobe,1),,1) # create nX1 matrix of 0's
for (i in 1:noOfprobe){
		redundancy[i,1]=1}

index=matrix(mat.or.vec(nrow(dataChip),1),,1)

####### the following block eleminate the redundancy

for (i in 1: nrow(probeName2)){
	vec <- probeName==probeName2[i,1]
	index[i]=colSums(vec)
	if(index[i]>1){
		pippo=(which(vec))
		redundancy[pippo[2:length(pippo)]]=0
		}
	}
########### end of the redundancy elemination block


#dataChip=dataChip(find(index),:);
#annota=annota(find(index));
#probeName2=probeName2(find(index));

dataChip=dataChip[which(index!=0),]
annota=annota[which(index!=0),]
probeName2=probeName2[which(index!=0),]
probeName2=matrix(probeName2,,1)

#probeName=probeName(find(redundancy));
#data=data(find(redundancy),:);
#vars=vars(find(redundancy),:);

probeName=probeName[which(redundancy!=0),]
probeName=matrix(probeName,,1)

data=data[which(redundancy!=0),]
vars=vars[which(redundancy!=0),]

###########
# index=zeros(size(data,1),1);
# preX=[];
# annotation=[];
# for i=1:size(data,1)
#    index(i)=sum(strcmp(probeName(i),probeName2));
#    if index(i)
#      preX=[preX;dataChip(find(strcmp(probeName(i), probeName2)),:)];
#      annotation=[annotation; annota(find(strcmp(probeName(i),probeName2)))];
#    end
# end

#preX=mat.or.vec(0,0)

preX=NULL
annotation=NULL
index=mat.or.vec(nrow(data),1)
for (i in 1: nrow(data)){
	index[i]=sum(probeName[i]==probeName2)
	if (index[i]==1)
		
		preX=rbind(preX,dataChip[which(probeName[i]==probeName2),])
		annotation=rbind(annotation,annota[which(probeName[i]==probeName2)])
}

############## Started again 14.12.2011

#data=data(find(index),:);
#vars=vars(find(index),:);
#probeName=probeName(find(index));
#X=zeros(size(preX,1),size(preX,2));
#I=find(preX<1e-3);
#X(I)=1;


data= data[which(index==1),]
vars= vars[which(index==1),]
probeName= probeName[which(index==1),]
probeName=matrix(probeName,,1)
X= mat.or.vec(nrow(preX),ncol(preX))
I <- c(which(preX<1e-3))
X[I] =1

#fakeX=sum(X,2);
#X=X(find(fakeX),:);
#annotation=annotation(find(fakeX));
#effectX=sum(X,1);
#TransNames=TransNames(find(effectX));
#X=X(:,find(effectX));
#data=data(find(fakeX),:);
#vars=vars(find(fakeX),:);

fakeX = rowSums(X)
X=X[which(fakeX!=0),]
annotation=annotation[which(fakeX!=0),]
annotation = matrix(annotation,,1)

effectX=colSums(X)
TransNames=TransNames[which(effectX!=0),]
TransNames=matrix(TransNames,,1)
X=X[,which(effectX!=0)]
data= data[which(fakeX!=0),]
vars= vars[which(fakeX!=0),]

data_vars_X_annotation_TransNames = list(data, vars, X, annotation, TransNames)

return(data_vars_X_annotation_TransNames)
}

