# CHIPDYNOLOADDATA loads Spellman Data with Lee et al ChIP data.
# CHIPDYNO toolbox
# chipDynoLoadData.R version 1.0.1
# FORMAT chipDynoLoadData <- function()
# DESC loads Spellman Data with Lee et al ChIP data.
# RETURN val : concatenated dataframe of (1) point estimate of the expression level,
# (2) connectivity measurement between genes and transcription factors
# (3), gene names and (4) transcription factors
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : 

chipDynoLoadData <- function() {

data_file="./data/SpellmanMicro.txt"
source("chipTextRead.R")
probeName_data = chipTextRead(data_file);

probeName = probeName_data[[1]] 
data = probeName_data [[2]]		

geneName_annotation_file = "./data/Yeast_Connectivity.txt"
connectivity_file="./data/Connectivity_Matrix.txt"
source("chipChipTextRead.R")
probeName2_annotation_dataChip = chipChipTextRead(geneName_annotation_file, connectivity_file);
probeName2 = probeName2_annotation_dataChip[[1]]
annotation = probeName2_annotation_dataChip[[2]]
dataChip = probeName2_annotation_dataChip[[3]]

transcription_file="./data/Trans_Names.txt"

TransNames <- read.table(transcription_file, header = TRUE, sep = "\t", quote = "", dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

TransNames=TransNames[4:(length(TransNames)-2)]
TransNames = matrix(t(TransNames),,1)
index = array(0,dim <- c(nrow(dataChip),1))
 
for (i in 1: nrow(dataChip)){
	vec <- probeName2[i]==probeName
	index[i]=colSums(vec)
}

dataChip=dataChip[which(index!=0),]
annotation=annotation[which(index!=0),]

index=array(0, dim<-c(nrow(data),1))

for (i in 1: nrow(data)){
	vec <- probeName[i]==probeName2
	index[i]=colSums(vec)
}
data=data[which(index!=0),];
probeName=probeName[which(index!=0),];

X=array(0, dim <-c(nrow(dataChip),ncol(dataChip)));
I <- c(which(dataChip<1e-3))
X[I]=1;

fakeX = rowSums(X)
X=X[which(fakeX!=0),]
annotation = matrix(annotation,,1)
annotation=annotation[which(fakeX!=0),]

effectX=colSums(X)
TransNames=TransNames[which(effectX!=0),]
TransNames=matrix(TransNames,,1)
X=X[,which(effectX!=0)]
data= data[which(fakeX!=0),]

val=list(data,X,annotation,TransNames)

return(val)
}