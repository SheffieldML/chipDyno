# CHIPDYNOTULOADDATA loads Tu Data with Lee et al ChIP data.
# CHIPDYNO toolbox
# chipDynoTuLoadData.R version 1.0.1
# FORMAT chipDynoTuLoadData <- function()
# DESC loads Tu Data with Lee et al ChIP data.
# ARG file_dictionary : List of gene names
# ARG file_probeIDTu : List of gesn's id
# ARG file_data : Point estimate for the expression level
# ARG file_vars : uncertainty of the expression level
# ARG file_dataChip : connectivity matrix between genes and transcription factors
# ARG file_annotation : Gene annotation
# ARG file_transNames : List of transcription factors.
# RETURN data_vars_X_annotation_TransNames : concatenated dataframe of 
# (1) point estimate of the expression level, 
# (2) uncertainty of the expression level,
# (3) connectivity measurement between genes and transcription factors
# (4) Gene names and (5) transcription factors
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : 

chipDynoTuLoadData <- function(file_dictionary, 
	file_probeIDTu, file_data, file_vars, file_dataChip, 
	file_annotation, file_transNames) {

#file_dictionary <- "./data/MetabolData/dictionary.txt";
#file_probeIDTu <- "./data/MetabolData/probeIDTu.txt";
#file_data <- "./data/MetabolData/YeastMetabolism_exprs.txt";
#file_vars <- "./data/MetabolData/YeastMetabolism_se.txt";

#file_dataChip <- "./data/Connectivity2.txt";
#file_annotation <- "./data/annotations2.txt"
#file_transNames <- "./data/Trans_Names2.txt"

source("chipTuTextRead.R")  # This will collect the values of variable
							# "ORF", "data", "vars" form "ChipTuTextRead.R" file.
ORF_data_vars = chipTuTextRead(file_dictionary, file_probeIDTu, file_data, file_vars);

ORF = ORF_data_vars[[1]]
data = ORF_data_vars[[2]]
vars = ORF_data_vars[[3]]

zeroValueRow = which(rowSums(vars)==0)
data = data[-zeroValueRow, ]
vars = vars[-zeroValueRow, ]

probeName = ORF # Rename ORF
probeName = matrix(probeName[-zeroValueRow, ],,1)
dataChip=as.matrix(read.table(file_dataChip))
probe_anno <- read.table(file_annotation, header = FALSE, sep = "\t", quote = "", dec = ".", col.names=c("prob", "anno"), na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = FALSE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)
probeName2=matrix(probe_anno$prob,,1) 
annota=matrix(probe_anno$anno,,1)
TransNames_tab <- read.table(file_transNames, header = FALSE, sep = "\t", quote = "", dec = ".", col.names=c("TN"), na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = FALSE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)
TransNames = matrix(TransNames_tab$TN,,1)
noOfprobe=nrow(probeName)
redundancy=matrix(mat.or.vec(noOfprobe,1),,1)

for (i in 1:noOfprobe){
		redundancy[i,1]=1}

index=matrix(mat.or.vec(nrow(dataChip),1),,1)

for (i in 1: nrow(probeName2)){
	vec <- probeName==probeName2[i,1]
	index[i]=colSums(vec)
	if(index[i]>1){
		pippo=(which(vec))
		redundancy[pippo[2:length(pippo)]]=0
		}
	}

dataChip=dataChip[which(index!=0),]
annota=annota[which(index!=0),]
probeName2=probeName2[which(index!=0),]
probeName2=matrix(probeName2,,1)

probeName=probeName[which(redundancy!=0),]
probeName=matrix(probeName,,1)

data=data[which(redundancy!=0),]
vars=vars[which(redundancy!=0),]

preX=NULL
annotation=NULL
index=mat.or.vec(nrow(data),1)
for (i in 1: nrow(data)){
	index[i]=sum(probeName[i]==probeName2)
	if (index[i]==1)
		
		preX=rbind(preX,dataChip[which(probeName[i]==probeName2),])
		annotation=rbind(annotation,annota[which(probeName[i]==probeName2)])
}

data= data[which(index==1),]
vars= vars[which(index==1),]
probeName= probeName[which(index==1),]
probeName=matrix(probeName,,1)
X= mat.or.vec(nrow(preX),ncol(preX))
I <- c(which(preX<1e-3))
X[I] =1

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