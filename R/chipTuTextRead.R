#	CHIPTUTEXTREAD assigns common names to probe IDs 
#	Read the experimental data and Standard errors
#
#	Description:
#	[ORF, data, vars]=chipTuTextRead(file_dictionary, file_probeIDTu, file_data, file_vars)
# 	rChipTuTextRead.R version 0.1.0
#	Written on 09.11.2011

chipTuTextRead <- function(file_dictionary, file_probeIDTu, file_data, file_vars){

# file_dictionary <- "./data/MetabolData/dictionary.txt";
# file_probeIDTu <- "./data/MetabolData/probeIDTu.txt";
# file_data <- "./data/MetabolData/YeastMetabolism_exprs.txt";
# file_vars <- "./data/MetabolData/YeastMetabolism_se.txt";

dictionary <- read.table(file_dictionary, header = FALSE, sep = "\t", quote = "", dec = ".", col.names=c("IDs_temp", "ORF", "ORF_temp","geneCN_temp"), na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = FALSE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

ORF= matrix(dictionary$ORF,,1) #Open Reading Frames
IDs= matrix(dictionary$IDs_temp,,1) # Genes Id

probeIDTu <- read.table(file_probeIDTu, sep=" ", fill= TRUE, col.names=c("IDSn"))

IDSnew <- matrix(probeIDTu$IDSn,,1) # Genes Id

data = as.matrix(read.table(file_data)) # Experimental Data
vars=as.matrix(read.table(file_vars)) # Standard Errors

ORF_data_vars= list(ORF, data, vars);

return(ORF_data_vars)
}

