# CHIPTUTEXTREAD assigns common names to probe IDs
# CHIPDYNO toolbox
# chipTuTextRead.R version 1.0.1
# FORMAT [ORF, data, vars]=chipTuTextRead()
# DESC assigns common names to probe IDs
# RETURN ORF : concatenated dataframe of (1) Gene names,
# RETURN data : (2) point estimate of the expression level, 
# (3) uncertainty of the expression level
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipTextRead

chipTuTextRead <- function(file_dictionary, file_probeIDTu, file_data, file_vars){

dictionary <- read.table(file_dictionary, header = FALSE, sep = "\t", 
				quote = "", dec = ".", col.names=c("IDs_temp", "ORF", 
				"ORF_temp","geneCN_temp"), na.strings = "NA", 
				colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, 
				fill = TRUE, strip.white = FALSE, blank.lines.skip = FALSE, 
				comment.char = "#", allowEscapes = FALSE, flush = FALSE)

ORF= matrix(dictionary$ORF,,1) 
IDs= matrix(dictionary$IDs_temp,,1)

probeIDTu <- read.table(file_probeIDTu, sep=" ", fill= TRUE, col.names=c("IDSn"))
IDSnew <- matrix(probeIDTu$IDSn,,1)
data = as.matrix(read.table(file_data))
vars=as.matrix(read.table(file_vars))
ORF_data_vars= list(ORF, data, vars)

return(ORF_data_vars)
}