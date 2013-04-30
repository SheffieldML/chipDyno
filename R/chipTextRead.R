# CHIPTEXTREAD reads TXT file for the Spellman data files.
# CHIPDYNO toolbox
# chipTextRead.R version 1.0.1
# FORMAT chipTextRead <- function(data_file)
# DESC reads TXT file for the Spellman data files.
# ARG file : the Spellman data file
# RETURN geneNames_data : concatenated dataframe of Gene names and 
# point estimate of the expression level
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : chipTuTextRead

chipTextRead <- function(data_file) {

dictionary <- read.table(data_file, header = TRUE, sep = "\t", quote = "", 
				dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, 
				skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, 
				blank.lines.skip = TRUE, comment.char = "#", 
				allowEscapes = FALSE, flush = FALSE)

dictionary[is.na(dictionary)] <- 0

geneName = as.matrix(dictionary[,1],,1)
data = dictionary [, 2: ncol(dictionary)] 

geneName_data = list(geneName, data)
return(geneName_data)
}