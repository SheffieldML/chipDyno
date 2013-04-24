#% CHIPCHIPTEXTREAD reads TXT file for the Lee ChIP data files.
#% CHIPDYNO toolbox
#% chipChipTextRead.m version 1.5
#% FORMAT function [geneName, annotation, data] = chipChipTextRead(file1, file2)
#% DESC reads TXT file for the Lee ChIP data files.
#% ARG geneName_annotation_file : the geneName and annotation file
#% ARG connectivity_file : the data file
#% RETURN f[[1]]: geneNames
#% RETURN f[[2]]: annotation of the geneNames
#% RETURN f[[3]] : TFA at different experimental point
#% COPYRIGHT : Neil D. Lawrence, 2005
#% COPYRIGHT : Guido Sanguinetti, 2005
#% MODIFICATIONS : Muhammad A. Rahman, 2013
#% SEEALSO : chipTextRead, chipTuTextRead

chipChipTextRead <- function(geneName_annotation_file, connectivity_file){

geneName_annotation  <- read.table(geneName_annotation_file, header = TRUE, sep = "\t", quote = "", dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, skip = 1, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

geneName = as.matrix(geneName_annotation[,1],,1)
annotation = as.matrix(geneName_annotation [,2],,1)

data  <- read.table(connectivity_file, header = FALSE, sep = "\t", quote = "", dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

f= list (geneName, annotation, data)
return (f)
}