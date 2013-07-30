# function [geneName, data] = chipTextRead(file)

#% CHIPTEXTREAD reads TXT file for the Spellman data files.
#%
#%	Description:
#%	[geneName, data] = chipTextRead(file)
#%% 	chipTextRead.R version 1.5

chipTextRead = function(data_file) {

#file = "./data/SpellmanMicro.txt"

###dictionary <- read.table("./data/SpellmanMicro.txt", header = TRUE, sep = "\t", quote = "", dec = ".", col.names=c("geneName","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","x18","x19","x20","x21","x22","x23","x24"), na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)
#dictionary$x1[is.na(dictionary$x1)] <- 0

dictionary <- read.table(data_file, header = TRUE, sep = "\t", quote = "", dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

dictionary[is.na(dictionary)] <- 0

geneName = as.matrix(dictionary[,1],,1)
data = dictionary [, 2: ncol(dictionary)] 

#% file ia a string containing the file name and the extension.
#% data is a matrix with 24 columns and N( number of genes) rows   
#[geneName,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14...
#    x15,x16,x17,x18,x19,x20,x21,x22,x23,x24]=...
#    textread(file,'%q %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
#    'headerlines',1,'whitespace','','delimiter','\t');
#data=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24];

geneName_data = list(geneName, data)
return(geneName_data)
}
