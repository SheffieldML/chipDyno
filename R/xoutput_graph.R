output_file = "xoutput_ana.txt"

output_col  <- read.table(output_file, header = FALSE, sep = " ", quote = "", dec = ".", , na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE)

cycle = as.matrix(output_col[[3]])
error = as.matrix(output_col[[6]])
scale = as.matrix(output_col[[9]])

M= array (1:2, dim <- c(2,1))
layout(M)

plot(cycle[150:200], scale[150:200], type='l', col='red')
plot(cycle[150:200], error[150:200], type='l', col='green4')


