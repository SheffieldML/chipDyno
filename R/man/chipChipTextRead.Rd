\name{chipChipTextRead}
\Rdversion{1.0.1}
\alias{chipChipTextRead}
\title{reads TXT file for the Lee ChIP data files.}
\description{
  reads TXT file for the Lee ChIP data files.
}
\usage{
values <- chipChipTextRead(file1, file2)
}
\arguments{
  \item{file1}{the geneName and annotation file}
  \item{file2}{the data file}
}
\value{
  \item{geneName}{geneNames}
  \item{annotation}{annotation of the geneNames}
  \item{data}{TFA at different experimental point}
}
\seealso{
\code{\link{chipTextRead}, \link{chipTuTextRead}}.
}
\examples{
## missing
}
\keyword{model}
