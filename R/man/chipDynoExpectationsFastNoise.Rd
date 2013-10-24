\name{chipDynoExpectationsFastNoise}
\Rdversion{1.0.1}
\alias{chipDynoExpectationsFastNoise}
\title{computes posterior expectations TFA}
\description{
  computes posterior expectations of TFA considering uncertainty of the expression level.
}
\usage{
values <- chipDynoExpectationsFastNoise(data, X, Sigma, beta,
                                        precs, gamma, mu,
                                        TransNames, annotations,
                                        transName, geneName)
}
\arguments{
  \item{data}{point estimate of the expression level}
  \item{X}{connectivity measurement between genes and transcription factors}
  \item{Sigma}{prior covariance matrix of TFA}
  \item{beta}{}
  \item{precs}{uncertainty of the expression level}
  \item{gamma}{degree of temporal continuity}
  \item{mu}{mean value of the transcription factor activity}
  \item{TransNames}{Transcription factors}
  \item{annotations}{Gene names}
  \item{transName}{specific transcription factor}
  \item{geneName}{specific gene name}
}
\value{
  \item{tf}{gene specific transcription factor activity}
  \item{tfErrors}{error in gene specific transcription factor activity}
  \item{tfErrorsDiffs}{}
}
\seealso{
\code{\link{chipDynoExpectationsFast}}.
}
\examples{
## missing
}
\keyword{model}
