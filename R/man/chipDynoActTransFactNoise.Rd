\name{chipDynoActTransFactNoise}
\Rdversion{1.0.1}
\alias{chipDynoActTransFactNoise}
\title{identifies significantly varying TFs with uncertainty of..}
\description{
  identifies significantly varying TFs.
}
\usage{
values <- chipDynoActTransFactNoise(data, X, Sigma, beta, precs,
                                    gamma, mu, TransNames,
                                    annotation, sigLev)
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
  \item{annotation}{Gene names}
  \item{sigLev}{threshold value}
}
\value{
  \item{list}{list of regulators for a specific gene}
  \item{newX}{}
  \item{newXVals}{}
}
\seealso{
\code{\link{chipDynoTransFact}, \link{chipDynoTransFactNoise}, \link{chipDynoActTransFact}}.
}
\examples{
## missing
}
\keyword{model}
