\name{chipDynoActTransFact}
\Rdversion{1.0.1}
\alias{chipDynoActTransFact}
\title{identifies significantly varying TFs.}
\description{
  identifies significantly varying TFs.
}
\usage{
values <- chipDynoActTransFact(data, X, Sigma, beta, gamma, mu,
                               TransNames, annotation, sigLev)
}
\arguments{
  \item{data}{point estimate of the expression level}
  \item{X}{connectivity measurement between genes and transcription factors}
  \item{Sigma}{prior covariance matrix}
  \item{beta}{}
  \item{gamma}{}
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
\code{\link{chipDynoActTransFact}}.
}
\examples{
## missing
}
\keyword{model}
