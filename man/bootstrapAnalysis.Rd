\name{bootstrapAnalysis}
\alias{bootstrapAnalysis}

\title{
  Perform the bootstrap analysis
}

\description{
  This function wraps the diffentent steps needed to do a
  bootstrap analysis on a matrix.
}

\usage{
  bootstrapAnalysis(
    currentMatrix,
    binSize,
    alpha,
    sampleSize,
    cores=1)
}

\arguments{
  \item{currentMatrix}{The matrix to use for the bootstrap
    analysis}

  \item{binSize}{The number of nucleotides in each bin for the
    bootstrap step.}

  \item{alpha}{Confidence interval.}

  \item{sampleSize}{Number of time each bin will be resampled (
    hould be at least 1000).}

  \item{cores}{Number of cores for parallel processing (require
    parallel package).}
}

\value{
  \code{bootstrapAnalysis} returns a \code{list} with the mean of the
    bootstraped vector and the quartile of order alpha/2 and
    (1-alpha/2).
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{MetaMineR:::bootstrapAnalysis(matrix(1:40, nrow=4, ncol=10), 4, alpha=0.05, sampleSize=100)}
}
