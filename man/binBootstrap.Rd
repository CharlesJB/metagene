\name{binBootstrap}
\alias{binBootstrap}

\title{
  Estimate mean and confidence interval of a column using bootstrap.
}

\usage{
  binBootstrap(
    data,
    alpha,
    sampleSize,
    cores=1)
}

\arguments{
  \item{data}{A vector representing a column from the binned
    matrix.}

  \item{alpha}{Confidence interval.}

  \item{sampleSize}{Number of time each bin will be resampled (
    hould be at least 1000).}

  \item{cores}{Number of cores for parallel processing (require
    parallel package).}
}

\value{
  \code{binBootstrap} returns a \code{list} with the mean of the bootstrapped
    \code{vector} and the quartile of order alpha/2 and (1-alpha/2)
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
