\name{convertReadsToDensity}
\alias{convertReadsToDensity}

\title{
  Convert a list of read in a vector of positions.
}

\description{
  This function takes a list of reads obtained with Rsamtools
  and convert them into a \code{vector} of density.
}

\usage{
  convertReadsToDensity(
    currentReads,
    currentFeature)
}

\arguments{
  \item{currentReads}{The list of read to parse.}

  \item{currentFeature}{The feature to parse.}
}

\value{
  \code{convertReadsToDensity} return a \code{vector} with the
    coverage of every positions calculated from the reads around
    the max distance from TSS.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
