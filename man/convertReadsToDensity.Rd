\name{convertReadsToDensity}
\alias{convertReadsToDensity}

\title{
  Convert a list of read in a vector of positions.
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
