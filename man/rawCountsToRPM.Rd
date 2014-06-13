\name{rawCountsToRPM}
\alias{rawCountsToRPM}

\title{
  Convert the raw counts into reads per million aligned (rpm)
}

\description{
  This function will convert raw counts into reads per million
  aligned (rpm) using the aligned reads counts obtained with the
  \code{prepareBamFiles} function.
}

\usage{
  rawCountsToRPM(
    rawCounts,
    bamFilesDescription,
    cores=1)
}

\arguments{
  \item{rawCounts}{The data structure returned by the
    parseBamFiles function.}

  \item{bamFilesDescription}{The data.frame obtained with the
    prepareBamFiles function.}

  \item{cores}{Number of cores for parallel processing (require
    parallel package).}
}

\value{
  \code{rawCountsToRPM} returns a \code{list} of \code{list} of \code{list}
    that contains the rpm for every features groups.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
