\name{parseBamFile}
\alias{parseBamFile}

\title{
  Parse a single bam file
}

\usage{
  parseBamFile(
    bamFile,
    features,
    cores=1)
}

\arguments{
  \item{bamFile}{The name of the bam file to parse. Must be sorted and indexed.}

  \item{features}{A data.frame with the infos for every features
    to parse
    Must have the folowing columns: feature, strand, space,
    start_position and end_position}
  \item{cores}{Number of cores for parallel processing (require parallel package).}
}

\value{
  \code{parseBamFile} returns a \code{list} with an element for
    every feature to parse. Each element contains a vector of
    reads expression.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
