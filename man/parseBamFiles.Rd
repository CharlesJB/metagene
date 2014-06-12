\name{parseBamFiles}
\alias{parseBamFiles}

\title{
  Parse multiple bamFiles
}

\usage{
  parseBamFiles(
    bamFiles,
    featuresGroups,
    cores=1)
}

\arguments{
  \item{bamFiles}{A vector of bam files}

  \item{featuresGroups}{A list of data.frame. One data.frame by
    group of features.
      The names of each element of the list correspond to the
      name of the group.}

  \item{cores}{Number of cores for parallel processing (require parallel package).}
}

\value{
  \code{parseBamFiles} returns a \code{list} of \code{list} of \code{list}
    that contains the raw counts for every features groups.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
