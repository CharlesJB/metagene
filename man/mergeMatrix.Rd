\name{mergeMatrix}
\alias{mergeMatrix}

\title{
  Convert list of vectors in matrix for a single group
}

\description{
  This function combine a list of list of vectors. The first
  step is to merge every list of vector in a group into a
  matrix. Then every matrices in the same group are combined.
}

\usage{
  mergeMatrix(
    group,
    level)
}

\arguments{
  \item{group}{Correspond to an element in the data subsection
    of the main data structure.}

  \item{level}{The names of the element to merge.}
}

\details{
  Pre-calculated features (such as transcription start sites of
  protein coding gene or enhancer) are available. Bam files are
  used to increase the resolution. Multiple combination of group
  of features and or group of bam files can be compared in a
  single analysis. Bootstraping analysis is used to compare the
  groups and locate regions with statistically different 
  enrichment profiles.
}

\value{
  \code{mergeMatrix} returns the \code{group} element from the
    arguments with an extra element named \code{matrix}.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
