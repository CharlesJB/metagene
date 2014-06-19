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

\examples{
  # Create the lists of vector
  vector1 <- list(c(1,2,3,4), c(5,6,7,8))
  vector2 <- list(c(9,10,11,12), c(13,14,15,16))

  # Add them to a list
  group <- list()
  group$vectors <- list(vector1, vector2)

  # We will modify the group objet to add an element names
  # matrix that contains the matrix that correspond to the
  # vectors binded by row
  MetaMineR:::mergeMatrix(group, "vectors")
}
