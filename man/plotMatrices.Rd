\name{plotMatrices}
\alias{plotMatrices}

\title{
  Create a graph
}

\description{
  This function wraps all the steps necessary to produce the
  final plot using the results of either \code{parseFeatures}
  or \code{parseRegions} functions.
}

\usage{
  plotMatrices(
    matricesGroups,
    data,
    binSize=100,
    alpha=0.05,
    sampleSize=1000,
    cores=1)
}

\arguments{
  \item{matricesGroups}{A list with every groupsto include in
    the graph. A group is one or more combination of
    featureGroup/designGroup. The names must correspond to the
    names of the matrix object returned from the parse function.}

  \item{data}{The object returned by the parse functions.}

  \item{binSize}{The number of nucleotides in each bin for the
    bootstrap step.}

  \item{alpha}{Confidence interval.}

  \item{sampleSize}{Number of time each bin will be resampled (
    hould be at least 1000).}

  \item{cores}{Number of cores for parallel processing (require
    parallel package).}
}

\value{
  \code{plotMatrices} return the \code{data.frame} used to
    produce the graph.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  # Prepare the groups
  bamFileName <- system.file("extdata/align1_rep1.bam", package="MetaFeatures")
  listFileName <- system.file("extdata/list1.txt", package="MetaFeatures")
  \dontrun{
  groups <- parseFeatures(bamFileName, listFileName, specie="mouse")

  # Do the graph
  DF <- plotMatrices(list(g1="list1"), groups)
  }
}
