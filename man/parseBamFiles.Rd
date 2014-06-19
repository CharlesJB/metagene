\name{parseBamFiles}
\alias{parseBamFiles}

\title{
  Parse multiple bamFiles
}

\description{
  This function will fetch the aligned reads from multiple bam
  files that overlap with a list of regions.
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

\examples{
  # Get a bam file name
  bamFileName <- system.file("extdata/align1_rep1.bam", package="MetaFeatures")

  # Load a feature file
  featuresFileName <- system.file("extdata/list1.txt", package="MetaFeatures")
  features <- MetaFeatures:::prepareFeatures(featuresFileName, specie="mouse")

  # Parse multiple bam files
  bam <- MetaFeatures:::parseBamFiles(bamFileName, features)
}
