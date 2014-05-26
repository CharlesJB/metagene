\name{applyOnGroups}

\alias{applyOnGroups}

\title{
  Apply a function on every groups of the main data structure
}

\usage{
  applyOnGroups(
    groups, 
    cores=1,
    FUN,
    ...)
}

\arguments{
  \item{groups}{The main data structure}

  \item{features}{Either a filename of a \code{vector} of filenames.
    Supported features: ensembl_gene_id. 
    File must contain a header that correspond to the name of the 
    group.}

  \item{cores}{Number of cores for parallel processing.
    Require parallel package.}

  \item{FUN}{The function to apply on every groups.}

  \item{...}{Extract arguments for FUN.}
}

\details{
  The function must be aware of the structure of the main data
  structure. See the \code{value} section of \code{plotFeatures}
  and \code{plotFeaturesByRegions} manuals.
}

\value{
  \code{plotFeatures} silently returns a the main data structure
    with the data used to produce the plot.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  # We use the \code{demo.dataStructure.Rda} file available with
  # this package:
  bamFile <- system.file("extdata/demo.dataStructure.Rda", 
                         package="MetaFeatures")

  # We create a simple function that will return the name of the
  # features for each groups
  extractFeatureNames <- function(group) {
    return(g$featureName)
  }
  featureNames <- applyOnGroups(groups, FUN=extractFeatureNames)
}
