\name{prepareFeatures}
\alias{prepareFeatures}

\title{
  Convert a list of IDs into a GRangesList
}

\description{
  This function will convert a list of Ensembl Gene IDs into
  a list of regions.
}

\usage{
  prepareFeatures(
    features,
    specie="human",
    maxDistance=5000,
    cores=1)
}

\arguments{
  \item{features}{Either a filename of a vector of filenames.
    Supported features: ensembl_gene_id
    If value is NULL: every ensembl_gene_id will be returned.}

  \item{specie}{human: Homo sapiens (default) / mouse: Mus musculus}

  \item{maxDistance}{The distance around features to include in
    the plot.}

  \item{cores}{Number of cores for parallel processing (require
    parallel package).}
}

\value{
  \code{prepareFeatures} returns a \code{list} of \code{GRangesList}.
    One \code{GRanges} by group of features. The names of
    each GRanges of the list correspond to the name of the
    group.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  # Load a feature file
  \dontrun{featuresFileName <- system.file("extdata/list1.txt", package="metagene")}
  \dontrun{features <- metagene:::prepareFeatures(featuresFileName, specie="mouse")}
}
