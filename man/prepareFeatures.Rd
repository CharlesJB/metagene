\name{prepareFeatures}
\alias{prepareFeatures}

\title{
  Convert a list of IDs into a list of genomic regions
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
    File must contain a header that correspond to the name of
    the group.
    If value is NULL: every ensembl_gene_id will be returned.}

  \item{specie}{human: Homo sapiens (default) / mouse: Mus musculus}

  \item{maxDistance}{The distance around feature to include in
    the plot.}

  \item{cores}{Number of cores for parallel processing (require
    parallel package).}
}

# Output:
#	A list of data.frame. One data.frame by group of features.
#	The names of each element of the list correspond to the name of the group.
\value{
  \code{prepareFeatures} returns a \code{list} of \code{data.frames}.
    One \code{data.frame} by group of features. The names of
    each element of the list correspond to the name of the
    group.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
