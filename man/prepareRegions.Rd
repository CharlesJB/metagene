\name{prepareRegions}
\alias{prepareRegions}

\title{
  Parse bed files and convert them in a list of GRangesList.
}

\description{
  This function will convert a bed files into a \code{GRangesList}.
}


\usage{
  prepareRegions(
    regions,
    cores=1)
}

\arguments{
  \item{regions}{A vector of bed file names corresponding to the
    regions to include in the analysis. The file name (minus the
    extension) will be used as the name of the region.}

  \item{cores}{Number of cores for parallel processing (require parallel package).}
}

\value{
  \code{prepareRegions} returns a \code{GRangesList}.
    One \code{GRanges} by group of features. The names of each
    element of the list correspond to the name of the group.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{regionsFileName <- system.file("extdata/list1.bed", package="metagene")}
  \dontrun{regions <- metagene:::prepareRegions(regionsFileName)}
}
