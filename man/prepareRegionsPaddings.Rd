\name{prepareRegionsPaddings}
\alias{prepareRegionsPaddings}

\title{
  Get the genomic position of the paddings
}

\description{
  This function will return the regions corresponding to the
  padding size on the side specified with the function call
}

\usage{
  prepareRegionsPaddings(
    regionsGroups,
    side,
    paddingSize=2000,
    cores=1)
}

\arguments{
  \item{regionsGroups}{The output of prepareRegions or
    prepareFeatures functions.}

  \item{side}{The side of the padding to prepare. Either 'left'
    or 'right'.}

  \item{paddingSize}{The length padding we want to add on each
    side of each regions.}

  \item{cores}{Number of cores for parallel processing (require
    parallel package).}
}

\value{
  \code{prepareRegionsPaddings} return a list of \code{data.frame}
    ,each \code{data.frame} corresponding to the paddings for a
    group of features.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{regionsFileName <- system.file("extdata/list1.bed", package="metagene")}
  \dontrun{regions <- metagene:::prepareRegions(regionsFileName)}
  \dontrun{regionsPaddingsLeft <- metagene:::prepareRegionsPaddings(regions, "left", paddingSize=200)}
  \dontrun{regionsPaddingsRight <- metagene:::prepareRegionsPaddings(regions, "right", paddingSize=200)}
}
