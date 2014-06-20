\name{parseBamFile}
\alias{parseBamFile}

\title{
  Parse a single bam file
}

\description{
  This function will fetch all the aligned reads from a bam file
  that overlap with a list of regions.
}

\usage{
  parseBamFile(
    bamFile,
    features,
    cores=1)
}

\arguments{
  \item{bamFile}{The name of the bam file to parse. Must be sorted and indexed.}

  \item{features}{A data.frame with the infos for every features
    to parse
    Must have the folowing columns: space, start_position, end_position
    start_position and end_position}
  \item{cores}{Number of cores for parallel processing (require parallel package).}
}

\value{
  \code{parseBamFile} returns a \code{list} with an element for
    every feature to parse. Each element contains a vector of
    reads expression.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  # Get bamFile
  \dontrun{bamFileName <- system.file("extdata/align1_rep1.bam", package="MetaMineR") }

  # Get features
  \dontrun{bedFileName <- system.file("extdata/list1.bed", package="MetaMineR")}
  \dontrun{features <- read.table(bedFileName, header=FALSE, stringsAsFactors=FALSE, nrow=10)}
  \dontrun{colnames(features) <- c("space", "start_position", "end_position")}
  \dontrun{features$end_position <- features$start_position + 200}
  \dontrun{features$start_position <- features$start_position - 200}

  # Parse the bam file
  \dontrun{bamDensity <- MetaMineR:::parseBamFile(bamFileName, features)}
}
