\name{extractReadsInRegion}
\alias{extractReadsInRegion}

\title{
  Extract reads from BAM file that overlap with a specified genomic region.
}

\description{
  This function uses Rsamtools to extract the reads that overlap
  a genomic regions.
}

\usage{
  extractReadsInRegion(
    bamFile,
    chr,
    start,
    end)
}

\arguments{
  \item{bamFile}{Path to the BAM file. The BAM file must exist.}

  \item{chr}{Current chromosome. The chromosome name must 
    be of character type.}

  \item{start}{Starting position of the current region. 
    The starting position has to be a positive integer.}

  \item{end}{Ending position of the current region.
    The ending position has to be a positive integer.}
}
\value{
  \code{extractReadsInRegion} return a \code{data.frame} containing
    every reads overlapping the current genomic regions.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  # Get a feature
  \dontrun{bedFileName <- system.file("extdata/list1.bed", package="metagene")}
  \dontrun{feature <- read.table(bedFileName, header=FALSE, stringsAsFactors=FALSE, nrow=1)}
  \dontrun{colnames(feature) <- c("space", "start_position", "end_position")}

  # Get a bam file
  \dontrun{bamFileName <- system.file("extdata/align1_rep1.bam", package="metagene") }

  # Extract reads
  \dontrun{reads <- metagene:::extractReadsInRegion(bamFileName, feature$space, feature$start_position, feature$end_position)}

  # Convert to density
  \dontrun{metagene:::convertReadsToDensity(reads, feature)}
}
