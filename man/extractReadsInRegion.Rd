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
  \item{bamFile}{Path to the bam file.}

  \item{chr}{Current chromosome.}

  \item{start}{Starting position of the current region.}

  \item{end}{Ending position of the current region.}
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
  bedFileName <- system.file("extdata/list1.bed", package="MetaFeatures")
  feature <- read.table(bedFileName, header=FALSE, stringsAsFactors=FALSE, nrow=1)
  colnames(feature) <- c("space", "start_position", "end_position")

  # Get a bam file
  bamFileName <- system.file("extdata/align1_rep1.bam", package="MetaFeatures") 

  # Extract reads
  reads <- MetaFeatures:::extractReadsInRegion(bamFileName, feature$space, feature$start_position, feature$end_position)

  # Convert to density
  MetaFeatures:::convertReadsToDensity(reads, feature)
}
