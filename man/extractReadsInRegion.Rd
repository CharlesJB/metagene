\name{extractReadsInRegion}
\alias{extractReadsInRegion}

\title{
  Extract reads from BAM file that overlap with a specified genomic region.
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
