\name{prepareBamFiles}
\alias{prepareBamFiles}

\title{
  Prepare bam files before parsing.
}

\description{
  Sort and index bam files, if necessary. Return the number of
  aligned reads for each bam file.
}

\usage{
  prepareBamFiles(
    bamFiles,
    cores=1)
}

\arguments{
  \item{bamFiles}{Vector containing the list of every bam
    filename to be included in the analysis.}

  \item{cores}{Number of cores used by the function.}
}

\details{
  This function will take a list of bam files and first check if
  they are indexed. If this is not the case, it will call the
  \code{Rsamtools} package to sort and index the files.

  The names of the orginal bam files will be stored in the
  \code{oldBam} column of the results and the sorted bam file
  name will be stored in the \code{bam} column. If the files
  already indexed, the names in both columns will be the same.

  The aligned read counts are then calculated for each files and
  added to the \code{alignedCount} column of the result.
}

\value{
  \code{prepareBamFiles} returns a \code{data.frame} with three
    columns and as many rows as there are bam files in the
    \code{bamFiles} argument. The \code{data.frame} has the
    following columns:
      * \code{oldBam}
      * \code{bam}
      * \code{alignedCount}
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  bamFileName <- system.file("extdata/align1_rep1.bam", package="MetaFeatures") 
  bamFileDescription <- MetaFeatures:::prepareBamFiles(bamFileName)
}
