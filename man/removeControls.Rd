\name{removeControls}
\alias{removeControls}

\title{
  Substract controls from a single group
}

\description{
  This function will substract controls for every replicates in
  the same group. If there are multiple controls, their values
  will be averaged before substracting them from the input.
  Values cannot be lower than 0.
}

\usage{
  removeControls(
    group,
    data.rpm,
    design,
    controlCores=1)
}

\arguments{
  \item{group}{A group extracted from the main data structure}

  \item{data.rpm}{The normalized data for every bam files.}

  \item{design}{The line from matrix explaining the relationship
    between current samples.}

  \item{controlCores}{Number of cores for parallel processing (
    require parallel package).}
}

\value{
  \code{removeControls} The group extracted from the main data
    structure from which the controls were substracted then
    deleted
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  # Prepare bam files
  bamFileInput <- system.file("extdata/align1_rep1.bam", package="MetaMineR")
  bamFileControl <- system.file("extdata/ctrl.bam", package="MetaMineR")
  bamFiles <- c(bamFileInput, bamFileControl)
  groups <- list()
  groups$bamFilesDescription <- MetaMineR:::prepareBamFiles(bamFiles)

  # Parse features and convert to RPM
  features <- system.file("extdata/list1.txt", package="MetaMineR")
  groups$design <- data.frame(Samples=c(bamFileInput, bamFileControl), Exp1=c(1,2))
  \dontrun{
  groups$regionsGroups <- MetaMineR:::prepareFeatures(features=features, specie="mouse", maxDistance=1000)
  groups$raw <- MetaMineR:::parseBamFiles(groups$bamFilesDescription$bam, groups$regionsGroups)
  groups$rpm <- MetaMineR:::rawCountsToRPM(groups$raw, groups$bamFilesDescription)
  groups$data <- MetaMineR:::prepareGroups(names(groups$regionsGroups), bamFiles=groups$bamFilesDescription, design=groups$design)

  # Remove controls
  groups$data <- MetaMineR:::applyOnGroups(groups=groups$data, cores=1, FUN=MetaMineR:::removeControls, data.rpm=groups$rpm, design=groups$design)
  }
}
