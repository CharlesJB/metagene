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
  \dontrun{bamFileInput <- system.file("extdata/align1_rep1.bam", package="metagene")}
  \dontrun{bamFileControl <- system.file("extdata/ctrl.bam", package="metagene")}
  \dontrun{bamFiles <- c(bamFileInput, bamFileControl)}
  \dontrun{groups <- list()}
  \dontrun{groups$bamFilesDescription <- metagene:::prepareBamFiles(bamFiles)}

  # Parse features and convert to RPM
  \dontrun{features <- system.file("extdata/list1.txt", package="metagene")}
  \dontrun{groups$design <- data.frame(Samples=c(bamFileInput, bamFileControl), Exp1=c(1,2))}
  \dontrun{groups$regionsGroups <- metagene:::prepareFeatures(features=features, specie="mouse", maxDistance=1000)}
  \dontrun{groups$raw <- metagene:::parseBamFiles(groups$bamFilesDescription$bam, groups$regionsGroups)}
  \dontrun{groups$rpm <- metagene:::rawCountsToRPM(groups$raw, groups$bamFilesDescription)}
  \dontrun{groups$data <- metagene:::prepareGroups(names(groups$regionsGroups), bamFiles=groups$bamFilesDescription, design=groups$design)}

  # Remove controls
  \dontrun{groups$data <- metagene:::applyOnGroups(groups=groups$data, cores=1, FUN=metagene:::removeControls, data.rpm=groups$rpm, design=groups$design)}
}
