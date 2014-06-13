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
