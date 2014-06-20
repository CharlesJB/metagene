\name{scaleVectors}
\alias{scaleVectors}

\title{
  Resize the vectors of every group so they have the same length.
}

\description{
  This function takes a list of list of vectors and scale them
  so they have the same size}

\usage{
  scaleVectors(
    group,
    level,
    domain,
    scaleCores=1)
}

\arguments{
  \item{group}{A \code{list} that contains a list of list of
    vectors}

  \item{level}{The names of the element of the group to scale}

  \item{domain}{The target length for the vectors}

  \item{scaleCores}{Number of cores for parallel processing 
    (require parallel package)}
}

\details{
  A group can contain multiple list of list of vectors (i.e.: if
  \code{debug == TRUE}). We need to to specify which one to use
  with the \code{level} argument.

  Each second level list can contain multiple list corresponding
  to each bam files associated with current group. Each vector
  of each bam files sublist will be scaled to the size specified
  with the \code{domain} argument.
}

\value{
  \code{scaleVectors} returns the same group that was used in
    input with an extra element named scaled.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{group <- list(l1=list(element1=list(1:10, 1:100)))}
  \dontrun{MetaMineR:::scaleVectors(group, "l1", 100)}
}
