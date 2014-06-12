\name{prepareGroups}
\alias{prepareGroups}

\title{
  Distribute the bam filenames in their respective groups.
}

\usage{
  prepareGroups(
    featuresGroupsNames,
    bamFiles,
    design=NULL)
}

\arguments{
  \item{featuresGroupsNames}{The name of the group of features
    in the current analysis.}

  \item{bamFiles}{A vector of bam filename(s).}

  \item{design}{A matrix explaining the relationship between
    multiple samples.
      * One line per samples.
      * One column per group of samples. For example, biological
        replicates and corresponding controls are in the same group.
          1: treatment file(s)
          2: control file(s)}
}
\value{
  \code{prepareGroups} returns a \code{list} of \code{list} of \code{list}.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
