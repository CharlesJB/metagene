\name{applyOnGroups}

\alias{applyOnGroups}

\title{
  Apply a function on every groups of the main data structure
}

\description{
  This function takes a \code{list} of elements and apply a
  user-specified function on them. Works like lapply, but can
  use multiple cores if asked.
}

\usage{
  applyOnGroups(
    groups,
    cores=1,
    FUN,
   ...)
}

\arguments{
  \item{groups}{The main data structure}

  \item{cores}{Number of cores for parallel processing.
    Require parallel package.}

  \item{FUN}{The function to apply on every groups.}

  \item{...}{Extract arguments for FUN.}
}
\value{
  The \code{applyOnGroups} return value will vary depending on
  the function specified by the user.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
