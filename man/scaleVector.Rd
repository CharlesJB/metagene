\name{scaleVector}
\alias{scaleVector}

\title{
  Scale the values of a vector to fit with predetermined size
}

\description{
  This function will resize the values of the vector so that
  they fit in a pre-defined range. The number of elements will
  changed, but the general distribution will remain the same.
}

\usage{
  scaleVector(
    values,
    domain)
}

\arguments{
  \item{values}{the values to scale}

  \item{domain}{the range to fit the value to}
}

\value{
  \code{scaleVector} returns a \code{vector} with the scaled data
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
