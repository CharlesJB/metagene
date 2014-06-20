\name{binMatrix}
\alias{binMatrix}

\title{
  Bin matrix columns
}


\description{
  This function will split a matrix into multiple segments.
}

\usage{
  binMatrix(
    data,
    binSize)
}

\arguments{
  \item{data}{The matrix to bin.}

  \item{binSize}{The number of nucleotides in each bin.}
}

\value{
  \code{binMatrix} returns a \code{matrix} with each columns
    representing the mean of \code{binSize} nucleotides.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{MetaMineR:::binMatrix(matrix(1:40, nrow=4, ncol=10), 4)}
}
