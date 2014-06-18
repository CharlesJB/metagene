\name{getDataFrame}
\alias{getDataFrame}

\title{
  Convert the bootstrapped data into a data.frame
}

\description{
  This function will convert the results of the bootstrap
  analysis into a \code{data.frame} that can be used directly
  with ggplot2 to produce the final plot.
}

\usage{
  getDataFrame(bootstrapData, range)
}

\arguments{
  \item{bootstrapData}{Data produced during the bootstrap analysis}

  \item{range}{The start and the end for the x-label}
}

\value{
  \code{getDataFrame} a \code{data.frame} with the condensed results
    from the main data structure.

    Columns:
      * Groups: name of current group
      * distances: the number of bin for each entry
      * means: the means to plot
      * qinf: the lower end of the confidence interval
      * qsup: the higher end of the confidence interval
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}
