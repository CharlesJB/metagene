\name{getDataFrame}
\alias{getDataFrame}

\title{
  Convert the bootstrapped data into a data.frame
}

\usage{
  getDataFrame(bootstrapData)
}

\arguments{
  \item{bootstapData}{Data produced during the bootstrap analysis}
}
#	A data.frame with the condensed results from the main data structure
#		Columns:
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
