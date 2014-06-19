\name{plotGraphic}
\alias{plotGraphic}

\title{
  Produce a plot with based on a data.frame
}

\description{
  This function will produce the final plot for an analysis.
}

\usage{
  plotGraphic(
    DF,
    title)
}

\arguments{
  \item{DF}{The data frame produced by the getDataFrame
    function}

  \item{title}{The title of the graph}
}

\value{
  \code{plotGraphic} plots a graph on the current device.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  DF <- data.frame(Groups=rep("group1",11), distances=seq(-5,5), means=rep(30, 11), qinf=rep(20, 11), qsup=rep(40, 11))
  MetaMineR:::plotGraphic(DF, "Example")
}
