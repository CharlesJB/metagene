\name{getGenes}
\alias{getGenes}

\title{
  Fetch the annotation of all genes.
}

\description{
  This function will fetch the positions of all known coding
  genes for a given specie. Currently supported species are: 
  \dQuote{mouse}, \dQuote{human} (default).
}

\usage{
  getGenes(
    specie="human")
}

\arguments{
  \item{specie}{human: Homo sapiens (default) / mouse: Mus
    musculus}
}

\details{
  This function will fetch all the ensembl_gene_id for a given
  specie ("human" or "mouse").
}

\value{
  \code{getGenes} return a \code{GRanges} object with a
  feature metadata that correspond to the ensembl_gene_id.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{knownGenes <- getGenes("human")}
}
