\name{getGenesBiomart}
\alias{getGenesBiomart}

\title{
  Fetch the annotation of all genes from biomart.
}

\description{
  This function will fetch the positions of all known coding
  genes for a given specie. Currently supported species are: 
  \dQuote{mouse}, \dQuote{human} (default).

  This function was used to create external datasets for the
  \code{getGenes} function.
}

\usage{
  getGenesBiomart(
    specie="human")
}

\arguments{
  \item{specie}{human: Homo sapiens (default) / mouse: Mus
    musculus}
}

\details{
  Using \code{biomaRt} package, this function will fetch all the
  ensembl_gene_id for a given specie ("human" or "mouse").
}

\value{
  \code{getGenesBiomart} return a \code{GRanges} object with a
  feature metadata that corresponds to emsembl_gene_id.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{knownGenes <- getGenesBiomart("human")}
}
