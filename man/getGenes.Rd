\name{getGenes}
\alias{getGenes}

\title{
  Fetch the annotation of all genes.
}

\description{
  This function will fetch the positions of all known coding
  genes for a given specie.
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
  Using \code{biomaRt} package, this function will fetch all the
  ensembl_gene_id for a given specie ("human" or "mouse").
}

\value{
  \code{getGenes} return a \code{data.frame} with five columns:
    1: feature -> Ensembl gene id
    2: strand -> -1 or 1
    3: space -> chromosome
    4: start_position -> position of the TSS
    5: end_position -> ending position of the last exon of the gene
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{
  knownGenes <- getGenes("human")
  }
}
