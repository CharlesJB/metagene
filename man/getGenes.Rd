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
  Using \code{biomaRt} package, this function will fetch all the
  ensembl_gene_id for a given specie ("human" or "mouse").
}

\value{
  \code{getGenes} return a \code{data.frame} with five columns:
   \enumerate{
    \item feature -> Ensembl gene id
    \item strand -> -1 or 1
    \item space -> chromosome
    \item start_position -> position of the TSS
    \item end_position -> ending position of the last exon of the gene
  }   
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  \dontrun{knownGenes <- getGenes("human")}
}
