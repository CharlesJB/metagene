\name{plotFeatures}
\alias{plotFeatures}

\title{
  Create a metagene plot based on a list of genomic features
}

\description{
  This function produces Metagene-like plots to compare the 
  behavior of DNA-interacting proteins at selected groups of 
  features.
}

\usage{
  plotFeatures(
    bamFiles, 
    features=NULL, 
    specie="hs", 
    maxDistance=5000, 
    design=NULL, 
    binSize=100)
}

\arguments{
  \item{bamFiles}{A \code{vector} of bamFile to plot.}

  \item{features}{Either a filename of a \code{vector} of filenames.
    Supported features: ensembl_gene_id. 
    File must contain a header that correspond to the name of the 
    group.}

  \item{specie}{hs: Homo sapiens (default). mm: Mus musculus.}

  \item{maxDistance}{The distance around feature to include 
    in the plot.}

  \item{design}{A \code{matrix} explaining the relationship 
    between multiple samples. One line per samples. One column
    per group of samples. For example, biological replicates and
    corresponding controls are in the same group.
    1: treatment file(s). 2: control file(s).}

  \item{binSize}{The number of nucleotides in each bin.}
}

\details{
  Pre-calculated features (such as transcription start sites of
  protein coding gene or enhancer) are available. Bam files are
  used to increase the resolution. Multiple combination of group
  of features and or group of bam files can be compared in a
  single analysis. Bootstraping analysis is used to compare the
  groups and locate regions with statistically different 
  enrichment profiles.
}

\value{
  \code{plotFeatures} silently returns a \code{data.frame} with
  the data used to produce the plot.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  # We use the \code{demo.bam} file available with this package:
  bamFile <- system.file("extdata/demo.bam", 
                         package="MetaFeatures")

  # The most basic usage is giving \code{plotFeatures} a single
  # bam file:
  plotFeatures(bamFile)

  # We can also give the function a list of features:
  features <- system.file("extdata/demoFeatures.txt", 
                         package="MetaFeatures")
  plotFeatures(bamFile, features=features)
}
