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
    specie="human", 
    maxDistance=5000, 
    design=NULL, 
    binSize=100,
    alpha=0.05,
    sampleSize=1000,
    cores=1)
}

\arguments{
  \item{bamFiles}{A \code{vector} of bamFile to plot.}

  \item{features}{Either a filename of a \code{vector} of filenames.
    Supported features: ensembl_gene_id. 
    File must contain a header that correspond to the name of the 
    group.}

  \item{specie}{human: Homo sapiens (default). mouse: Mus musculus.}

  \item{maxDistance}{The distance around feature to include 
    in the plot.}

  \item{design}{A \code{matrix} explaining the relationship 
    between multiple samples. One line per samples. One column
    per group of samples. For example, biological replicates and
    corresponding controls are in the same group.
    1: treatment file(s). 2: control file(s).}

  \item{binSize}{The number of nucleotides in each bin.}

  \item{alpha}{Confidence interval to include in the plot. For
    a value of 0.05, the confidence interval between 0.05 and
    0.95 will be shown in the plot.}

  \item{sampleSize}{Number of time each bin will be resampled.
    Should be at least 1000.}

  \item{cores}{Number of cores for parallel processing.
    Require parallel package.}
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
  \code{plotFeatures} silently returns an instance of the main
  data structure used to produce the plot.

  The main data structure is a list of list.
  The first level contains all the groups to be plotted.
  The second level contains an element for each bam file for the
    current group.

  When the data structure is completed, each group contains the
    following elements:
      * featureName
      * designName
      * bamFiles
      * matrix
      * bootstrap
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
