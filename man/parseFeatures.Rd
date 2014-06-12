\name{parseFeatures}
\alias{parseFeatures}
\title{
  Parse an experiment using a list of features
}

\description{
  This function produces the list object that contains all the
  information necessary to produce a metagene-like plot with the
  \code{plot.graphic} function.
}

\usage{
  parseFeatures(
    bamFiles, 
    features=NULL, 
    specie="human", 
    maxDistance=5000, 
    design=NULL, 
    cores=1,
    debug=FALSE)
}

\arguments{
  \item{bamFiles}{A \code{vector} of bamFile to plot.}

  \item{features}{Either a filename of a \code{vector} of filenames.
    Supported features: ensembl_gene_id. 
    File must contain a header that correspond to the name of the 
    group.
    If value is NULL, all known RefSeq genes will be used.}

  \item{specie}{human: Homo sapiens (default). mouse: Mus musculus.}

  \item{maxDistance}{The distance around feature to include 
    in the plot.}

  \item{design}{A \code{matrix} explaining the relationship 
    between multiple samples. One line per samples. One column
    per group of samples. For example, biological replicates and
    corresponding controls are in the same group.
    1: treatment file(s). 2: control file(s).}

  \item{cores}{Number of cores for parallel processing.
    Require parallel package.}

  \item{debug}{Keep the intermediate files (can use a lot of memory).
    TRUE or FALSE.}
}

\details{
  This function will extract the read density from alignments
  files (bam) in the viscinity of transcription start sites of
  one or multiple list of genes.

  The values are normalized as read per millions aligned (RPM).

  It is possible to parse multiple groups of gene by saving each
  list in a separate file and by listing the file names in a
  vector as the \code{features} parameter.

  By using the \code{design} parameter, the \code{parseFeatures}
  function will deal with more complex experimental design such
  as the use of replicates and/or controls. The values of
  controls are substracted from every replicates.
}

\value{
  \code{plotFeatures} returns a list that contains the data
  necessary to produce a plot.

  The data structure is a list of lists.

  The first level contain the following fields:
    * design: The information from the design file.
    * param: The values of the argument used with parseFeatures.
    * bamFilesDescription: A \code{data.frame} with the
      following columns;
        * The names of the original bam files, the names of the
        * sorted bam files and the number of aligned reads.
    * matrix: A list of matrix that will be used to produce the
        plot. One element by combination of features/design
        groups.
}

\author{
  Charles Joly Beauparlant
  <Charles.Joly-Beauparlant@crchul.ulaval.ca>
}

\examples{
  # We use the \code{demo.bam} file available with this package:
  \dontrun{bamFile <- system.file("extdata/demo.bam", 
                         package="MetaFeatures")}

  # The most basic usage is giving \code{plotFeatures} a single
  # bam file:
  \dontrun{groups <- parseFeatures(bamFile)}

  # We can also give the function a list of features:
  \dontrun{features <- system.file("extdata/demoFeatures.txt", 
                         package="MetaFeatures")}
  \dontrun{parseFeatures(bamFile, features=features)}
}
