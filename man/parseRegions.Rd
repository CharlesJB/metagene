\name{parseRegions}
\alias{parseRegions}
\title{
  Parse an experiment using regions that can be of different length.
}

\description{
  This function produces the list object that contains all the
  information necessary to produce a metagene-like plot with the
  \code{plot.graphic} function.
}

\usage{
  parseRegions(
    regions, 
    bamFiles, 
    specie="human", 
    design=NULL, 
    paddingSize,
    cores=1,
    debug=FALSE)
}

\arguments{
  \item{regions}{A \code{vector} of bed file names corresponding
    to the regions to include in the analysis.}

  \item{bamFiles}{A \code{vector} of bamFile to plot.}

  \item{specie}{human: Homo sapiens (default). mouse: Mus musculus.}

  \item{design}{A \code{matrix} explaining the relationship 
    between multiple samples. One line per samples. One column
    per group of samples. For example, biological replicates and
    corresponding controls are in the same group.
    1: treatment file(s). 2: control file(s).}

  \item{paddingSize}{The length of padding we want to add on
    each side of each regions.}

  \item{cores}{Number of cores for parallel processing.
    Require parallel package.}

  \item{debug}{Keep the intermediate files (can use a lot of memory).
    TRUE or FALSE.}
}

\details{
  This function will extract the read density from alignments
  files (bam) in one or multiple list of regions.

  The values are normalized as read per millions aligned (RPM).

  It is possible to parse multiple groups of regions by saving
  each regions in a separate bed file and by listing the file
  names in a vector as the \code{regions} parameter.

  By using the \code{design} parameter, the \code{parseRegions}
  function will deal with more complex experimental design such
  as the use of replicates and/or controls. The values of
  controls are substracted from every replicates.
}

\value{
  \code{plotRegions} returns a list that contains the data
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
  # We can also give the function a list of features:
  \dontrun{regions <- system.file("extdata/demoRegions.txt", 
                         package="MetaFeatures")}
  \dontrun{parseRegions(bamFile, features=features)}
}
