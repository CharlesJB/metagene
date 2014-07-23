# Created by Charles Joly Beauparlant
# 2013-11-26

# Parse an experiment using a list of features
#
# Input:
#    bamFiles:      A vector of bamFile to plot. TODO: Should also accept a list
#                   of bamfiles where each elements would be grouped together
#    features:      Either a filename of a vector of filenames.
#                   Supported features:
#                       * ensembl_gene_id
#                   File must contain a header that correspond to the name of
#                   the group.
#    specie:        human: Homo sapiens (default) / mouse: Mus musculus
#    maxDistance:   The distance around feature to include in the plot.
#    design:        A data.frame explaining the relationship between multiple
#                   samples.
#                       * One line per samples.
#                       * One column per group of samples. For example,
#                         biological replicates and corresponding controls
#                         are in the same group.
#                             1: treatment file(s)
#                             2: control file(s)
#    cores:         Number of cores for parallel processing (require parallel
#                   package).
#    debug:         Keep intermediate files (can use a lot of memory).
#                   TRUE or FALSE.
#
# Output:
#    A list of list of list of list.
#    First level is the name of the group.
#    Second level is the metadata of the group:
#        featureName:    The name of the group of regions.
#        designName:     The name of the column in the design file.
#        bamFiles:       The list of bam files associated with this combination
#                        of featureName/designName.
#    * The third level is the bamFiles, which are a list of the normalized
#      expression of every regions.
#    * There are as many elements in the list as there are regions to parse.
#    * Each bam files has a list containing every regions.
parseFeatures <- function(bamFiles, features=NULL, specie="human", maxDistance=5000, design=NULL, cores=1, debug=FALSE) {
    groups <- list()
    if (! is.null(design)) {
        groups$design <- design
    }

    # 0. Check if params are valid
    groups$param <- list()
    groups$param$specie <- specie
    groups$param$maxDistance <- maxDistance
    groups$param$cores <- cores
    groups$range <- c(-1*maxDistance, maxDistance)

    # 1. Prepare bam files
    cat("Step 1: Prepare bam files...")
    groups$bamFilesDescription <- prepareBamFiles(bamFiles, cores=cores)
    cat(" Done!\n")

    # 2. Prepare regions
    cat("Step 2: Prepare regions...")
    groups$regionsGroups <- prepareFeatures(features=features, specie=specie, maxDistance=maxDistance, cores=cores)
    cat(" Done!\n")

    # 3. Parse bam files
    cat("Step 3: Parse bam files...\n")
    groups$raw <- parseBamFiles(groups$bamFilesDescription$bam, groups$regionsGroups, cores=cores)
    groups$rpm <- rawCountsToRPM(groups$raw, groups$bamFilesDescription)
    if (debug == FALSE) groups$raw <- NULL
    groups$data <- prepareGroups(names(groups$regionsGroups), bamFiles=groups$bamFilesDescription, design=groups$design)
    groups$data <- applyOnGroups(groups=groups$data, cores=1, FUN=removeControls, data.rpm=groups$rpm, design=groups$design, controlCores=cores)
    if (debug == FALSE) groups$rpm <- NULL
    cat("Step 3: Parse bam files... Done!\n")

    # 4. Merge matrix
    cat("Step 4: Merge matrix...")
    groups$data <- applyOnGroups(groups=groups$data, cores=cores, FUN=mergeMatrix, level="noCTRL")
    if (debug == FALSE) {
        groups$data <- applyOnGroups(groups=groups$data, cores=cores, FUN=function(x) { x$noCTRL <- NULL; return(x) })
    }
    groups$matrix <- lapply(1:length(groups$data), function(x) cbind(groups$data.leftPaddings[[x]]$matrix, groups$data[[x]]$matrix, groups$data.rightPaddings[[x]]$matrix))
    names(groups$matrix) <- names(groups$data)
    if (debug == FALSE) {
        groups$data <- NULL
        groups$regionsGroups <- NULL
    }
    cat(" Done!\n")
    return(groups)
}
