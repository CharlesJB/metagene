# Created by Charles Joly Beauparlant
# 2013-11-26

# Parse an experiment using regions that can be of different length.
#
# Input:
#    regions:        A vector of bed file names corresponding to the regions to
#                    include in the analysis.
#                    The file name (minus the extension) will be used as the
#                    name of the region.
#    bamFiles:       A vector of bamFile to plot. TODO: Should also accept a
#                    list of bamfiles where each elements would be grouped
#                    together
#    specie:         human: Homo sapiens (default) / mouse: Mus musculus
#    design:         A data.frame explaining the relationship between
#                    multiple samples.
#                        * One line per samples.
#                        * One column per group of samples. For example,
#                          biological replicates and corresponding controls are
#                          in the same group.
#                            1: treatment file(s)
#                            2: control file(s)
#    paddingSize:    The length padding we want to add on each side of each
#                    regions.
#                    The padding is added after the scaling step.
#    cores:          Number of cores for parallel processing (require
#                    parallel package).
#    debug:          Keep intermediate files (can use a lot of memory).
#                    TRUE or FALSE.
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
parseRegions <- function(regions, bamFiles, specie="human", design=NULL, paddingSize=2000, cores=1, debug=FALSE) {
    groups <- list()
    if (! is.null(design)) {
        groups$design <- design
    }

    # 0. Check if params are valid
    groups$param <- list()
    groups$param$specie <- specie
    groups$param$paddingSize <- paddingSize
    groups$param$cores <- cores
    groups$range <- c(-1,1) # TODO: find something better

    # 1. Prepare bam files
    cat("Step 1: Prepare bam files...")
    groups$bamFilesDescription <- prepareBamFiles(bamFiles, cores=cores)
    cat(" Done!\n")

    # 2. Prepare regions
    cat("Step 2: Prepare regions...")
    #regionsGroups <- prepareRegions(regions, cores=cores)
    groups$regionsGroups <- prepareRegions(regions, cores=cores)
    groups$regionsGroups.leftPaddings <- prepareRegionsPaddings(groups$regionsGroups, side="left", paddingSize=paddingSize, cores=cores)
    groups$regionsGroups.rightPaddings <- prepareRegionsPaddings(groups$regionsGroups, side="right", paddingSize=paddingSize, cores=cores)
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

    # 4. Scale vectors
    cat("Step 4: Scale vectors...")
    # Get the length of every bamFiles of a group
    getLengthsGroup <- function(group) {
        getLengthsBamFile <- function(bamFile) {
            return(sapply(bamFile, length))
        }
        return(unlist(lapply(group$noCTRL, getLengthsBamFile)))
    }
    # Get the median length of every regions groups in the current experiment
    medianLength <- median(unlist(applyOnGroups(groups$data, cores=cores, getLengthsGroup)))
    groups$data <- applyOnGroups(groups=groups$data, cores=1, FUN=scaleVectors, domain=medianLength, level="noCTRL", scaleCores=cores)
    if (debug == FALSE) {
        groups$data <- applyOnGroups(groups=groups$data, cores=cores, FUN=function(x) { x$noCTRL <- NULL; return(x) })
    }
    cat(" Done!\n")

    # 5. Parse padding
    cat("Step 5: Parse padding...\n")
    # TODO: put this in a function to call only once for step3, once for left and once for right.
    groups$raw.leftPaddings <- parseBamFiles(groups$bamFilesDescription$bam, groups$regionsGroups.leftPaddings, cores=cores)
    groups$raw.rightPaddings <- parseBamFiles(groups$bamFilesDescription$bam, groups$regionsGroups.rightPaddings, cores=cores)
    groups$rpm.leftPaddings <- rawCountsToRPM(groups$raw.leftPaddings, groups$bamFilesDescription)
    if (debug == FALSE) groups$raw.leftPaddings <- NULL
    groups$rpm.rightPaddings <- rawCountsToRPM(groups$raw.rightPaddings, groups$bamFilesDescription)
    if (debug == FALSE) groups$raw.rightPaddings <- NULL
    groups$data.leftPaddings <- prepareGroups(names(groups$regionsGroups.leftPaddings), bamFiles=groups$bamFilesDescription, design=groups$design)
    groups$data.rightPaddings <- prepareGroups(names(groups$regionsGroups.rightPaddings), bamFiles=groups$bamFilesDescription, design=groups$design)
    groups$data.leftPaddings <- applyOnGroups(groups=groups$data.leftPaddings, cores=1, FUN=removeControls, data.rpm=groups$rpm.leftPaddings, design=groups$design, controlCores=cores)
    if (debug == FALSE) groups$rpm.leftPaddings <- NULL
    groups$data.rightPaddings <- applyOnGroups(groups=groups$data.rightPaddings, cores=1, FUN=removeControls, data.rpm=groups$rpm.rightPaddings, design=groups$design, controlCores=cores)
    if (debug == FALSE) groups$rpm.rightPaddings <- NULL
    cat(" Done!\n")

    # 6. Merge matrix
    cat("Step 6: Merge matrix...")
    groups$data <- applyOnGroups(groups=groups$data, cores=cores, FUN=mergeMatrix, level="scaled")
    groups$data.leftPaddings <- applyOnGroups(groups=groups$data.leftPaddings, cores=cores, FUN=mergeMatrix, level="noCTRL")
    groups$data.rightPaddings <- applyOnGroups(groups=groups$data.rightPaddings, cores=cores, FUN=mergeMatrix, level="noCTRL")
    if (debug == FALSE) {
        groups$data <- applyOnGroups(groups=groups$data, cores=cores, FUN=function(x) { x$scaled <- NULL; return(x) })
        groups$data.leftPaddings <- applyOnGroups(groups=groups$data.leftPaddings, cores=cores, FUN=function(x) { x$noCTRL <- NULL; return(x) })
        groups$data.rightPaddings <- applyOnGroups(groups=groups$data.rightPaddings, cores=cores, FUN=function(x) { x$noCTRL <- NULL; return(x) })
    }
    groups$matrix <- lapply(1:length(groups$data), function(x) cbind(groups$data.leftPaddings[[x]]$matrix, groups$data[[x]]$matrix, groups$data.rightPaddings[[x]]$matrix))
    names(groups$matrix) <- names(groups$data)
    if (debug == FALSE) {
        groups$data <- NULL
        groups$data.leftPaddings <- NULL
        groups$data.rightPaddings <- NULL
        groups$regionsGroups <- NULL
        groups$regionsGroups.leftPaddings <- NULL
        groups$regionsGroups.rightPaddings <- NULL
    }
    cat(" Done!\n")
    return(groups)
}
