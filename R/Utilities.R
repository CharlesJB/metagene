# Created by Charles Joly Beauparlant
# 2013-11-26

# Sort and index bam files, if necessary. Return the number of aligned reads for
# each bam file.
#
# Input:
#    bamFiles:    Vector containing the list of every bam filename to be
#                 included in the analysis.
#    cores:       Number of cores used by the function.
#
# Prerequisites:
#     The number of cores has to be a positive integer.
#     All BAM files must exist.
#
# Output:
#    A data.frame containing the indexed bam filename and number of aligned
#    reads for each bam file.
#    Column names: bamFiles and alignedCount
prepareBamFiles <- function(bamFiles, cores = 1) {
    # Check prerequisites

    # The number of cores has to be a positive integer
    if(!is.numeric(cores) || as.integer(cores) != cores || cores <= 0) {
            stop("The number of cores has to be a positive integer.")
    }

    # All BAM file names must be of string type
    if (sum(unlist(lapply(bamFiles, is.character))) != length(bamFiles)) {
        stop("At least one BAM file name is not a valid name (a character string).")
    }

    # All BAM files must exist
     if (sum(unlist(lapply(bamFiles, file.exists))) != length(bamFiles)) {
         stop("At least one BAM file does not exist.")
     }

    # This function will only index a file if there is no index file
    indexBamFiles <- function(bamFile) {
        if (file.exists(paste(bamFile, ".bai", sep=""))  == FALSE) {
            # If there is no index file, we sort and index the current bam file
            # TODO: we need to check if the sorted file was previously produced before
            #       doing the costly sort operation
            sortedBamFile <- basename(bamFile)
            sortedBamFile <- paste(sortedBamFile, ".sorted", sep="")
            sortBam(bamFile, sortedBamFile)
            sortedBamFile <- paste(sortedBamFile, ".bam", sep="")
            indexBam(sortedBamFile)
            bamFile <- sortedBamFile
            cat(bamFile)
        }
        return(bamFile)
    }
    results <- as.data.frame(unlist(lapply(bamFiles, indexBamFiles)))
    colnames(results) <- "bam"
    results$oldBam <- bamFiles
    results$bam <- as.character(results$bam)
    results$oldBam <- as.character(results$oldBam)

    # This function will calculate the number of aligned read in each bam file
    countAlignedReads <- function(bamFile) {
        return(countBam(bamFile, param=ScanBamParam(flag = scanBamFlag(isUnmappedQuery=FALSE)))$records)
    }
    if (cores > 1) {
        library(parallel)
        results$alignedCount <- unlist(mclapply(bamFiles, countAlignedReads, mc.cores=cores))
    } else {
        results$alignedCount <- unlist(lapply(bamFiles, countAlignedReads))
    }

    return(results)
}

# Convert a list of IDs into a list of genomic regions
#
# Input:
#    features:       Either a filename of a vector of filenames.
#                    Supported features: ensembl_gene_id
#                    File must contain a header that correspond to the name
#                    of the group
#    specie:         human: Homo sapiens (default) / mouse: Mus musculus
#    maxDistance:    The distance around feature to include in the plot.
#    cores:          Number of cores for parallel processing (require
#                    parallel package).
#
# Prerequisites:
# All features files must exist.
# The specie has to be either "mouse" or "human" (default).
# The maximum distance has to be a positive integer.
# The number of cores has to be a positive integer.
#
# Output:
#    A list of data.frame. One data.frame by group of features.
#    The names of each element of the list correspond to the name of the group.
prepareFeatures <- function(features, specie="human", maxDistance=5000, cores=1) {

    # Check prerequisites

    # All features file names must be of string type
    if (!is.null(features) && sum(unlist(lapply(features, is.character))) != length(features)) {
        stop("At least one features file name is not a valid name (a character string).")
    }

    # All features files must exist
    if (!is.null(features) && sum(unlist(lapply(features, file.exists))) != length(features)) {
        stop("At least one features file does not exist.")
    }

    # The specie argument has to a valid specie
    if (! specie %in% get_valid_species()){
        stop(paste("Incorrect parameter for specie name.\nCurrently supported species are: \"", paste(get_valid_species(), collapse="\", \""),  "\".", collapse="", sep=""))
    }

    # The maximum dsitance has to be a positive integer
    if (!is.numeric(maxDistance) || maxDistance <= 0 || maxDistance %% 1 != 0) {
        stop("The maximum distance has to be a positive numeric with no decimals.")
    }

    # The number of cores has to be a positive integer
    if(!is.numeric(cores) || cores <= 0 || cores %% 1 != 0) {
        stop("The number of cores has to be a positive numeric with no decimals.")
    }

    knownGenes <- getGenes(specie)
    extractFeatures <- function(filename) {
        currentFeatures <- read.table(filename, header = TRUE)
        #currentGroup <- colnames(currentFeatures)[1]
        currentFeature <- knownGenes[knownGenes$feature %in% as.character(currentFeatures[,]),]
        # We want to return regions at +- maxDistance from starting position of current feature
        currentFeature$end_position <- currentFeature$start_position
        currentFeature$start_position <- currentFeature$start_position - maxDistance
        currentFeature$end_position <- currentFeature$end_position + maxDistance
        return(currentFeature)
    }
    featuresGroups <- list()
    if (!is.null(features)) {
        if (cores > 1) {
            library(parallel)
            featuresGroups <- mclapply(features, extractFeatures, mc.cores=cores)
        } else {
            featuresGroups <- lapply(features, extractFeatures)
        }
        names(featuresGroups) <- unlist(lapply(features, function(x) as.character(read.table(x, nrows=1)[1,])))
    } else {
        knownGenes$end_position <- knownGenes$start_position
        knownGenes$start_position <- knownGenes$start_position - maxDistance
        knownGenes$end_position <- knownGenes$end_position + maxDistance
        featuresGroups$allTSS <- knownGenes
    }
    return(featuresGroups)
}

# Fetch the annotation of all genes
#
# Input:
#    specie:    human: Homo sapiens (default) / mouse: Mus musculus
#
# Prerequisites:
# The specie has to be either "mouse" or "human" (default).
#
# Output:
#    A data.frame with 5 columns:
#        1: feature -> Ensembl gene id
#        2: strand -> -1 or 1
#        3: space -> chromosome
#        4: start_position -> position of the TSS
#        5: end_position -> ending position of the last exon of the gene
getGenes <- function(specie="human") {
    # Check prerequisites

    # The specie argument has to a valid specie
    if (! specie %in% get_valid_species()){
        stop(paste("Incorrect parameter for specie name.\nCurrently supported species are: \"", paste(get_valid_species(), collapse="\", \""),  "\".", collapse="", sep=""))
    }

    # Set the correct specie
    if (specie == "human") {
        chrom <- c(as.character(seq(1,21)),"X","Y")
        ensmart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    } else {
        chrom <- c(as.character(seq(1,19)),"X","Y")
        ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    }

    # Fetch the data
    attributes <- c("ensembl_gene_id","strand", "chromosome_name","start_position","end_position")
    filters <- c("chromosome_name")
    sub.ensmart <- getBM(attributes=attributes,filters=filters,values=chrom, mart=ensmart)
    colnames(sub.ensmart) <- c("feature", "strand", "space", "start_position", "end_position")
    sub.ensmart$space <- paste("chr", sub.ensmart$space, sep="")

    # If a gene is on the negative strand, we want the start position value to be greater than the end position
    tmp <- sub.ensmart$start_position
    sub.ensmart[sub.ensmart$strand == -1,]$start_position <- sub.ensmart[sub.ensmart$strand == -1,]$end_position
    sub.ensmart[sub.ensmart$strand == -1,]$end_position <- tmp[sub.ensmart$strand == -1]

    return(sub.ensmart)
}


# Get list of currently supported species
#
# Input:
#   None
#
# Output:
#   List of currently supported species
get_valid_species<-function() {
	# Return list of valid species
	return(c("mouse", "human"))
}

