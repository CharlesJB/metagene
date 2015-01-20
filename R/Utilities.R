# Created by Charles Joly Beauparlant
# 2013-11-26

# Fetch the annotation of all genes
#
# Input:
#    specie:    human: Homo sapiens (default) / mouse: Mus musculus
#
# Prerequisites:
# The specie has to be either "mouse" or "human" (default).
#
# Output:
#    A GRanges object with a ensembl_gene_id columns
getGenes <- function(specie="human") {
    # Check prerequisites

    # The specie argument has to a valid specie
    if (! specie %in% get_valid_species()){
        message <- "Incorrect parameter for specie name.\n"
        message <- paste0(message, "Currently supported species are: \"")
        message <- paste0(message, paste(get_valid_species(), collapse="\", \""))
        message <- paste0(message, "\".")

        stop(paste("Incorrect parameter for specie name.\nCurrently supported species are: \"", paste(get_valid_species(), collapse="\", \""),  "\".", collapse="", sep=""))
    }

    # Set the correct specie
    if (specie == "human") {
        load(system.file("extdata/TSS.human", package="metagene"))
        toReturn <- TSS.human
    } else {
        load(system.file("extdata/TSS.mouse", package="metagene"))
        toReturn <- TSS.mouse
    }

    # Fetch the data
    return(toReturn)
}

# Fetch the annotation of all genes from biomaRt
#
# Input:
#    specie:    human: Homo sapiens (default) / mouse: Mus musculus
#
# Prerequisites:
# The specie has to be either "mouse" or "human" (default).
#
# Output:
#    A GRanges object with a feature columns corresponding to
#    ensembl_gene_id
getGenesBiomart <- function(specie="human") {
    library(biomaRt)

    # Check prerequisites

    # The specie argument has to a valid specie
    if (! specie %in% get_valid_species()){
        message <- "Incorrect parameter for specie name.\n"
        message <- paste0(message, "Currently supported species are: \"")
        message <- paste0(message, paste(get_valid_species(), collapse="\", \""))
        message <- paste0(message, "\".")

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
    colnames(sub.ensmart) <- c("feature", "strand", "seqnames", "start", "end")
    sub.ensmart$seqnames <- paste0("chr", sub.ensmart$seqnames)

    # Center the range on the TSS
    positive <- sub.ensmart$strand == 1
    negative <- sub.ensmart$strand == -1
    sub.ensmart[positive,]$end <- sub.ensmart[positive,]$start
    sub.ensmart[negative,]$start <- sub.ensmart[negative,]$end

    return(as(sub.ensmart, "GRanges"))
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
