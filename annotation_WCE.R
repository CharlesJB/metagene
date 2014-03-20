# Created by Charles Joly Beauparlant
# 2013-11-25

# Fetch the annotation of all mouse genes
#
# INPUT:
#
# OUTPUT:
#	A data.frame with 5 columns:
#		1: chr
#		2: start
#		3: end
#		4: gene name (ENSMUG)
#		5: strand
getMouseGenes <- function(maxDistance) {
	suppressMessages(library(ChIPpeakAnno))
	suppressMessages(library(org.Mm.eg.db))
	data(TSS.mouse.NCBIM37)
	selected <- c(as.character(1:19), "X", "Y")
	anno <- as.data.frame(TSS.mouse.NCBIM37[TSS.mouse.NCBIM37$space %in% selected,])
	anno <- anno[,c(1,2,3,5,6)]
	anno$space <- as.character(anno$space)
#	anno$space <- paste("chr", anno$space, sep="")
	anno$distancetoFeature <- 0
	colnames(anno)[4] <- "feature"
	anno$start_position <- anno$start
	anno$end <- anno$start + maxDistance
	anno$start <- anno$start - maxDistance
	print(maxDistance)
	return(anno)
}

# Add gene annotation to a bed file. Currently hard coded for mus musculus data.
#
# INPUT:
#	BEDfile: Path to the bed file to load.
#
# OUTPUT:
#	A data.frame with the annotation obtained with ChIPpeakAnnoPackage.
annotateBED <- function(BEDfile) {
	suppressMessages(library(ChIPpeakAnno))
	suppressMessages(library(org.Mm.eg.db))
	data(TSS.mouse.NCBIM37)
	myPeaks <- BED2RangedData(BEDfile)
	annotatedPeaks <- annotatePeakInBatch(myPeaks, AnnotationData=TSS.mouse.NCBIM37)
	IDs_df <- as.data.frame(addGeneIDs(annotatedPeaks, "org.Mm.eg.db", IDs2Add="ensembltrans"))
	IDs_df$space <- paste("chr", IDs_df$space, sep="")
	return(as.data.frame(IDs_df))
}

# Convert a ensembl_transcript_id to ensembl_gene_id. 
#
# INPUT: 
#	TranscriptVector: A vector with the ensembl_transcript_id to convert
#
# OUPUT:
#	A data.frame with 2 columns:
#		* ensembl_transcript_id
#		* ensembl_gene_id
convertTranscriptIDtoEnsemblID <- function(TranscriptVector) {
	suppressMessages(library(biomaRt))
	ensembl <- useMart("ensembl")
	ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
	return(getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id'), filters='ensembl_transcript_id', values=TranscriptVector, mart=ensembl))
}

# This function will load all the ensembl_transcript_id from every groups.
#
# INPUT:
#	group_filenames: A vector with all the file names for the groups.
#		Each file has a header, which correspond to the unique name of the group.
#		Each file contains a single column with the ensembl_transcript_id of the genes.
#
# OUPUT:
#	A data.frame with 3 or more columns:
#		* ensembl_transcript_id
#		* ensembl_gene_id
#		* a column for each group with a boolean value for each ensembl_transcript_id
#
#	Note: there is as many columns as there are unique ensembl_transcript_id among every groups.
loadGroups <- function(groups) {
	# 0. Load group filenames
	group_filenames <- read.table(groups, stringsAsFactors=FALSE)[,1]
	# 1. Get a list of all unique id
	list_unique_id <- vector()
	for (filename in group_filenames) {
		current_group <- read.table(filename, header=TRUE)
		list_unique_id <- unique(c(list_unique_id, as.vector(current_group[[1]])))
	}
	# 2. Obtain the ensembl_id for each id
#	list_unique_id <- convertTranscriptIDtoEnsemblID(list_unique_id)
	list_unique_id <- as.data.frame(list_unique_id)
	colnames(list_unique_id) <- c("ensembl_gene_id")
	# 3. Add the boolean for each groups
	for (filename in group_filenames) {
		base_name <- basename(sub("^([^.]*).*", "\\1", filename))
		current_group <- read.table(filename, header=TRUE)
		list_unique_id[,base_name] <- list_unique_id[[1]] %in% current_group[[1]]
	}
	return(list_unique_id)
}
