# Created by Charles Joly Beauparlant
# 2013-11-26

# Create a metagene plot based on a list of genomic features
#
# Input:
#	bamFiles:	A vector of bamFile to plot. TODO: Should also accept a list of bamfiles where each elements would be grouped together
#	features:	Either a filename of a vector of filenames.
#			Supported features: ensembl_gene_id
#			File must contain a header that correspond to the name of the group
#	specie:		human: Homo sapiens (default) / mouse: Mus musculus
#	maxDistance:	The distance around feature to include in the plot.
#	design:		A matrix explaining the relationship between multiple samples.
#			One line per samples.
#			One column per group of samples. For example, biological replicates and corresponding controls are in the same group.
#			1: treatment file(s)
#			2: control file(s)
#	binSize:	The number of nucleotides in each bin for the bootstrap step.
#	alpha: 		Confidence interval.
# 	sampleSize:	Number of time each bin will be resampled (should be at least 1000).
#	cores:		Number of cores for parallel processing (require parallel package).
# TODO: Add group of bam files in design file
plotFeatures <- function(bamFiles, features=NULL, specie="human", maxDistance=5000, design=NULL, binSize=100, alpha=0.05, sampleSize=1000, cores=1) {
	# 0. Check if params are valid

	# 1. Prepare bam files
	cat("Step 1: Prepare bam files...")
	bamFiles <- prepareBamFiles(bamFiles, cores=cores)
	cat(" Done!\n")

	# 2. Prepare regions
	cat("Step 2: Prepare regions...")
	featuresGroups <- prepareFeatures(features=features, specie=specie, cores=cores)
	cat(" Done!\n")

	# 3. Parse bam files
	cat("Step 3: Parse bam files...\n")
	groupNames <- prepareGroups(names(featuresGroups), bamFiles=bamFiles, design=design)
	groups <- parseBamFiles.old(bamFiles, featuresGroups, groups=groupNames, design=design, cores=cores)
	cat("Step 3: Parse bam files... Done!\n")

	# 4. Merge matrix
	cat("Step 4: Merge matrix...")
	groups <- applyOnGroups(groups=groups, cores=cores, FUN=mergeMatrix.old)
	cat(" Done!\n")

	# 5. Bootstrap
	cat("Step 5: Bootstrap...")
	bootstrap <- function(x) {
		x$bootstrap <- bootstrapAnalysis.old(x, binSize=binSize, alpha=alpha, sampleSize=sampleSize, cores=cores)
		return(x)
	}
	groups <- applyOnGroups(groups=groups, cores=1, FUN=bootstrap)
	cat(" Done!\n")

	# TODO: Add previous params to plotFeatures function
	# 6. Plot
	cat("Step 6: Plot...")
	groups$graphData <- plot.getDataFrame.old(groups)
	plot.graphic(groups$graphData, paste(names(groups)))
	cat(" Done!\n")
	return(groups)
}

# Convert list of vectors in matrix for a single group
#
# Input:
#	group:	Correspond to an element in the data subsection of the main data structure.
#	level: 	The names of the element to merge.
#
# Output:
#	The group in input with an extra element named matrix (group$matrix).
mergeMatrix <- function(group, level) {
	mergeMatrixFixedLengthFeatures <- function(currentGroup) {
		#currentBamFiles <- currentGroup$bamFiles
		#return(lapply(1:length(currentBamFiles), function(x) currentBamFiles[[group]] <- do.call(rbind, currentBamFiles[[group]])))
		#return(lapply(1:length(currentBamFiles), function(x) currentGroup$noCTRL[[group]] <- do.call(rbind, currentGroup$noCTRL[[group]])))
		return(lapply(currentGroup[[level]], function(x) do.call(rbind, x)))
	}
	group$matrix <- do.call(rbind, mergeMatrixFixedLengthFeatures(group))
	return(group)
}

# Convert list of vectors in matrix for a single group
#
# Input:
#	group:	Correspond to an element in the data subsection of the main data structure.
#
# Output:
#	The group in input with an extra element named matrix (group$matrix).
mergeMatrix.old <- function(group) {
       mergeMatrixFixedLengthFeatures <- function(currentExperiment) {
               currentBamFiles <- currentExperiment$bamFiles
               return(lapply(1:length(currentBamFiles), function(x) currentBamFiles[[group]] <- do.call(rbind, currentBamFiles[[group]])))
        }
	group$matrix <- do.call(rbind, mergeMatrixFixedLengthFeatures(group))
        return(group)
}

# Parse an experiment using regions that can be of different length.
#
# Input:
#	regions:	A vector of bed file names corresponding to the regions to include in the analysis.
#			The file name (minus the extension) will be used as the name of the region.
#	bamFiles:	A vector of bamFile to plot. TODO: Should also accept a list of bamfiles where each elements would be grouped together
#	specie:		human: Homo sapiens (default) / mouse: Mus musculus
#	design:		A matrix explaining the relationship between multiple samples.
#			One line per samples.
#			One column per group of samples. For example, biological replicates and corresponding controls are in the same group.
#			1: treatment file(s)
#			2: control file(s)
#	paddingSize:	The length padding we want to add on each side of each regions.
#			The padding is added after the scaling step.
#	cores:		Number of cores for parallel processing (require parallel package).
#	debug:		Keep intermediate files (can use a lot of memory). TRUE or FALSE.
#
# Output:
#	A list of list of list of list.
#	First level is the name of the group.
#	Second level is the metadata of the group:
#		featureName:	The name of the group of regions.
#		designName:	The name of the column in the design file.
#		bamFiles:	The list of bam files associated with this combination of featureName/designName
#	The third level is the bamFiles, which are a list of the normalized expression of every regions.
#		There are as many elements in the list as there are regions to parse.
#		Each bam files has a list containing every regions.
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

	# 1. Prepare bam files
	cat("Step 1: Prepare bam files...")
	#bamFiles <- prepareBamFiles(bamFiles, cores=cores)
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
	groups$data <- applyOnGroups(groups=groups$data, cores=1, FUN=scaleVectors, medianLength=medianLength, scaleCores=cores)
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

# Resize the vectors of every group so they have the same length
# Note: For this function to work correctly, we must have substracted and removed the controls.
#	This is because the controls can be redundant between groups, so the median could be influenced.
scaleVectors <- function(group, medianLength, scaleCores=1) {
	scaleBamFile <- function(bamFile, domain, scaleCores=1) {
		if (scaleCores > 1) {
			return(mclapply(bamFile, scaleVector, domain, mc.cores=scaleCores))
		} else {
			return(lapply(bamFile, scaleVector, domain))
		}
	}
	#group$bamFilesBeforeScaling <- group$bamFiles
	group$scaled <- lapply(group$noCTRL, scaleBamFile, domain=medianLength, scaleCores=scaleCores)
	names(group$bamFiles) <- names(group$bamFilesBeforeScaling)
	return(group)
}

# Apply a function on every groups of the main data structure
#
# Input:
#	groups:		The main data structure
#	cores:		Number of cores for parallel processing (require parallel package).
#	FUN:		The function to apply on every groups
#	...:		Extract arguments for FUN
#
# Output:
#	The result of the function applied
applyOnGroups <- function(groups, cores=1, FUN, ...) {
	if (cores > 1) {
		library(parallel)
		return(mclapply(groups, function(x) FUN(x, ...), mc.cores=cores))
	} else {
		return(lapply(groups, function(x) FUN(x, ...)))
	}
}

# Sort and index bam files, if necessary. Return the number of aligned reads for each bam file.
#
# Input:
#	bamFiles:	Vector containing the list of every bam filename to be included in the analysis.
# 	cores: 		Number of cores used by the function.
#
# Prerequisites:
# The number of cores has to be a positive integer.
# All BAM files must exist.
#
# Output:
#	A data.frame containing the indexed bam filename and number of aligned reads for each bam file.
#	Column names: bamFiles and alignedCount
prepareBamFiles <- function(bamFiles, cores = 1) {
	library(Rsamtools)

	# Check prerequisites
	
   	# The number of cores has to be a positive integer
	if(!is.numeric(cores) || cores <= 0) {
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
			sortedBamFile <- sub("^([^.]*).*", "\\1", bamFile)
			sortedBamFile <- paste(sortedBamFile, ".sorted", sep="")
			cat(sortedBamFile)
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

# Fetch the annotation of all genes
#
# INPUT:
#	
#	specie:		human: Homo sapiens (default) / mouse: Mus musculus
#
# Prerequisites:
# The specie has to be either "mouse" or "human" (default).
#
# OUTPUT:
#	A data.frame with 5 columns:
#		1: feature -> Ensembl gene id
#		2: strand -> -1 or 1
#		3: space -> chromosome
#		4: start_position -> position of the TSS
#		5: end_position -> ending position of the last exon of the gene
getGenes <- function(specie="human") {
	require(biomaRt)
	
	# Check prerequisites
	
	# The specie argument has only two valid possibilities
	if (! specie %in% c("mouse", "human")){
		stop("Incorrect parameter for specie name.\nCurrently supported species are \"human\" and \"mouse\".")	
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

# Convert a list of IDs into a list of genomic regions
#
# Input:
#	features:	Either a filename of a vector of filenames.
#			Supported features: ensembl_gene_id
#			File must contain a header that correspond to the name of the group
#	specie:		human: Homo sapiens (default) / mouse: Mus musculus
#	maxDistance:	The distance around feature to include in the plot.
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Prerequisites:
# All features files must exist.
# The specie has to be either "mouse" or "human" (default).
# The maximum distance has to be a positive integer.
# The number of cores has to be a positive integer.
#
# Output:
#	A list of data.frame. One data.frame by group of features.
#	The names of each element of the list correspond to the name of the group.
prepareFeatures <- function(features, specie="human", maxDistance=5000, cores=1) {
	
	# Check prerequisites
	
	# All features file names must be of string type 
	if (sum(unlist(lapply(features, is.character))) != length(features)) {
		stop("At least one features file name is not a valid name (a character string).")
	}
	
	# All festures files must exist
	if (sum(unlist(lapply(features, file.exists))) != length(features)) {
		stop("At least one features file does not exist.")
	}
	
	# The specie argument has only two valid possibilities
	if (! specie %in% c("mouse", "human")){
		stop("Incorrect parameter for specie name.\nCurrently supported species are \"human\" and \"mouse\".")	
	}
	
	# The maximum dsitance has to be a positive integer
	if(!is.integer(maxDistance) || maxDistance <= 0) {
		stop("The maximum distance has to be a positive integer.")
	}	
	
	# The number of cores has to be a positive integer
	if(!is.integer(cores) || cores <= 0) {
		stop("The number of cores has to be a positive integer.")
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
	if (!is.null(features)) {
		if (cores > 1) {
			library(parallel)
			featuresGroups <- mclapply(features, extractFeatures, mc.cores=cores)
		} else {
			featuresGroups <- lapply(features, extractFeatures)
		}
		names(featuresGroups) <- unlist(lapply(features, function(x) as.character(read.table(x, nrow=1)[1,])))
	} else {
		knownGenes$end_position <- knownGenes$start_position
		knownGenes$start_position <- knownGenes$start_position - maxDistance
		knownGenes$end_position <- knownGenes$end_position + maxDistance
		featuresGroups$allTSS <- knownGenes
	}
	return(featuresGroups)
}

# Parse bed files and convert them in a list of data.frames
#	regions:	A vector of bed file names corresponding to the regions to include in the analysis.
#			The file name (minus the extension) will be used as the name of the region.
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list of data.frame. One data.frame by group of features.
#	The names of each element of the list correspond to the name of the group.
prepareRegions <- function(regions, cores=1) {
	# 1. Parse the bed files
	readBedFile <- function(bedFileName) {
		currentBed <- read.table(bedFileName, header=FALSE)
		# We only need the infos of the first three columns
		colnames(currentBed)[1:3] <- c("space", "start_position", "end_position")
		# TODO: should we check is start_position is always smaller than end_position?
		return(currentBed)
	}
	toReturn <- applyOnGroups(regions, cores=cores, FUN=readBedFile)

	# 2. Add the names
	getBaseName <- function(bedFileName) {
		return(sub("^([^.]*).*", "\\1", basename(bedFileName)))
	}
	names(toReturn) <- unlist(applyOnGroups(regions, cores=cores, FUN=getBaseName))

	return(toReturn)
}

# Get the genomic position of the paddings
#
# Input:
#	regionsGroups:	The output of prepareRegions or prepareFeatures functions.
#	side:		The side of the padding to prepare. Either 'left' or 'right'.
#	paddingSize:	The length padding we want to add on each side of each regions.
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list of data.frames, each data.frame corresponding to the paddings for a group
#	of features.
prepareRegionsPaddings <- function(regionsGroups, side, paddingSize=2000, cores=1) {
	# Check prerequisites

	# Side must be either 'left' or 'right'
	if (! side %in% c("left", "right")) {
		stop("Side argument must be either 'left' or 'right'.")
	}

	# The number of cores has to be a positive integer
	if(!is.numeric(cores) || cores <= 0) {
		stop("The number of cores has to be a positive integer.")
	}

	getPaddings <- function(currentRegions) {
		currentPaddings <- currentRegions
		if (side == "left") {
			currentPaddings$start_position <- currentRegions$start_position - paddingSize
			currentPaddings$end_position <- currentRegions$start_position
			currentPaddings$start_position[currentPaddings$start_position < 0] <- 0
		} else {
			currentPaddings$start_position <- currentRegions$end_position
			currentPaddings$end_position <- currentRegions$end_position + paddingSize
		}
		return(currentPaddings)
	}

	toReturn <- applyOnGroups(regionsGroups, cores=cores, FUN=getPaddings)
	names(toReturn) <- names(regionsGroups)
	return(toReturn)
}


# Distribute the bam filenames in their respective groups.
#
# Input:
#	featuresGroupsNames:	The name of the group of features in the current analysis.
#	bamFiles:		A vector of bam filename(s).
#	design:			A matrix explaining the relationship between multiple samples.
#				One line per samples.
#				One column per group of samples. For example, biological replicates and corresponding controls are in the same group.
#				1: treatment file(s)
#				2: control file(s)
#
# Output:
#	A list of list of list.
#	First level is the name of the group.
#	Second level is the metadata of the group:
#		featureName:	The name of the group of features.
#		designName:	The name of the column in the design file.
#		bamFiles:	The list of bam files associated with this combination of featureName/designName
#	Third level is the actual list of bam files.
# 	Note1: If there are groups of feature and groups of bam file, every possible combination of feature groups / bam groups will be created.
#              Otherwise, the is a group for each group of features.
# 	Note2: The names of the elements the bam list are the same as the ones in the design file, which are also the original bam file names.
#	       The content of the each element of the list of bam file is the sorted bam file name, which should be use when parsing bam.
prepareGroups <- function(featuresGroupsNames, bamFiles, design=NULL) {
	allGroups <- list()
	for (featureName in featuresGroupsNames) {
		# If there is no design file, we put every bamFile in each groups
		if (is.null(design)) {
			allGroups[[featureName]] <- list()
			allGroups[[featureName]]$featureName <- featureName
			allGroups[[featureName]]$bamFiles <- as.list(as.character(bamFiles$bam))
			names(allGroups[[featureName]]$bamFiles) <- as.character(bamFiles$oldBam)
		} else {
			for (i in 2:ncol(design)) {
				name <- paste(featureName, colnames(design)[i], sep="_")
				bam <- design[design[,i] != 0,][,1]
				designName <- colnames(design)[i]
				# We want the sorted name of the bam files as the content of the list
				j <- numeric()
				for (bamFile in bam) {
					j <- append(j, which(bamFiles$oldBam == bamFile))
				}
				sortedBamNames <- as.character(bamFiles[j,]$bam)
				# Create the actual element of the list
				allGroups[[name]] <- list()
				allGroups[[name]]$featureName <- featureName
				allGroups[[name]]$designName <- designName
				allGroups[[name]]$bamFiles <- as.list(sortedBamNames)
				names(allGroups[[name]]$bamFiles) <- bam
			}
		}
	}
	return(allGroups)
}

# Parse multiple bamFiles
#
# Input:
#	bamFiles:	The data.frame obtained with the prepareBamFiles function.
#	featuresGroups:	A list of data.frame. One data.frame by group of features.
#			The names of each element of the list correspond to the name of the group.
#	design:		A matrix explaining the relationship between multiple samples.
#			If there is a design, the controls will automatically be substracted from every
#			replicates of the input files.
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list of list of list of list.
#	First level is the name of the group.
#	Second level is the metadata of the group:
#		featureName:	The name of the group of features.
#		designName:	The name of the column in the design file.
#		bamFiles:	The list of bam files associated with this combination of featureName/designName
#	The third level is the bamFiles, which are a list of the normalized expression of every features.
#		There are as many elements in the list as there are features to parse.
#		Each bam files has a list containing every features.
parseBamFiles.old <- function(bamFiles, featuresGroups, groups, design=NULL, cores=1) {
	parseGroup <- function(groupName, featuresGroups) {
		print(paste("parseGroup:", groupName))
		# Get bam files
		currentBamFiles <- groups[[groupName]]$bamFiles

		# Get features
		currentFeatureName <- groups[[groupName]]$featureName
		currentFeatures <- featuresGroups[[currentFeatureName]]

		# Parse each bam
		currentBamFiles <- lapply(currentBamFiles, function(x) parseBamFile.old(x, bamFiles[bamFiles$bam == x,]$alignedCount, currentFeatures, cores))
		# Add the names
		currentGroup <- list()
		currentGroup$featureName <- groups[[groupName]]$featureName
		currentGroup$designName <- groups[[groupName]]$designName
		currentGroup$bamFiles <- currentBamFiles
		names(currentGroup$bamFiles) <- names(groups[[groupName]]$bamFiles)
		currentGroup$bamFilesWithCTRL <- currentGroup$bamFiles
		return(currentGroup)
	}
	# 1 Parse every bam file individually
	groupNames <- names(groups)
	groups <- lapply(groupNames, function(x) parseGroup(x, featuresGroups))
	names(groups) <- groupNames
	# 2 Remove control
	if (!is.null(design)) {
		groups <- applyOnGroups(groups=groups, cores=cores, FUN=removeControls.old, design=design)
	}
	return(groups)
}

# Parse multiple bamFiles
#
# Input:
#	bamFiles:	A vector of bam files
#	featuresGroups:	A list of data.frame. One data.frame by group of features.
#			The names of each element of the list correspond to the name of the group.
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list of list of list that contains the raw counts for every features groups:
#	First level:	One entry by features group
#	Second level:	One entry per bamFile in current features group
#	Third level:	One entry per feature in current features group
parseBamFiles <- function(bamFiles, featuresGroups, cores=1) {
	parseFeatures <- function(features, bamFiles, cores) {
		raw.counts <- lapply(bamFiles, parseBamFile, features=features, cores=cores)
		names(raw.counts) <- bamFiles
		return(raw.counts)
	}
	raw.data <- lapply(featuresGroups, parseFeatures, bamFiles=bamFiles, cores=cores)
	names(raw.data) <- names(featuresGroups)
	return(raw.data)
}

# Convert the raw counts into reads per million aligned (rpm)
#
# Input:
#	rawCounts:		The data structure returned by the parseBamFiles function.
#	bamFilesDescription:	The data.frame obtained with the prepareBamFiles function.
#	cores:			Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list of list of list that contains the rpm for every features groups:
#	First level:	One entry by features group
#	Second level:	One entry per bamFile in current features group
#	Third level:	One entry per feature in current features group
rawCountsToRPM <- function(rawCounts, bamFilesDescription, cores=1) {
	convertToRPM <- function(bamFile, bamFilesDescription, cores) {
		alignedCount <- bamFilesDescription[bamFilesDescription$bam == names(bamFile),]$alignedCount
		bamFile <- bamFile[[1]]
		return(applyOnGroups(bamFile, cores, FUN=function(x) { x / (alignedCount / 1000000) }))
	}
	rpmCounts <- list()
	for (i in 1:length(rawCounts)) {
		currentFeaturesGroupName <- names(rawCounts)[i]
		featuresGroup <- rawCounts[[i]]
		rpmCounts[[currentFeaturesGroupName]] <-
			lapply(1:length(featuresGroup), function(x) convertToRPM(featuresGroup[x], bamFilesDescription=bamFilesDescription, cores=cores))
		names(rpmCounts[[currentFeaturesGroupName]]) <- names(featuresGroup)
	}
	return(rpmCounts)
}

# Substract controls from a single group
#
# Input:
#	group:		A group extracted from the main data structure
#	currentDesign:	The line from matrix explaining the relationship between current samples.
#
# Output:
#	The group extracted from the main data structure from which the controls were substracted then deleted
removeControls.old <- function(group, design) {
	# 1. Extract relevant design columns
	currentDesign <- design[,c(1,which(colnames(design) == group$designName))]
	treatmentNames <- as.character(currentDesign[currentDesign[,2] == 1, 1])
	controlNames <- as.character(currentDesign[currentDesign[,2] == 2, 1])

	# 2. Merge controls
	# TODO
	if (length(controlNames) > 1) {
		print("TODO")
	} else {
		mergedControls <- group$bamFiles[[controlNames[1]]]
		group$bamFiles[[controlNames[1]]] <- NULL
	}

	# 3. Substract merged control from every treatment samples
	substractFeature <- function(treatment, i) {
		currentFeatureTreatment <- unlist(treatment[i])
		currentFeatureControl <- unlist(mergedControls[i])
		newValues <- currentFeatureTreatment - currentFeatureControl
		newValues[newValues < 0] <- 0
		return(newValues)
	}
	#substractControl <- function(treatmentName) {
	substractControl <- function(currentTreatment) {
		newTreatment <- lapply(1:length(currentTreatment), function(x) substractFeature(currentTreatment, x))
		return(newTreatment)
	}
	group$bamFiles <- lapply(group$bamFiles, substractControl)
	names(group$bamFiles) <- treatmentNames
	return(group)
}

# Substract controls from a single group
#
# Input:
#	group:		A group extracted from the main data structure
#	currentDesign:	The line from matrix explaining the relationship between current samples.
#
# Output:
#	The group extracted from the main data structure from which the controls were substracted then deleted
removeControls <- function(group, data.rpm, design=NULL, controlCores=1) {
	if (!is.null(design)) {
		designName <- group$designName

		# 1. Extract relevant design columns
		currentDesign <- design[,c(1,which(colnames(design) == designName))]
		treatmentNames <- as.character(currentDesign[currentDesign[,2] == 1, 1])
		controlNames <- as.character(currentDesign[currentDesign[,2] == 2, 1])

		# 2 Get the relevant group of rpm
		current.rpm <- data.rpm[[group$featureName]]

		# 3. Merge controls
		# TODO
		if (length(controlNames) > 1) {
			print("TODO")
			# uuu <- data.frame(zzz <- c(1,2,3), xxx <- c(4,5,6), yyy <- c(7,8,9))
		} else {
			mergedControls <- current.rpm[[controlNames[1]]]
		}

		# 4. Substract merged control from every treatment samples
		substractFeature <- function(treatment, i) {
			currentFeatureTreatment <- unlist(treatment[i])
			currentFeatureControl <- unlist(mergedControls[i])
			newValues <- currentFeatureTreatment - currentFeatureControl
			newValues[newValues < 0] <- 0
			return(newValues)
		}
		substractControl <- function(currentTreatment) {
			newTreatment <- applyOnGroups(1:length(currentTreatment), cores=controlCores, FUN=function(x) substractFeature(currentTreatment, x))
			return(newTreatment)
		}
		treatment.rpm <- current.rpm[names(current.rpm) %in% treatmentNames]
		group$noCTRL <- lapply(treatment.rpm, substractControl)
		names(group$noCTRL) <- treatmentNames
		return(group)
	}
	return(group)
}

# Parse a single bam file
#
# Input:
#	bamFile:	The name of the bam file to parse. Must be sorted and indexed.
#	alignedCount:	The number of aligned reads for the current sample.
# 			TODO: if NULL, it should be calculated in the function.
#	features:	A data.frame with the infos for every features to parse
#			Must have the folowing columns: feature, strand, space, start_position and end_position
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list with an element for every feature to parse. Each element contains a vector of normalized reads expression.
parseBamFile.old <- function(bamFile, alignedCount, features, cores=1) {
	print(paste("Current bam:", bamFile))
	extractReadsDensity <- function(feature, bamFile) {
		# Extract raw counts
		currentReads <- extractReadsInRegion(bamFile, feature$space, feature$start_position, feature$end_position)
		vectorResult <- convertReadsToDensity(currentReads, feature)

		# If on negative strand, invert the current vector
		#if (feature$strand == "-1" |  feature$strand == -1 | feature$strand == "-") {
		if (feature$start > feature$end) {
			vectorResult <- rev(vectorResult)
		}

		# Convert to RPM
		vectorResult <- vectorResult / (alignedCount / 1000000)
		return(vectorResult)

	}
	if (cores > 1) {
		library(parallel)
		return(mclapply(1:nrow(features), function(x) extractReadsDensity(features[x,],  bamFile=as.character(bamFile)), mc.cores=cores))
	} else {
		return(lapply(1:nrow(features), function(x) extractReadsDensity(features[x,],  bamFile=as.character(bamFile))))
	}
}

parseBamFile <- function(bamFile, features, cores=1) {
	print(paste("Current bam:", bamFile))
	extractReadsDensity <- function(feature, bamFile) {
		# Extract raw counts
		currentReads <- extractReadsInRegion(bamFile, feature$space, feature$start_position, feature$end_position)
		vectorResult <- convertReadsToDensity(currentReads, feature)

		# If on negative strand, invert the current vector
		if (feature$start > feature$end) {
			vectorResult <- rev(vectorResult)
		}
		return(vectorResult)
	}
	return(applyOnGroups(1:nrow(features), cores=cores, FUN=function(x) extractReadsDensity(features[x,], bamFile=bamFile)))
}

# Extract reads from BAM file that overlap with a specified genomic region.
#
# INPUT:
#	bamFile:	Path to the bam file.
#	chr:		Current chromosome.
#	start:		Starting position of the current region.
#	end:		Ending position of the current region.
#
# Prerequisites:
# The BAM file must exist.
# The chromosome name must be in character format.
# The starting position has to be a positive integer.
# The ending position has to be a positive integer.
#
# OUTPUT:
#	A data.frame containing every reads overlapping the current genomic region:
#		* rname
#		* pos
#		* qwidth
extractReadsInRegion <- function(bamFile, chr, start, end) {
	
	# Check prerequisites
	
	# The BAM file name must be of string type 
	if (!is.character(bamFile)) {
		stop("The BAM file name is not a valid name (a character string).")
	}
	
	# The BAM file name must exist
	if (! file.exists(bamFile)) {
		stop("The BAM file does not exist.")
	}
	
	# TODO : VALIDER SI A FAIRE OU NON
	#The chromosome name must be of string type 
	#if (!is.character(chr)) {
	#	stop("The chromosome name is not a valid name (a character string).")
	#}
		
	# The starting position has to be a positive integer
	if(!is.numeric(start) || start <= 0) {
		stop("The starting position has to be a positive integer.")
	}	
	
	# The ending position has to be a positive integer
	if(!is.numeric(end) || end <= 0) {
		stop("The ending position has to be a positive integer.")
	}
	
	suppressMessages(library(Rsamtools))
	if (start > end) {
		tmp <- start
		start <- end
		end <- tmp
	}
	df <- data.frame()
	which <- GRanges(seqnames=Rle(chr), ranges=IRanges(start, end))
	what <- c("rname", "pos", "qwidth")
	param <- ScanBamParam(which=which, what=what)
	bam <- scanBam(as.character(bamFile), param=param)
	bam <- unname(bam)
	return(do.call("DataFrame", bam))
}

# Convert a list of read in a vector of positions.
#
# INPUT:
# 	currentReads:		The list of read to parse.
#	currentFeature:		The feature to parse.
#
# OUPUT:
#	A vector of with the coverage of every positions calculated from the reads around the max distance from TSS.
convertReadsToDensity <- function(currentReads, currentFeature) {
	maxSize <- abs(currentFeature$end_position - currentFeature$start_position)
	start <- min(currentFeature$start_position, currentFeature$end_position)
	vectorResult <- numeric(maxSize)
	if (nrow(currentReads) > 0) {
		positions <- unlist(mapply(function(x,y) seq(x, x+y), currentReads$pos - currentFeature$start_position, currentReads$qwidth-1))
		# TODO: add unit test -> first position < 0 used to return FALSE and skip current feature (when it was &&)
		positions <- positions[positions > 0 & positions <= maxSize]
		vectorResult <- tabulate(positions, nbins=maxSize)
	}
	return(vectorResult)
}

# Scale the values of a vector to fit with predetermined size
#
# INPUT:
#       values: the values to scale
#       domain: the range to fit the value to
#
# OUTPUT:
#       A vector with the scaled data
scaleVector <- function(values, domain) {
        to_return <- numeric(domain)
        if (length(values) < domain) {
                ratio <- domain / length(values)
                for (i in seq(1, length(values))) {
                        current_start <- round((i - 1) * ratio + 1)
                        current_end <- round(i * ratio)
                        to_return[current_start:current_end] <- values[i]
                }
        }
        else if (length(values) > domain) {
                ratio <- length(values) / domain
                for (i in seq(1, domain)) {
                        current_start <- round((i - 1) * ratio + 1)
                        current_end <- round(i * ratio)
                        to_return[i] <- mean(values[current_start:current_end])
                }

        }
        else {
                return(values)
        }
        return(to_return)
}

# Create a graph
#
# Input:
#	matricesGroups:	A list with every groupsto include in the graph. A group is one or more
#			combination of featureGroup/designGroup. The names must correspond to the
#			names of the matrix object returned from the parse function.
#	data:		The object returned by the parse functions.
#	binSize:	The number of nucleotides in each bin for the bootstrap step.
#	alpha: 		Confidence interval.
# 	sampleSize:	Number of time each bin will be resampled (should be at least 1000).
#	cores:		Number of cores for parallel processing (require parallel package).
# Ouput:
#	The data.frame used to produce the graph
plot.matrices <- function(matricesGroups, data, binSize=100, alpha=0.05, sampleSize=1000, cores=1) {
	# 1. Extract and/or combine relevant matrices
	for (matrixName in names(matricesGroups)) {
		conditions <- unlist(names(data$matrix) %in% matricesGroups[[matrixName]])
		matricesGroups[[matrixName]] <- do.call(rbind, data$matrix[conditions])
	}
	# 2. Bootstrap
	bootstrap <- function(x, bootstrapCores=1) {
		return(bootstrapAnalysis(x, binSize=binSize, alpha=alpha, sampleSize=sampleSize, cores=bootstrapCores))
	}
	bootstrapResults <- applyOnGroups(matricesGroups, cores=1, FUN=bootstrap, bootstrapCores=cores)

	# 3. Prepare data.frame
	DF <- getDataFrame(bootstrapResults)

	# 4. Create graph
	plot.graphic(DF, paste(names(matricesGroups)))
	return(DF)
}

# Perform the bootstrap analysis
#
# Input:
#	currentGroup:	A group extracted from the main data structure
#	binSize:	The number of nucleotides in each bin for the bootstrap step.
#	alpha: 		Confidence interval.
# 	sampleSize:	Number of time each bin will be resampled (should be at least 1000).
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list with the mean of the bootstraped vector and the quartile of order alpha/2 and (1-alpha/2)
#
bootstrapAnalysis.old <- function(currentGroup, binSize, alpha, sampleSize, cores=1) {
	binnedMatrix <- binMatrix(currentGroup$matrix, binSize)
	if (cores > 1) {
		bootResults <- mclapply(1:ncol(binnedMatrix), function(x) binBootstrap.old(binnedMatrix[,x], alpha=alpha, sampleSize=sampleSize), mc.cores=cores)
	} else {
		bootResults <- lapply(1:ncol(binnedMatrix), function(x) binBootstrap.old(binnedMatrix[,x], alpha=alpha, sampleSize=sampleSize))
	}
	bootResults <- do.call(rbind, bootResults)
	toReturn <- list()
	toReturn$mean <- unlist(bootResults[,1])
	toReturn$qinf <- as.numeric(unlist(bootResults[,2]))
	toReturn$qsup <- as.numeric(unlist(bootResults[,3]))
	return(toReturn)
}

# Perform the bootstrap analysis
#
# Input:
#	currentMatrix:	The matrix to use for the bootstrap analysis
#	binSize:	The number of nucleotides in each bin for the bootstrap step.
#	alpha: 		Confidence interval.
# 	sampleSize:	Number of time each bin will be resampled (should be at least 1000).
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list with the mean of the bootstraped vector and the quartile of order alpha/2 and (1-alpha/2)
#
bootstrapAnalysis <- function(currentMatrix, binSize, alpha, sampleSize, cores=1) {
	binnedMatrix <- binMatrix(currentMatrix, binSize)
	if (cores > 1) {
		bootResults <- mclapply(1:ncol(binnedMatrix), function(x) binBootstrap(binnedMatrix[,x], alpha=alpha, sampleSize=sampleSize), mc.cores=cores)
	} else {
		bootResults <- lapply(1:ncol(binnedMatrix), function(x) binBootstrap(binnedMatrix[,x], alpha=alpha, sampleSize=sampleSize))
	}
	bootResults <- do.call(rbind, bootResults)
	toReturn <- list()
	toReturn$mean <- unlist(bootResults[,1])
	toReturn$qinf <- as.numeric(unlist(bootResults[,2]))
	toReturn$qsup <- as.numeric(unlist(bootResults[,3]))
	return(toReturn)
}
# Bin matrix columns
#
# INPUT:
# 	data:		The matrix to bin.
#	binSize:	The number of nucleotides in each bin.
#
# OUPUT:
#	A matrix with each column representing the mean of binSize nucleotides.
binMatrix <- function(data,binSize)
{
	n <- ((ncol(data)-1)/binSize + 1)
	a <- sapply(1:n, function(j){(j-1)*binSize+1})
	newdata <- matrix(0, nrow=nrow(data), ncol=(ncol(data)-1)/binSize)
	for (j in 1:(n-1)) {
		newdata[,j] <- sapply(1:nrow(data), function(i){mean(data[i,a[j]:a[j+1]])}) }
	return(newdata)
}

# Estimate mean and confidence interval of a column using bootstrap.
#
# INPUT:
# 	data:		A vector representing a column from the binned matrix.
#	alpha: 		Confidence interval.
# 	sampleSize:	Number of time each bin will be resampled (should be at least 1000).
#
# OUPUT:
#	A list with the mean of the bootstraped vector and the quartile of order alpha/2 and (1-alpha/2)
binBootstrap.old <- function(data, alpha, sampleSize)
{
	# TODO: it would probably be more efficient to parallelize here
	size <- length(data)
	# TODO: try to sample data directly
	S <- matrix(replicate(sampleSize, data[sample(1:length(data),size,replace=TRUE)]), nrow=sampleSize)
	# TODO: We should also be able to plot non-bootstrapped mean
	# TODO: Try to parallelize at this level instead to lower memory consumption
	mean <- mean(sapply(1:sampleSize, function(i){mean(S[i,])}))
	qinf <- quantile(sapply(1:sampleSize, function(i){mean(S[i,])}), prob=alpha/2)
	qsup <- quantile(sapply(1:sampleSize, function(i){mean(S[i,])}), prob=(1-alpha/2))
	liste <- list(mean=mean, qinf=qinf, qsup=qsup)
	return(liste)
}

binBootstrap <- function(data, alpha, sampleSize, cores=cores)
{
	# TODO: it would probably be more efficient to parallelize here
	size <- length(data)
	# TODO: try to sample data directly
	S <- matrix(replicate(sampleSize, data[sample(1:length(data),size,replace=TRUE)]), nrow=sampleSize)
	# TODO: We should also be able to plot non-bootstrapped mean
	mean <- mean(sapply(1:sampleSize, function(i){mean(S[i,])}))
	qinf <- quantile(sapply(1:sampleSize, function(i){mean(S[i,])}), prob=alpha/2)
	qsup <- quantile(sapply(1:sampleSize, function(i){mean(S[i,])}), prob=(1-alpha/2))
	liste <- list(mean=mean, qinf=qinf, qsup=qsup)
	return(liste)
}

# Convert the main data structure into a data.frame
#
# Input:
#	groups:	The main data structure
#
# Output:
#	A data.frame with the condensed results from the main data structure
#		Columns:
#		 * Groups: name of current group
#		 * distances: the number of bin for each entry
#		 * means: the means to plot
#		 * qinf: the lower end of the confidence interval
#		 * qsup: the higher end of the confidence interval
plot.getDataFrame.old <- function(groups) {
	lists <- lapply(groups, function(x) return(x$bootstrap))
	grid = seq(-5000,5000,length=length(lists[[1]]$mean))
	DF= data.frame (
		Groups <- factor(rep(names(groups), each=length(grid))),
		distances <- rep(grid, length(lists)),
		means <- c(sapply(1:length(lists), function(x) lists[[x]]$mean)),
		qinf <-  c(sapply(1:length(lists), function(x) lists[[x]]$qinf)),
		qsup <-  c(sapply(1:length(lists), function(x) lists[[x]]$qsup))
	)
	colnames(DF) <- c("Groups", "distances", "means", "qinf", "qsup")
	return(DF)
}

# Convert the bootstrapped data into a data.frame
#
# Input:
#	bootstapData:	Data produced during the bootstrap analysis
#
# Output:
#	A data.frame with the condensed results from the main data structure
#		Columns:
#		 * Groups: name of current group
#		 * distances: the number of bin for each entry
#		 * means: the means to plot
#		 * qinf: the lower end of the confidence interval
#		 * qsup: the higher end of the confidence interval
getDataFrame <- function(bootstrapData) {
	#lists <- lapply(groups, function(x) return(x$bootstrap))
	#grid = seq(-5000,5000,length=length(bootstrapData[[1]]$mean))
	grid = 1:length(bootstrapData[[1]]$mean)
	DF = data.frame (
		Groups <- factor(rep(names(bootstrapData), each=length(grid))),
		distances <- rep(grid, length(bootstrapData)),
		means <- c(sapply(1:length(bootstrapData), function(x) bootstrapData[[x]]$mean)),
		qinf <-  c(sapply(1:length(bootstrapData), function(x) bootstrapData[[x]]$qinf)),
		qsup <-  c(sapply(1:length(bootstrapData), function(x) bootstrapData[[x]]$qsup))
	)
	colnames(DF) <- c("Groups", "distances", "means", "qinf", "qsup")
	return(DF)
}

# Produce a plot with based on a data.frame
#
# Input:
#	DF:	The data frame produced by the plot.getDataFrame function
#	title:	The title of the graph
#
# Ouput:
#	The graph that is printed on the screen.
plot.graphic <- function(DF, title) {
	library(ggplot2)
	p <- ggplot(DF, aes(x=distances, y=means, ymin=qinf, ymax=qsup)) +
	geom_ribbon(aes(fill=Groups), alpha=0.3) +
	geom_line(aes(color=Groups),size=1,litype=1,bg="transparent")+
	theme(panel.grid.major = element_line())+
	theme(panel.grid.minor = element_line())+
	theme(panel.background = element_blank())+
	theme(panel.background = element_rect())+
	theme_bw()+
	theme(axis.title.x = element_blank())+
	ylab("Mean of RPM for each 100 bps")+
	ggtitle(title)
	print(p)
}
