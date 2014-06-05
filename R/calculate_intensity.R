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
	featuresGroups <- prepareRegions(features=features, specie=specie, cores=cores)
	cat(" Done!\n")

	# 3. Parse bam files
	cat("Step 3: Parse bam files...\n")
	groupNames <- prepareGroups(names(featuresGroups), bamFiles=bamFiles, design=design)
	groups <- parseBamFiles(bamFiles, featuresGroups, groups=groupNames, design=design, cores=cores)
	cat("Step 3: Parse bam files... Done!\n")

	# 4. Merge matrix
	cat("Step 4: Merge matrix...")
	# Convert list of vectors in matrix for a single group
	mergeMatrixFixedLengthFeatures <- function(currentExperiment) {
		currentBamFiles <- currentExperiment$bamFiles
		return(lapply(1:length(currentBamFiles), function(x) currentBamFiles[[x]] <- do.call(rbind, currentBamFiles[[x]])))
	}
	mergeMatrix <- function(x) {
		x$matrix <- do.call(rbind, mergeMatrixFixedLengthFeatures(x))
		return(x)
	}
	groups <- applyOnGroups(groups=groups, cores=cores, FUN=mergeMatrix)
	cat(" Done!\n")

	# 5. Bootstrap
	cat("Step 5: Bootstrap...")
	bootstrap <- function(x) {
		x$bootstrap <- bootstrapAnalysis(x, binSize=binSize, alpha=alpha, sampleSize=sampleSize, cores=cores)
		return(x)
	}
	groups <- applyOnGroups(groups=groups, cores=1, FUN=bootstrap)
	cat(" Done!\n")

	# TODO: Add previous params to plotFeatures function
	# 6. Plot
	cat("Step 6: Plot...")
	groups$graphData <- plot.getDataFrame(groups)
	plot.graphic(groups$graphData, paste(names(groups)))
	cat(" Done!\n")
	return(groups)
}

# Create a metagene plot based on a list of genomic regions
#
# Input:
#	bamFiles:	A vector of bamFile to plot. TODO: Should also accept a list of bamfiles where each elements would be grouped together
#	file:		The name of the output file. In pdf format. TODO: also send to screen
#	ranges:		List of genomic regions. Path to a bed file. TODO: add genomic ranges class.
#			TODO: Hard coder FANTOM enhancer -> default
# 	design:		A matrix explaining the relationship between multiple samples.
#			One line per samples.
#			One column per group of samples. For example, biological replicates and corresponding controls are in the same group.
#			1: treatment file(s)
#			2: control file(s)
#	scaling:	The way the regions should be scaled. Cannot be used with filling.
#			"median": Default. Every regions are resized to fit the median size of the regions.
#			"average": Every regions are resized to fit the average size of the regions.
#			"largest": Every regions are resized to fit the largest size of the regions.
#			"smallest": Every regions are resized to fit the smallest size of the regions.
#			NULL: No scaling.
#	filling:	Define how are smaller regions are fitted to the largest. Cannot be used with scaling.
#			"extend": Resize all regions to fit largest.
#			"NA": Fill smaller regions with NA.
#			NULL: Default. No filling.
#	padding:	Add nucleotides to the size of every regions.
#			TODO: 2 arguments, type de padding (absolute ou relative) et taille du padding (paddingSize)
#			If the value is greater than one, the rounded value will be added in both directions.
#			If the value is smaller than one, a relative padding corresponding to the value will be added in both directions.
#			NULL: Default. No padding
#	centering:	A function to define how are the regions centered.
#			NULL: Default. Regions are centered at the middle.
#	...:		Extra arguments for centering function.
#
# TODO: nouvel argument pour normalisation qui fit avec edgeR
# TODO: nouvel argument pour overlapping regions (merge, mean: on fait la moyenne)
#
# output:
#
plotFeaturesByRegions <- function(bamFiles, file="regions.pdf", ranges=NULL, design=NULL, scaling="median", filling=NULL, padding=NULL, centering=NULL, ...) {
	# 0. Check if params are valid
	# 1. Prepare bam files
	# 2. Parse regions
	# 	For loop for now. Eventually, this is the place that will be parallelized.
	# 3. Realign regions
}

# Extract the names from the main data structure
#
# Input:
#	groups:	The main data structure
#
# Output:
#	A data structure similar to the main data structure of MetaFeatures, but with no data
extractNames <- function(groups) {
	toReturn <- list()
	for (groupName in names(groups)) {
		toReturn[[groupName]] <- list()
		toReturn[[groupName]]$featureName <- groups[[groupName]]$featureName
		toReturn[[groupName]]$designName <- groups[[groupName]]$designName
		toReturn[[groupName]]$bamFiles <- list()
		for (bamFile in groups[[groupName]]$bamFiles) {
			toReturn[[groupName]]$bamFiles[[bamFile]] <- ""
		}
	}
	return(toReturn)
}

# Copy extracted names to the main data structures
#
# Input:
#	newNames:	The output of extractNames function
#	groups:		The main data structure of MetaFeatures
#
# Output:
#	The group input with new names
copyNames <- function(newNames, groups) {
	names(groups) <- names(newNames)
	for (groupName in names(newNames)) {
		for (i in 1:length(groups[[groupName]])) {
			if (class(groups[[groupName]][[i]]) == "list") {
				names(groups[[groupName]])[i] <- "bamFiles"
			}
		}
		groups[[groupName]]$featureName <- newNames[[groupName]]$featureName
		groups[[groupName]]$designName <- newNames[[groupName]]$designName
		names(groups[[groupName]]$bamFiles) <- names(newNames[[groupName]]$bamFiles)
	}
	return(groups)
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

# Apply a function on every bam files of the main data structure
#
# Input:
#	groups:		The main data structure
#	cores:		Number of cores for parallel processing (require parallel package).
#	FUN:		The function to apply on every groups
#	...:		Extract arguments for FUN
#
# Output:
#	The result of the function applied
# TODO: This function should be made more general and be debugged
#applyOnBamFiles <- function(groups, cores=1, FUN, ...) {
	#oldNames <- extractNames(groups)
	#for (group in groups) {
		#for (experiment in group$bamFiles) {
			#if (cores > 1) {
				#library(parallel)
				#groups[[group]][[experiment]] <- mclapply(groups[[group]][[experiment]], function(x) FUN(x, ...), mc.cores=cores)
			#} else {
				#groups[[group]] <- lapply(groups[[group]][[experiment]], function(x) FUN(x, ...))
			#}
		#}
	#}
	#return(copyNames(oldNames, groups))
#}


# Check parameters for the plot functions
checkParams <- function(bamfiles, features=NULL, maxDistance=NULL, ranges=NULL, design=NULL, scaling="median", filling=NULL, padding=NULL, centering=NULL) {

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
		stop("At least one BAM file name is not a valid name (not a character string).")
	}
	
	# All BAM files must exist
 	if (sum(unlist(lapply(bamFiles, file.exists))) != length(bamFiles)) {
 		stop("At least one BAM file does not exist.")
 	}
	
	# This function will only index a file if there is no index file
	indexBamFiles <- function(bamFile) {
		if (file.exists(paste(bamFile, ".bai", sep=""))  == FALSE) {
			# If there is no index file, we sort and index the current bam file
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
# Output:
#	A list of data.frame. One data.frame by group of features.
#	The names of each element of the list correspond to the name of the group.
prepareRegions <- function(features, specie="human", maxDistance=5000, cores=1) {
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
parseBamFiles <- function(bamFiles, featuresGroups, groups, design=NULL, cores=1) {
	parseGroup <- function(groupName, featuresGroups) {
		print(paste("parseGroup:", groupName))
		# Get bam files
		currentBamFiles <- groups[[groupName]]$bamFiles

		# Get features
		currentFeatureName <- groups[[groupName]]$featureName
		currentFeatures <- featuresGroups[[currentFeatureName]]

		# Parse each bam
		currentBamFiles <- lapply(currentBamFiles, function(x) parseBamFile(x, bamFiles[bamFiles$bam == x,]$alignedCount, currentFeatures, cores))
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
		groups <- removeControlsFromGroups(groups=groups, design=design, cores=cores)
	}
	return(groups)
}

# Substract controls from multiple groups
#
# Input:
#	groups:	The main data structure
#	design:	A matrix explaining the relationship between multiple samples.
#	cores:	Number of cores for parallel processing (require parallel package).
#
# Output:
#	The main data structure from which the controls were substracted then deleted
removeControlsFromGroups <- function(groups, design, cores=1) {
	newGroups <- lapply(groups, function(x) removeControls(currentGroup=x, design=design, cores=cores))
	names(newGroups) <- names(groups)
	return(newGroups)
}

# Substract controls from a single group
#
# Input:
#	currentGroup:	A group extracted from the main data structure
#	design:		A matrix explaining the relationship between multiple samples.
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	The group extracted from the main data structure from which the controls were substracted then deleted
removeControls <- function(currentGroup, design, cores=1) {
	# 1. Extract relevant design columns
	currentDesign <- design[,c(1,which(colnames(design) == currentGroup$designName))]
	treatmentNames <- as.character(currentDesign[currentDesign[,2] == 1, 1])
	controlNames <- as.character(currentDesign[currentDesign[,2] == 2, 1])

	# 2. Merge controls
	# TODO
	if (length(controlNames) > 1) {
		print("TODO")
	} else {
		mergedControls <- currentGroup$bamFiles[[controlNames[1]]]
		currentGroup$bamFiles[[controlNames[1]]] <- NULL
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
	currentGroup$bamFiles <- lapply(currentGroup$bamFiles, substractControl)
	names(currentGroup$bamFiles) <- treatmentNames
	return(currentGroup)
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
parseBamFile <- function(bamFile, alignedCount, features, cores=1) {
	print(paste("Current bam:", bamFile))
	extractReadsDensity <- function(feature, bamFile) {
		# Extract raw counts
		currentReads <- extractReadsInRegion(bamFile, feature$space, feature$start_position, feature$end_position)
		vectorResult <- convertReadsToDensity(currentReads, feature)

		# If on negative strand, invert the current vector
		if (feature$strand == "-1" |  feature$strand == -1 | feature$strand == "-") {
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

# Extract reads from bam file that overlap with a genomic region
#
# INPUT:
#	bamFile:	Path to the bam file.
#	chr:		Current chromosome.
#	start:		Starting position of the current region.
#	end:		Ending position of the current region.
#
# OUTPUT:
#	A data.frame containing every reads overlapping the current genomic region:
#		* rname
#		* pos
#		* qwidth
extractReadsInRegion <- function(bamFile, chr, start, end) {
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
bootstrapAnalysis <- function(currentGroup, binSize, alpha, sampleSize, cores=1) {
	binnedMatrix <- binMatrix(currentGroup$matrix, binSize)
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
binBootstrap <- function(data, alpha, sampleSize)
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
plot.getDataFrame <- function(groups) {
	lists <- lapply(groups, function(x) return(x$bootstrap))
	grid = seq(-5000,5000,length=length(lists[[1]]$mean))
	DF = data.frame (
		Groups <- factor(rep(names(groups), each=length(grid))), 
		distances <- rep(grid, length(lists)),
		means <- c(sapply(1:length(lists), function(x) lists[[x]]$mean)),
		qinf <-  c(sapply(1:length(lists), function(x) lists[[x]]$qinf)),
		qsup <-  c(sapply(1:length(lists), function(x) lists[[x]]$qsup))
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
