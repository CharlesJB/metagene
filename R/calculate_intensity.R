# Created by Charles Joly Beauparlant
# 2013-11-26

# Create a metagene plot based on a list of genomic features
#
# Input:
#	bamFiles:	A vector of bamFile to plot. TODO: Should also accept a list of bamfiles where each elements would be grouped together
#	features:	Either a filename of a vector of filenames.
#			Supported features: ensembl_gene_id
#			File must contain a header that correspond to the name of the group
#	specie:		hs: Homo sapiens (default) / mm: Mus musculus
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
plotFeatures <- function(bamFiles, features=NULL, specie="hs", maxDistance=5000, design=NULL, binSize=100, alpha=0.05, sampleSize=1000, cores=1) {
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
	groups <- prepareGroups(names(featuresGroups), bamFiles=bamFiles, design=design)
	# TODO: When groups of bam files are implemented in design, we need to change the following function
	#	that will return an element in the list for each combination of group of gene and group of bam
	#	I.E.: 2 bam groups and 2 feature groups will mean 4 groups in total
	# TODO: Make a function to prepare and parse groups
	parseGroup <- function(currentGroup) {
		# Extract the data.frame corresponding the current group in the list of groups
		print(paste("Current group:", currentGroup))
		currentFeatures <- featuresGroups[[which(names(featuresGroups) == currentGroup)]]
		listMatrix <- parseBamFiles(bamFiles, currentFeatures, cores=cores)
		if (!is.null(design)) {
			listMatrix <- mergeDesign(listMatrix, design)
		}
		return(listMatrix)
	}
	listMatrixByGroup <- lapply(names(featuresGroups), parseGroup)
	cat("Step 3: Parse bam files... Done!\n")

	# 4. Merge matrix
	cat("Step 4: Merge matrix...")
	mergeMatrix <- function(bamList) {
		lapply(1:length(bamList), function(x) bamList[[x]] <- do.call(rbind,bamList[[x]]))
	}
	mergedMatrix <- lapply(1:length(listMatrixByGroup), function(x) do.call(rbind, mergeMatrix(listMatrixByGroup[[x]])))
	cat(" Done!\n")

	# 5. Bootstrap
	cat("Step 5: Bootstrap...")
	bootstrapData <- function(i) {
		binnedMatrix <- binMatrix(mergedMatrix[[i]],binSize)
		if (cores > 1) {
			bootstrapedData <- mclapply(1:ncol(binnedMatrix), function(x) binBootstrap(binnedMatrix[,x], alpha=alpha, sampleSize=sampleSize), mc.cores=cores)
		} else {
			bootstrapedData <- lapply(1:ncol(binnedMatrix), function(x) binBootstrap(binnedMatrix[,x], alpha=alpha, sampleSize=sampleSize))
		}
		return(bootstrapedData)
	}
	bootstrapedDataList <- lapply(1:length(mergedMatrix), bootstrapData)
	cat(" Done!\n")
	return(bootstrapedDataList)
	# TODO: Check param with Rawane
	# TODO: Add previous params to plotFeatures function
	# 6. Plot
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

# Check parameters for the plot functions
checkParams <- function(bamfiles, features=NULL, maxDistance=NULL, ranges=NULL, design=NULL, scaling="median", filling=NULL, padding=NULL, centering=NULL) {

}

# Sort and index bam files, if necessary. Return the number of aligned reads for each bam file.
#
# Input:
#	bamFiles:	Vector containing the list of every bam filename to be included in the analysis.
#
# Output:
#	A data.frame containing the indexed bam filename and number of aligned reads for each bam file.
#	Column names: bamFiles and alignedCount
prepareBamFiles <- function(bamFiles, cores = 1) {
	library(Rsamtools)

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
#	specie: "mouse" or "human"
#
# OUTPUT:
#	A data.frame with 5 columns:
#		1: feature -> Ensembl gene id
#		2: strand -> -1 or 1
#		3: space -> chromosome
#		4: start_position -> position of the TSS
#		5: end_position -> ending position of the last exon of the gene
getGenes <- function(specie) {
	require(biomaRt)
	# Set the correct specie
	if (specie == "hs") {
		chrom <- c(as.character(seq(1,21)),"X","Y")
		ensmart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
	} else if (specie == "mm") {
		chrom <- c(as.character(seq(1,19)),"X","Y")
		ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
	} else {
		print("Incorrect parameter for specie name")
		print("Currently supported specie are human and mouse")
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
#	specie:		hs: Homo sapiens (default) / mm: Mus musculus
#	maxDistance:	The distance around feature to include in the plot.
#	cores:		Number of cores for parallel processing (require parallel package).
#
# Output:
#	A list of data.frame. One data.frame by group of features.
#	The names of each element of the list correspond to the name of the group.
prepareRegions <- function(features, specie="hs", maxDistance=5000, cores=1) {
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
	if (cores > 1) {
		library(parallel)
		featuresGroups <- mclapply(features, extractFeatures, mc.cores=cores)
	} else {
		featuresGroups <- lapply(features, extractFeatures)
	}
	names(featuresGroups) <- unlist(lapply(features, function(x) as.character(read.table(x, nrow=1)[1,])))
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
#	A list of list.
#	First level is the name of the group.
#	Second level is the correct bam file for each group.
# 	Note1: If there are groups of feature and groups of bam file, every possible combination of feature groups / bam groups will be created.
# 	Note2: The names of the elements the bam list are the same as the ones in the design file, which are also the original bam file names.
#	       The content of the each element of the list of bam file is the sorted bam file name, which should be use when parsing bam.
prepareGroups <- function(featuresGroupsNames, bamFiles, design=NULL) {
	allGroups <- list()
	for (featureName in featuresGroupsNames) {
		# If there is no design file, we put every bamFile in each groups
		if (is.null(design)) {
			allGroups[[featureName]] <- as.list(bamFiles$bam)
			names(allGroups[[featureName]]) <- bamFiles$oldBam
		} else {
			for (i in 2:ncol(design)) {
				name <- paste(featureName, colnames(design)[i], sep="_")
				bam <- design[design[,i] != 0,]$Sample
				# We want the sorted name of the bam files as the content of the list
				j <- numeric()
				for (bamFile in bam) {
					j <- append(j, which(bamFiles$oldBam == bamFile))
				}
				sortedBamNames <- as.character(bamFiles[j,]$bam)
				allGroups[[name]] <- as.list(sortedBamNames)
				names(allGroups[[name]]) <- bam
			}
		}
	}
	return(allGroups)
}


#parseBamFiles <- function(bamFiles, features, cores=1) {
parseBamFiles <- function(features, groups, cores=1) {
	extractReadsDensity <- function(feature, bamFile) {
		# Extract raw counts
		currentReads <- extractReadsInRegion(bamFile, feature$space, feature$start_position, feature$end_position)
		vectorResult <- convertReadsToDensity(currentReads, feature)

		# If on negative strand, invert the current vector
		if (feature$strand == "-1" |  feature$strand == -1 | feature$strand == "-") {
			vectorResult <- rev(vectorResult)
		}

		# Convert to RPM
		currentAlignedCount <- bamFiles[bamFiles$bam == bamFile,]$alignedCount
		vectorResult <- vectorResult / (currentAlignedCount / 1000000)
		return(vectorResult)

	}
	parseBam <- function(bamFile, features) {
		print(paste("Current bam:", bamFile))
		if (cores > 1) {
			library(parallel)
			return(mclapply(1:nrow(features), function(x) extractReadsDensity(features[x,],  bamFile=as.character(bamFile)), mc.cores=cores))
		} else {
			return(lapply(1:nrow(features), function(x) extractReadsDensity(features[x,],  bamFile=as.character(bamFile))))
		}
	}

	return(lapply(bamFiles$bam, parseBam, features=features))
}

# Extract read density from a region
#
# Input:
#	geneID:		The current ensembl gene id to parse
#	knownGenes:	The annotation to translate ID into genomic positions
#	bamFiles:	The data.frame obtained with the prepareBamFiles function.
#	maxDistance:	The distance on each side of the beginning of feature to include in the analysis.
#
# Output:
# 	A matrix with the read density for the current regions with as many line as there are bam files.
# TODO: return count as read per million aligned read (RPM)
# TODO: replace knownGenes by GRanges
getRegionReadDensity <- function(geneID, knownGenes, bamFiles, maxDistance) {
	extractReadsDensity <- function(bamfile) {
		# Fetch infos from current feature
		currentFeature <- knownGenes[knownGenes$feature == geneID,]
		currentSpace <- currentFeature$space
		currentStart <- currentFeature$start - maxDistance
		currentEnd <- currentFeature$start + maxDistance
		currentStrand <- currentFeature$strand

		# Extract raw counts
		currentReads <- extractReadsInRegion(bamFile, currentSpace, currentStart, currentEnd)
		vectorResult <- convertReadsToDensity(currentReads, currentFeature$start, maxDistance)

		# If on negative strand, invert the current vector
		if (currentStrand == "-1" |  currentStrand == -1 | currentStrand == "-") {
			vectorResult <- rev(vectorResult)
		}

		# Convert to RPM
		currentAlignedCount <- bamFiles[bamFiles$bam == bamFile,]$alignedCount
		vectorResult <- vectorResult / (currentAlignedCount / 1000000)
		return(vectorResult)
	}

	# TODO: use parallel with next line (?)
	listResults <- lapply(bamFiles$bam, extractReadsDensity)
	#rawMatrix <- lapply(nrow(featuresGroups), extractReadsDensity)
	#rawMatrix <- do.call(rbind, rawMatrix)
	#return(rawMatrix)
}

# Extract read density information from a list of regions and convert them into a matrix
#
# Input:
#	bamFile:	The bam file to parse.
#	regions:	GRanges object representing the list of genomic ranges for the current group of file.
#
# Output:
#	A matrix with as much columns as the largest region and as many line as the number of regions.
parseRegions <- function(bamFile, regions) {

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
		positions <- positions[positions > 0 && positions <= maxSize]
		vectorResult <- tabulate(positions, nbins=maxSize)
	}
	return(vectorResult)
}

# Bin matrix columns
#
# INPUT:
# 	data:		The matrix to bin.
#	binSize:	The number of nucleotides in each bin.
#
# OUPUT:
#	A matrix with each column representing the mean of binSize nucleotides.
binMatrix <- function(data,bin)
{
	n <- ((ncol(data)-1)/bin + 1)
	a <- sapply(1:n, function(j){(j-1)*bin+1})
	newdata <- matrix(0, nrow=nrow(data), ncol=(ncol(data)-1)/bin)
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
	S <- matrix(replicate(sampleSize, data[sample(1:length(data),size,replace=TRUE)]), nrow=sampleSize)
	mean <- mean(sapply(1:sampleSize, function(i){mean(S[i,])}))
	qinf <- quantile(sapply(1:sampleSize, function(i){mean(S[i,])}), prob=alpha/2)
	qsup <- quantile(sapply(1:sampleSize, function(i){mean(S[i,])}), prob=(1-alpha/2))
	liste <- list(mean=mean, qinf=qinf, qsup=qsup)
	return(liste)
}
