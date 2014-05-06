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
#	binSize:	The number of nucleotides in each bin.
plotFeatures <- function(bamFiles, features=NULL, specie="hs", maxDistance=5000, design=NULL, binSize=100) {
	# 0. Check if params are valid
	# 1. Prepare bam files
	# 2. Prepare regions
	# TODO: We need to iterate on experiment in design files and produce a matrix for each experiment (output: a list of matrices?)
	#	The loop should start here and stop before Bootstrap
	knownGenes <- getGenes(specie)
	allFeatures <- data.frame()
	for (filename in features) {
		currentFeatures <- read.table(filename, header=TRUE)
		currentGroup <- colnames(currentFeatures)[1]
		currentFeatures$group <- currentGroup
		colnames(currentFeatures) <- c("feature", "group")
		allFeatures <- rbind(allFeatures, currentFeatures)
	}
	# 3. Parse regions
	rawMatrix <- lapply(nrow(allFeatures), getRegionReadDensity(allFeatures[x,]$feature), knownGenes=knowGenes, bamFiles=bamFiles)
	rawMatrix <- do.call(rbind, rawMatrix)
	normalizedMatrix <- rawMatrix # TODO: change this when the design analysis is implemented
	# 4. Bootstrap
	#bootstrapedMatrix <- binBootstrap(normalizedMatrix, binSize, alpha=0.05, nech=1000, size=???)
	# TODO: Check param with Rawane
	# TODO: Add previous params to plotFeatures function
	# 5. Plot
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
prepareBamFiles <- function(bamFiles) {
	library(Rsamtools)
	results <- data.frame()

	# This function will only index a file if there is no index file
	indexBamFiles <- function(bamFile) {
		if (file.exists(paste(bamFile, ".bai", sep=""))  == FALSE) {
			# If there is no index file, we sort and index the current bam file
			sortedBamFile <- sub("^([^.]*).*", "\\1", bamFile)
			sortedBamFile <- paste(sortedBamFile, ".sorted", sep="")
			sortBam(bamFile, sortedBamFile)
			sortedBamFile <- paste(sortedBamFile, ".bam", sep="")
			indexBam(sortedBamFile)
			bamFile <- sortedBamFile
		}
		return(bamFile)
	}
	results$bamFiles <- sapply(bamFiles, indexBamFiles)

	# This function will calculate the number of aligned read in each bam file
	countAlignedReads <- function(bamFile) {
		return(countBam(bamFile, param=ScanBamParam(flag = scanBamFlag(isUnmappedQuery=FALSE)))$records)
	}
	results$alignedCount <- sapply(bamFiles, countAlignedReads)

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

# Extract read density from a region
#
# Input:
#	ensembl_gene_id:	The current ensembl gene id to parse
#	knownGenes:		The annotation to translate ID into genomic positions
#	bamFiles:		A vector of bam file to parse.
#
# Output:
# 	A matrix with the read density for the current regions with as many line as there are bam files.
# TODO: return count as read per million aligned read (RPM)
getRegionReadDensity <- function(geneID, knowGenes, bamFiles, maxDistance) {
	extractReadsDensity <- function(bamfile) {
		currentFeature <- knowGenes[knowGenes$feature == geneID,]
		currentReads <- extractsReadsInRegion(bamFile, currentFeature$space, currentFeature$start, currentFeature$end)
		vectorResult <- convertReadsToPosVector(currentReads, currentFeature, maxDistance)
		# If on negative strand, invert the current vector
		currentStrand <- currentFeature$strand
		if (currentStrand == "-1" |  currentStrand == -1 | currentStrand == "-") {
			vectorResult <- rev(vectorResult)
		}
		return(vectorResult)
	}

	rawMatrix <- lapply(nrow(allFeatures), extractReadsDensity)
	rawMatrix <- do.call(rbind, rawMatrix)
	return(rawMatrix)
}

# Convert a list of feature into genomic regions
#
# Input:
#	features:	A vector of Ensembl gene IDs.
#	knownGenes:	A data.frame with known genes
#			Must contain the following columns:
#			 feature -> Ensembl gene id
#			 strand -> -1 or 1
#			 space -> chromosome
#			 start_position -> position of the TSS
#			 end_position -> ending position of the last exon of the gene
#	maxDistance:	The distance around feature to include in the plot.
#
# Output:
#	A data.frame with the genomic positions (4 columns)
#		1: feature -> Ensembl gene id
#		2: strand -> -1 or 1
#		3: space -> chromosome
#		4: start -> position of the region
#		5: end -> ending position of the region
prepareRegions <- function(features, knownGenes, maxDistance) {
	# Keep only ensembl id in features
	if (is.null(features)) {
		result <- knownGenes
	} else {
		positions <- data.frame()
		for (feature in features$feature) {
			currentIndex <- which() # TODO: remove prepareRegion function. Replace with a apply function that takes a single feature id and return a single line of the matrix
		}
	}
	# Expand the regions
	result[result$start < result$end]$start <- result[result$start > result$end]$start - maxDistance
	result[result$start < result$end]$end <- result[result$start > result$end]$end + maxDistance
	result[result$start > result$end]$start <- result[result$start > result$end]$start + maxDistance
	result[result$start > result$end]$end <- result[result$start > result$end]$end - maxDistance
	return(result)
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
#	bam_file:	Path to the bam file.
#	chr:		Current chromosome.
#	start:		Starting position of the current region.
#	end:		Ending position of the current region.
#
# OUTPUT:
#	A data.frame containing every reads overlapping the current genomic region:
#		* rname
#		* pos
#		* qwidth
extractsReadsInRegion <- function(bam_file, chr, start, end) {
	suppressMessages(library(Rsamtools))
	df <- data.frame()
	#which <- GRanges(seqnames=Rle(chr), ranges=IRanges(start, end),strand=Rle(strand))
	which <- GRanges(seqnames=Rle(chr), ranges=IRanges(start, end))
	what <- c("rname", "pos", "qwidth")
	param <- ScanBamParam(which=which, what=what)
	bam <- scanBam(bam_file, param=param)
	.unlist <- function(x) {
		## do.call(c, ...) coerces factor to integer, which is undesired
		x1 <- x[[1L]]
		if (is.factor(x1)) {
			structure(unlist(x), class = "factor", levels = levels(x1))
		} else {
			do.call(c, x)
		}
	}
	bam <- unname(bam)
	elts <- setNames(bamWhat(param), bamWhat(param))
	lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))
	return(do.call("DataFrame", lst))
}

# Convert a list of read in a vector of positions.
#
# INPUT:
# 	current_reads:		The list of read to parse.
#	current_TSS_offset:	The offset of the reads from the TSS.
# 	max_distance:		Maximal distance from TSS.
#
# OUPUT:
#	A vector of with the coverage of every positions calculated from the reads around the max distance from TSS.
convertReadsToPosVector <- function(current_reads, current_TSS_offset, max_distance) {
	vector_result <- numeric(max_distance*2+1)
	if (nrow(current_reads) > 0) {
		positions <- unlist(mapply(function(x,y) seq(x, x+y), current_reads$pos - current_TSS_offset, current_reads$qwidth-1))
		positions <- positions[abs(positions)<=max_distance] # to remove reads beyond max distance
		positions <- positions + max_distance
		vector_result <- tabulate(positions, nbins=max_distance*2+1)
		#if (length(positions) > 0) { # TODO: Change for tabulate?
			#for (i in positions) {
				#vector_result[i] <- vector_result[i] + 1
			#}
		#}
	}
	return(vector_result)
}

###########################################################################################################################################
############## Fonctions utilisées pour estimer les moyennes et les bandes de confiances ##################################################
###########################################################################################################################################


### La fonction binBootstrap ci-après estime par boostrap les paramètres (moyenne, quartiles) d'une colonne d'une matrice issue d'un binage
### Elle prend six agrguments:
### data: une matrice de données,
### bin: le nombre de colonnes qui caracetérisent le binage,
### column: la colonne de la matrice binée sur laquelle s'opère le boostrap,
### alpha: (par défaut fixé à 0.05) définit le seuil des bandes de confiance,
### nech: le nombre d'échantillons boostrap (au minimum 1000),
### size: la taille des échantillons boostrap (par défaut size est égal à la taille de l'échantillon initial).
### La fonction retourne un élément de type liste dont les éléments sont la moyenne et les quartiles d'ordre alpha/2 et (1-alpha/2).
binBootstrap <- function(data,bin,column,alpha,nech,size)
{
	### La fonction newData ci-dessous a pour but de regrouper et de moyenner des colonnes d'une matrice de données.
	### Elle prend deux arguments:
	### data: une matrice de données normalisées
	### bin: le nombre de colonnes qui caracetérisent un groupe.
	### Cette fonction retourne une nouvelle matrice dont chacune des colonnes correspond à une moyenne sur bin colonnes de l'ancienne matrice
	newData <- function(data,bin)
	{
		n <- ((ncol(data)-1)/bin + 1)
		a <- sapply(1:n, function(j){(j-1)*bin+1})
		newdata <- matrix(0, nrow=nrow(data), ncol=(ncol(data)-1)/bin)
		for (j in 1:(n-1)) {
			newdata[,j] <- sapply(1:nrow(data), function(i){mean(data[i,a[j]:a[j+1]])}) }
		return(newdata)
	}
	data.binage <- newData(data,bin)
	X <- data.binage[,column]
	S <- matrix(replicate(nech, X[sample(1:length(X),size,replace=TRUE)]), nrow=nech)
	mean <- mean(sapply(1:nech, function(i){mean(S[i,])}))
	qinf <- quantile(sapply(1:nech, function(i){mean(S[i,])}), prob=alpha/2)
	qsup <- quantile(sapply(1:nech, function(i){mean(S[i,])}), prob=(1-alpha/2))
	liste <- list(mean=mean, qinf=qinf, qsup=qsup)
	return(liste)
}

###########################################################################################################################################
############## Fonctions utilisées pour tester l'égalité des moyennes d'enrichissement de deux groupes ##################################################
###########################################################################################################################################

### La fonction pvalue calcule la p-valeur boostrap du test bilatéral sur deux groupes.
### Elle prend quatre arguments:
### X1: un échantillon issu du groupe 1,
### X2: un échantillon du groupe 2,
### nech: le nombre d'échantillons boostrap de chaque groupe,
### size: la taille des échantillons boostrap.
pvalue <- function(X1, X2, nech, size) {
	### La fonction statistic ci-dessous est utilisée pour tester l'égalité entre les moyennes de deux groupes.
	### Elle prend deux arguments:
	### X1: un échantillon issu du groupe 1,
	### X2: un échantillon du groupe 2.
	### Elle retourne la valeur de la statistique du test bilatéral sur les deux groupes.
	statistic <- function(X1,X2)
	{
		var.pool = ((length(X1)-1)*var(X1)+(length(X2)-1)*var(X2))/(length(X1)+length(X2)-2)
		stat =   (mean(X1)-mean(X2))/sqrt(var.pool(X1,X2)*(1/length(X1)+1/length(X2)))
		return(stat)
	}
	data1 <- matrix(replicate(nech, X1[sample(1:length(X1),size,replace=TRUE)]), nrow=nech)
	data2 <- matrix(replicate(nech, X2[sample(1:length(X2),size,replace=TRUE)]), nrow=nech)
	stat.boot <- sapply(1:nech,function(i){statistic(data1[i,],data2[i,])})
	stat.obs <- abs(statistic(X1,X2))
	pval <- length(stat.boot[stat.boot >= stat.obs])/nech  + length(stat.boot[stat.boot <= -stat.obs])/nech
	return(pval)
}
