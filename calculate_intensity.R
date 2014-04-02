# Created by Charles Joly Beauparlant
# 2013-11-26

# Extract reads from bam file that overlap with enriched peaks.
#
# INPUT:
#	bam_file: 		Path to the bam file.
# 	annotated_peaks:	Annoted bed file with ensembltrans added with ChIPpeakAnno::addGeneannotated_peaks
#
# OUTPUT:
#	A data.frame containing every reads overlapping a list of enriched peaks with 3 columns:
#		* rname
#		* pos
#		* qwidth
extractsReadsInPeaks <- function(bam_file, annotated_peaks) {
	suppressMessages(library(Rsamtools))
	df <- data.frame()
	which <- GRanges(seqnames=Rle(annotated_peaks$space), ranges=IRanges(annotated_peaks$start, annotated_peaks$end),strand=Rle(annotated_peaks$strand))
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

# Extract reads from bam file that overlap with a genomic region
#
# INPUT:
#	bam_file:	Path to the bam file.
#	chr:		Current chromosome.
#	start:		Starting position of the current region.
#	end:		Ending position of the current region.
#	strand:		Strand of the current peak
#
# OUTPUT:
#	A data.frame containing every reads overlapping a list of enriched peaks with 3 columns:
#		* rname
#		* pos
#		* qwidth
extractsReadsInPeaksNoAnnotation <- function(bam_file, chr, start, end) {
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

# Create the base data.frame containing all the groups/id combinations. Only regions within max_distance are kept.
#
# INPUT:
# 	annotated_peaks:	Annoted bed file with ensembltrans added with ChIPpeakAnno::addGeneannotated_peaks
# 	groups: 		Data frame with every ensembl TSS and their status (TRUE or FALSE) for each group
# 	max_distance:		Maximal distance from TSS.
#
# OUTPUT:
#	A data.frame with one row per group/id combination.
initializeDf <- function(annotated_peaks, groups, max_distance) {
	all_groups <- data.frame()
	# 1. Remove annotated_peaks outside of max_distance
	distance_annotated_peaks <- annotated_peaks[which(abs(annotated_peaks["distancetoFeature"]) <= max_distance),]
	for(group in names(groups)[2:length(names(groups))]) {
		# 2. Keep only annotated_peaks found in current group
		current_group_annotated_peaks <- groups["ensembl_gene_id"][which(groups[group]==TRUE),]
#		print(length(current_group_annotated_peaks))
		indexes <- which(as.vector(distance_annotated_peaks$feature) %in% current_group_annotated_peaks)
#		print(length(indexes))
		if (any(indexes)) {
			current_annotated_peaks <- as.data.frame(distance_annotated_peaks[indexes,]$feature)
			colnames(current_annotated_peaks) <- c("ensembl_gene_id")
			current_annotated_peaks$group <- group
			all_groups <- rbind(all_groups, current_annotated_peaks)
		}
	}
	return(all_groups)
}

# Calculate the read intentity around TSS of every group/id combination.
#
# INPUT:
# 	annotated_peaks:	Annoted bed file with ensembltrans added with ChIPpeakAnno::addGeneannotated_peaks
#	bam_file: 		Path to the bam file.
#	initialized_df:		A data.frame with one row per group/id combination.
# 	max_distance:	Maximal distance from TSS.
#
# OUPUT:
#	A data.frame with the intentisties around TSS of every group/TSS combination.
parseBam <- function(annotated_peaks, bam_file, initialized_df, max_distance) {
	# 1. Initialize result data.frame
	result <- matrix(0, ncol=max_distance*2+1, nrow=nrow(initialized_df))
	# 2. For every TSS...
	for (gene_id in initialized_df$ensembl_gene_id) {
		print(paste("gene_id:", gene_id))
		result_index <- which(initialized_df$ensembl_gene_id == gene_id)
		current_peaks_indexes <- which(annotated_peaks$feature == gene_id)
		# 3. For every peaks in current TSS
		for (i in current_peaks_indexes) {
			# 4. Get current offset info
			current_TSS_offset <- annotated_peaks[i,]$start_position
			# 5. Extract reads 
			current_reads <- extractsReadsInPeaks(bam_file, annotated_peaks[i,])
			# 6. Increment the result data.frame for the current gene_id index
			vector_result <- convertReadsToPosVector(current_reads, current_TSS_offset, max_distance)
			result[result_index,] <- result[result_index,] + vector_result
		}
	}
	result <- as.data.frame(result)
	colnames(result) <- as.character(seq(-max_distance,max_distance))
	return(result)
}

# Calculate the read intentity around a list of regions
#
# INPUT:
#	chromosomes:	Vector with the names of the chromosomes
#	starts: 	Vector with starting values of the regions.
#	ends:		Vector with the ending of the regions.
#	bam_file: 	Path to the bam file.
#
# OUPUT:
#	A data.frame with the intentisties around features
#parseRegions <- function(chromosomes, starts, ends, bam_file) {
parseRegionsWithPadding <- function(chromosomes, starts, ends, bam_file) {
	percent_padding <- 0.2
	# 0. Find the median region size
	median_value <- median(mapply(function(x,y) max(x,y)-min(x,y), starts, ends))
	#domain_size <- median_value + median_value * percent_padding
	domain_size <- median_value + 2 * (median_value * percent_padding)
	# 1. Initialize result data.frame
	result <- matrix(0, ncol=domain_size, nrow=length(starts))

	# 2. For every region
	for(i in seq(1, length(starts))) {
		# 2.1 Parse values
		current_start <- min(starts[i], ends[i])
		current_start <- current_start - (current_start * percent_padding) # Add padding
		current_end <- max(starts[i], ends[i])
		current_end <- current_end + (current_end * percent_padding) # Add padding
		current_chr <- chromosomes[i]
		print(paste(current_chr, current_start, current_end, sep=' '))
		# 2.2 Extracts reads
		current_reads <- extractsReadsInPeaksNoAnnotation(bam_file, current_chr, current_start, current_end)
		# 2.3 Convert to relative position
		# Offset is the center position of the current region
		current_offset <- current_start + ((current_end - current_start) / 2)
		max_distance <- (current_end - current_start) / 2
		vector_result <- convertReadsToPosVector(current_reads, current_offset, max_distance)
		# 2.4 Scale the vector
		vector_result <- scaleVector(vector_result, domain_size)
		result[i,] <- vector_result
	}
	return(result)
}

# Calculate the read intentity around a list of regions (invert values when strand is negative)
# All regions must have the same size
#
# INPUT:
#	chromosomes:	Vector with the names of the chromosomes
#	starts: 	Vector with starting values of the regions.
#	ends:		Vector with the ending of the regions.
#	bam_file: 	Path to the bam file.
#
# OUPUT:
#	A data.frame with the intentisties around features
parseRegionsStrand <- function(chromosomes, starts, ends, strand, bam_file) {
	# 1. Initialize result data.frame
	result <- matrix(0, ncol=(ends[1]-starts[1])+1, nrow=length(starts))

	# 2. For every region
	for(i in seq(1, length(starts))) {
		# 2.1 Parse values
		current_start <- min(starts[i], ends[i])
		current_end <- max(starts[i], ends[i])
		current_chr <- chromosomes[i]
		current_strand <- strand[i]
		print(paste(current_chr, current_start, current_end, sep=' '))
		# 2.2 Extracts reads
		current_reads <- extractsReadsInPeaksNoAnnotation(bam_file, current_chr, current_start, current_end)
		if (length(current_reads) > 0) {
			# 2.3 Convert to relative position
			# Offset is the center position of the current region
			current_offset <- current_start + ((current_end - current_start) / 2)
			max_distance <- (current_end - current_start) / 2
			vector_result <- convertReadsToPosVector(current_reads, current_offset, max_distance)
			if (current_strand == "-1" | current_strand == -1 | current_strand == "-") {
				result[i,] <- rev(vector_result)
			} else {
				result[i,] <- vector_result
			}
		}
	}
	return(result)
}

parseRegionsAlternatePadding <- function(chromosomes, starts, ends, bam_file) {
	padding_value <- 2000
	# 0. Find the median region size
	median_value <- median(mapply(function(x,y) max(x,y)-min(x,y), starts, ends))
	domain_size <- median_value + padding_value * 2
	# 1. Initialize result data.frame
	result <- matrix(0, ncol=domain_size, nrow=length(starts))

	# 2. For every region
	for(i in seq(1, length(starts))) {
		# 2.1 Parse values
		current_start <- min(starts[i], ends[i])
		current_start <- current_start - padding_value # Add padding
		current_end <- max(starts[i], ends[i])
		current_end <- current_end + padding_value # Add padding
		current_chr <- chromosomes[i]
		print(paste(current_chr, current_start, current_end, sep=' '))
		# 2.2 Extracts reads
		current_reads <- extractsReadsInPeaksNoAnnotation(bam_file, current_chr, current_start, current_end)
		# 2.3 Convert to relative position
		# Offset is the center position of the current region
		current_offset <- current_start + ((current_end - current_start) / 2)
		max_distance <- (current_end - current_start) / 2
		vector_result <- convertReadsToPosVector(current_reads, current_offset, max_distance)
		# 2.4 Scale the vector
		vector_result <- scaleVector(vector_result, domain_size)
		result[i,] <- vector_result
	}
	return(result)
}

# Scale the values of a vector to fit with predetermined size
#
# INPUT:
#	values:	the values to scale
#	domain:	the range to fit the value to
#
# OUTPUT:
# 	A vector with the scaled data
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

# Calculate the read intentity around a list of regions
#
# INPUT:
# 	regions:	A data.frame with all the regions to parse.
#			Require the following column names:
#				* space:	name of the chromosome
#				* start: 	start of the region.
#				* end:		end of the region.
#			The following column name is facultative:
#				* strand:	strand of the region (default "+")
#	bam_file: 	Path to the bam file.
# 	max_distance:	Maximal distance from feature.
#
# OUPUT:
#	A data.frame with the intentisties around features
# TODO: Fill padding with NA instead of looking for the reads
parseBamNoAnnotation <- function(regions, bam_file, max_distance) {
	# 1. Initialize result data.frame
	result <- matrix(0, ncol=max_distance*2+1, nrow=nrow(regions))

	# 2. For every region
	for(i in seq(1, nrow(regions))) {
		# 2.1 Parse values
		current_start <- regions[i,]$start
		current_end <- regions[i,]$end
		current_chr <- regions[i,]$space
		current_strand <- regions[i,]$strand
		if (current_start > current_end) {
			tmp <- current_start
			current_start <- current_end
			current_end <- tmp
		}
		if (is.null(current_strand)) {
			current_strand <- "+"
		}
		# 2.2 Extracts reads
		padding <- ((max_distance*2+1) - (current_start+current_end)) / 2
		current_reads <- extractsReadsInPeaksNoAnnotation(bam_file, current_chr, current_start, current_end)
		# 2.3 Convert to relative position
		# Offset is the center position of the current region
		current_offset <- current_start + ((current_end - current_start) / 2)
		vector_result <- convertReadsToPosVector(current_reads, current_offset, max_distance)
		result[i,] <- vector_result
	}
	colnames(result) <- as.character(seq(-max_distance,max_distance))
	return(result)
}

# Calculate the read intensity around a list of regions. Regions are scaled to fit with the median.
#
# INPUT:
# 	regions:	A data.frame with all the regions to parse.
#			Require the following column names:
#				* space:	name of the chromosome
#				* start: 	start of the region.
#				* end:		end of the region.
#			The following column name is facultative:
#				* strand:	strand of the region (default "+")
#	bam_file: 	Path to the bam file.
#
# OUPUT:
#	A scaled data.frame with the intentisties around features
parseBamNoAnnotationScaled <- function(regions, bam_file) {
	# 1. Initialize result data.frame
	result <- matrix(0, ncol=max_distance*2+1, nrow=nrow(regions))

	# 2. For every region
	for(i in seq(1, nrow(regions))) {
		# 2.1 Parse values
		current_start <- regions[i,]$start
		current_end <- regions[i,]$end
		current_chr <- regions[i,]$space
		current_strand <- regions[i,]$strand
		if (current_start > current_end) {
			tmp <- current_start
			current_start <- current_end
			current_end <- tmp
		}
		if (is.null(current_strand)) {
			current_strand <- "+"
		}
		# 2.2 Extracts reads
		padding <- ((max_distance*2+1) - (current_start+current_end)) / 2
		current_reads <- extractsReadsInPeaksNoAnnotation(bam_file, current_chr, current_start, current_end)
		# 2.3 Convert to relative position
		# Offset is the center position of the current region
		current_offset <- current_start + ((current_end - current_start) / 2)
		vector_result <- convertReadsToPosVector(current_reads, current_offset, max_distance)
		result[i,] <- vector_result
	}
	colnames(result) <- as.character(seq(-max_distance,max_distance))
	return(result)
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

# Calculate the intensities of the reads in the enriched peaks for multiple groups.
#
# INPUT:
#	bam_file: 	Path to the bam file.
#	bed_file: 	Path to the bed file (obtained by doing peak calling on the bam file).
# 	max_distance:	Maximal distance from TSS.
#	group_filenames: A vector with all the file names for the groups.
#		Each file has a header, which correspond to the unique name of the group.
#		Each file contains a single column with the ensembl_transcript_id of the genes.
#
# OUPUT:
#	A data frame with every combination of group/id observed with the intensities of each
#	combination around the TSS.
calculateIntensities <- function(bam_file, bed_file, max_distance, group_filenames) {

}

#
## annotated_peaks: 	Annoted bed file with ensembltrans added with addGeneannotated_peaks
## gene_lists: 	Data frame with as many columns as there are groups to compare.
##		Header of each columns must be the name of the group.
## bam_file:	List of bam file, one per group
#calculate_intensity <- function(annotated_peaks, gene_lists, bam_file) {
#	# For each entry in annotated_peaks...
#	for (i in 1:nrow(annotated_peaks)) {
#		current_ensembl_id <- annotated_peaks$ensembl[i]
#		groups <- fetch_groups(current_ensembl_id, gene_lists)
#		if (length(groups) > 0) {
#		# If nearest TSS is in one of the groups
#			# If within 5kb of nearest TSS
#				# Extract reads from bam file in relevent groups
#		}
#	}
#}

addGroupInfos <- function(annotated_peaks, groups) {
#	c["ensembl_gene_id"][which(c["high_level_V8"]==TRUE),]
	# For every group in group list
	for(group in names(groups)[3:length(names(groups))]) {
		# Keep annotated_peaks that are found in current group
		current_annotation <- annotated_peaks
		current_group_annotated_peaks <- groups["ensembl_gene_id"][which(groups[group]==TRUE),]
		
		# Check valid ID
#		for(i in 1:nrow(annotated_peaks)) {
#			if (annotated_peaks[["ensembl_gene_id"]][i] %in% current_group_annotated_peaks) {
#				
#			}
	}
}
