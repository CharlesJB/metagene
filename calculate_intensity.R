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
	library(Rsamtools)
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
#d[as.character(c(-1,0,1))][1,] <- d[as.character(c(-1,0,1))][1,] + 1	
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
			vector_result <- numeric(max_distance*2+1)
			# 4. Get current offset info
			current_TSS_offset <- annotated_peaks[i,]$start_position
			# 5. Extract reads 
			current_reads <- extractsReadsInPeaks(bam_file, annotated_peaks[i,])
			# 6. Increment the result data.frame for the current gene_id index
			if (nrow(current_reads) > 0) {
				positions <- unlist(mapply(function(x,y) seq(x, x+y), current_reads$pos - current_TSS_offset, current_reads$qwidth))
				positions <- positions[abs(positions)<=max_distance] # to remove reads beyond max distance
				positions <- positions + max_distance
				# TODO: faire les operations sur un vecteur puis incrementer la matrice en une seule etape
#				vector_result <- vector_result + tabulate(positions, nbins=max_distance*2+1)
				if (length(positions) > 0) {
					for (i in positions) {
						vector_result[i] <- vector_result[i] + 1
					}
				}
			}
			result[result_index,] <- result[result_index,] + vector_result
                }
	}
	result <- as.data.frame(result)
	colnames(result) <- as.character(seq(-max_distance,max_distance))
	return(result)
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
