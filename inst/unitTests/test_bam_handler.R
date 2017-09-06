## Test functions present in the bam_handler.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
    library( "rtracklayer" )
    library( "DBChIP" )
}

### }}}

bam_files <- get_demo_bam_files()
named_bam_files <- bam_files
names(named_bam_files) <- paste("file", seq(1, length(bam_files)), sep = "_")
not_indexed_bam_file <- metagene:::get_not_indexed_bam_file()
different_seqnames <- c(bam_files[1],
                        metagene:::get_different_seqnames_bam_file())
regions <- lapply(metagene:::get_demo_regions(), rtracklayer::import)

demo_bh <- metagene:::Bam_Handler$new(bam_files)
demo_bh_one <- metagene:::Bam_Handler$new(bam_files[1])

###################################################
## Test the Bam_Handler$new() function (initialize)
###################################################

## Single valid bam file
test.bam_handler_single_valid_file <- function() {
    bam_handler <- demo_bh_one$clone()
    checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")))
}

## Different seqnames bam file warning
test.bam_handler_different_seqnames_bam_file_warning <- function() {
    bam_handler <- demo_bh_one$clone()
	msg <- "\n\nSome bam files have discrepancies in their "
    msg <- paste0(msg, "seqnames.")
    msg <- paste0(msg, "\n\n")
    msg <- paste0(msg, "This could be caused by chromosome names")
    msg <- paste0(msg, " present only in a subset of the bam ")
    msg <- paste0(msg, "files (i.e.: chrY in some bam files, but ")
    msg <- paste0(msg, "absent in others.\n\n")
    msg <- paste0(msg, "This could also be caused by ")
    msg <- paste0(msg, "discrepancies in the seqlevels style")
    msg <- paste0(msg, " (i.e.: UCSC:chr1 versus NCBI:1)\n\n")
    obs <- tryCatch(metagene:::Bam_Handler$new(different_seqnames),
                    warning = conditionMessage)
    checkIdentical(obs, exp)
}

## Invalid bam file - not indexed
test.bam_handler_not_indexed_single_bam_file <- function() {
    obs <- tryCatch(metagene:::Bam_Handler$new(not_indexed_bam_file),
                    error = conditionMessage)
    exp <- "All BAM files must be indexed"
    checkIdentical(obs, exp)
}

## Multiple bam files, one not indexed
test.bam_handler_multiple_bam_file_one_not_indexed <- function() {
    one_bam_file_not_indexed <- c(bam_files, not_indexed_bam_file)
    obs <- tryCatch(metagene:::Bam_Handler$new(one_bam_file_not_indexed),
                    error = conditionMessage)
    exp <- "All BAM files must be indexed"
    checkIdentical(obs, exp)
}

## Unnamed bam files
test.bam_handler_unamed_bam_files <- function() {
    bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
    obs <- rownames(bam_handler$get_bam_files())
    exp <- tools::file_path_sans_ext(basename(bam_files))
    checkIdentical(obs, exp)
}

## Named bam files
test.bam_handler_named_bam_files <- function() {
    bam_handler <- metagene:::Bam_Handler$new(bam_files = named_bam_files)
    obs <- rownames(bam_handler$get_bam_files())
    exp <- paste("file", seq(1, length(bam_files)), sep = "_")
    checkIdentical(obs, exp)
}

## Valid bam files, numeric cores
test.bam_handler_valid_files_numeric_cores <- function() {
    bam_handler <- metagene:::demo_bh_multicore$clone()
    checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")))
}

## Valid bam files, bpparam cores
test.bam_handler_valid_files_bpparam_cores <- function() {
    cores <- BiocParallel::SnowParam(workers = 2)
    bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files[1],
                                              cores = cores)
    checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")))
}

## Zero core should not be accepted as an argument
test.bam_handler_initialize_zero_core_number<- function() {
    obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files,
                                               cores = 0),
                    error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    checkIdentical(obs, exp)
}

## Negative integer core should not be accepted as an argument
test.bam_handler_initialize_negative_core_number<- function() {
    obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files,
                                               cores = -1),
                    error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    checkIdentical(obs, exp)
}

## Something other than an integer number should not be accepted
## as an core number argument
test.bam_handler_initialize_not_integer_core_number<- function() {
    obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files,
                                               cores = 2.22),
                    error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    checkIdentical(obs, exp)
}

## Something other than an integer number should not be accepted as
## an core number argument
test.bam_handler_initialize_string_core_number<- function() {
    obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files,
                                               cores ="NotAInteger"),
                    error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    checkIdentical(obs, exp)
}

## All BAM files must be in string format
test.bam_handler_initialize_file_name_not_in_string_format<- function() {
    obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = c(1,2)),
                    error = conditionMessage)
    exp <- "bam_files must be a vector of BAM filenames"
    checkEquals(obs, exp)
}

## All bam files must exist
test.bam_handler_initialize_with_not_existing_files<- function() {
    bam_files <- c("NotExistingFile", "NotExistingFile2")
    obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files),
                    error = conditionMessage)
    exp <- "At least one BAM file does not exist"
    checkEquals(obs, exp)
}

###################################################
## Test the bam_handler$get_aligned_count()
###################################################

## Valid use case
test.bam_handler_get_aligned_count_valid_case <- function() {
    # Note, count were obtained with "samtools view -c -F0x4 ${file}"
    exp <- list(4635, 1896, 956, 1999, 6025)
    names(exp) <- bam_files
    bam_handler <- demo_bh$clone()
    obs <- lapply(bam_files, bam_handler$get_aligned_count)
    names(obs) <- bam_files
    checkTrue(all(mapply("==", obs, exp)))
}

## Valid use case, multicore
test.bam_handler_get_aligned_count_valid_case_multicore <- function() {
    # Note, count were obtained with "samtools view -c -F0x4 ${file}"
    exp <- list(4635, 1896, 956, 1999, 6025)
    names(exp) <- bam_files
    bam_handler <- metagene:::demo_bh_multicore$clone()
    obs <- lapply(bam_files, bam_handler$get_aligned_count)
    names(obs) <- bam_files
    checkTrue(all(mapply("==", obs, exp)))
}

## Invalid bam file
test.bam_handler_get_aligned_count_invalid_bam_file <- function() {
    bam_handler <- demo_bh$clone()
    bam_file <- "not_a_valid_bam_file"
    obs <- tryCatch(bam_handler$get_aligned_count(bam_file = bam_file),
                    error=conditionMessage)
    exp <- "Bam file not_a_valid_bam_file not found."
    checkEquals(obs, exp)
}

###################################################
## Test the bam_handler$get_rpm_coefficient()
###################################################

## Valid use case
test.bam_handler_get_rpm_coefficient_valid_case <- function() {
    # Note, count were obtained with "samtools view -c -F0x4 ${file}"
    exp <- list(4635/1000000, 1896/1000000, 956/1000000, 1999/1000000,
                6025/1000000)
    names(exp) <- bam_files
    bam_handler <- demo_bh$clone()
    obs <- lapply(bam_files, bam_handler$get_rpm_coefficient)
    names(obs) <- bam_files
    checkTrue(all(mapply("==", obs, exp)))
}

## Valid use case, multicore
test.bam_handler_get_rpm_coefficient_valid_case_multicore <- function() {
    # Note, count were obtained with "samtools view -c -F0x4 ${file}"
    exp <- list(4635/1000000, 1896/1000000, 956/1000000, 1999/1000000,
                6025/1000000)
    names(exp) <- bam_files
    bam_handler <- metagene:::demo_bh_multicore$clone()
    obs <- lapply(bam_files, bam_handler$get_rpm_coefficient)
    names(obs) <- bam_files
    checkTrue(all(mapply("==", obs, exp)))
}

## Invalid bam file
test.bam_handler_get_rpm_coefficient_invalid_bam_file <- function() {
    bam_handler <- demo_bh$clone()
    bam_file <- "not_a_valid_bam_file"
    obs <- tryCatch(bam_handler$get_rpm_coefficient(bam_file = bam_file),
                    error=conditionMessage)
    exp <- "Bam file not_a_valid_bam_file not found."
    checkEquals(obs, exp)
}

###################################################
## Test the bam_handler$get_coverage()
###################################################

## Valid use
test.bam_handler_get_coverage_valid_use <- function() {
    bam_handler <- demo_bh$clone()
    region <- regions[[1]][1]
    bam_file <- bam_files[1]
    obs <- bam_handler$get_coverage(bam_file, region)
    param <- Rsamtools::ScanBamParam(which = reduce(region))
    exp <- GenomicAlignments::readGAlignments(bam_file, param = param)
    exp <- GenomicAlignments::coverage(exp)
    checkIdentical(obs, exp)
}

## Multicore
test.bam_handler_get_coverage_multicore <- function() {
    bam_file <- bam_files[1]
    region <- regions[[1]]
    bam_handler <- metagene:::demo_bh_multicore$clone()
    coverages <- bam_handler$get_coverage(bam_file, region)
    checkEquals(length(coverages), 22)
}

## Multiple chromosomes
test.bam_handler_get_coverage_multiple_chromosomes <- function() {
    bam_file <- metagene:::get_coverage_bam_file()
    region <- rtracklayer::import(metagene:::get_coverage_region())
    param <- Rsamtools::ScanBamParam(which = GenomicRanges::reduce(region))
    exp <- GenomicAlignments::readGAlignments(bam_file, param = param)
    count <- Rsamtools::countBam(bam_file)$records
    exp <- GenomicAlignments::coverage(exp)
    bam_handler <- metagene:::Bam_Handler$new(bam_file)
    obs <- bam_handler$get_coverage(bam_file, region)
    checkTrue(all(sum(exp - obs) == 0))
    checkTrue(identical(names(exp), names(obs)))
}

## Duplicated regions
test.bam_handler_get_coverage_duplicated_regions <- function() {
    bam_handler <- demo_bh$clone()
    region <- regions[[1]]
    bam_file <- bam_files[1]
    obs <- bam_handler$get_coverage(bam_file, region)
    param <- Rsamtools::ScanBamParam(which = reduce(region))
    exp <- GenomicAlignments::readGAlignments(bam_file, param = param)
    exp <- GenomicAlignments::coverage(exp)
    checkIdentical(obs, exp)
    # Sanity check
    param <- Rsamtools::ScanBamParam(which = region)
    sane <- GenomicAlignments::readGAlignments(bam_file, param = param)
    sane <- GenomicAlignments::coverage(sane)
    checkTrue(!identical(reduce(region), region))
    checkTrue(!identical(obs, sane))
}

## Negative coverage
test.bam_handler_get_coverage_negative_coverage <- function() {
    bam_handler <- demo_bh$clone()
    region <- regions[[1]][1]
    bam_file <- bam_files[4]
    obs <- bam_handler$get_coverage(bam_file, region)
    checkTrue(all(vapply(obs, function(x) all(runValue(x) >= 0), logical(1))))
}

## Invalid bam file
test.bam_handler_get_coverage_invalid_bam_file <- function() {
    region <- regions[1]
    bam_handler <- demo_bh$clone()
    obs <- tryCatch(bam_handler$get_coverage(bam_file = "not_a_valid_bam_file",
                                             regions = region),
                    error = conditionMessage)
    exp <- "Bam file not_a_valid_bam_file not found."
    checkIdentical(obs, exp)
}

## Invalid regions class
test.bam_handler_get_coverage_invalid_regions_class <- function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- "not_a_valid_region"
    obs <- tryCatch(bam_handler$get_coverage(bam_file = bam_file,
                                             regions = region),
                    error = conditionMessage)
    exp <- "Parameter regions must be a GRanges object."
    checkIdentical(obs, exp)
}

## Invalid regions length
test.bam_handler_get_coverage_invalid_regions_length <- function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- GenomicRanges::GRanges()
    obs <- tryCatch(bam_handler$get_coverage(bam_file = bam_file,
                                             regions = region),
                    error = conditionMessage)
    exp <- "Parameter regions must not be an empty GRanges object"
    checkIdentical(obs, exp)
}

## All seqnames not in bam
test.bam_handler_get_coverage_invalid_regions_all_seqnames_not_in_bam <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    # All seqlevels
    seqlevels(region) <- "invalid_level"
    obs <- tryCatch(bam_handler$get_coverage(bam_file = bam_file,
                                             regions = region),
                    error = conditionMessage)
    exp <- "Some seqlevels of regions are absent in bam_file"
    checkIdentical(obs, exp)
}

## All seqnames not in bam force seqlevels
test.bam_handler_get_coverage_all_seqnames_not_in_bam_force_seqlevels <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    # All seqlevels
    seqlevels(region) <- "invalid_level"
    obs <- tryCatch(bam_handler$get_coverage(bam_file = bam_file,
                                             regions = region,
                                             force_seqlevels = TRUE),
                    error = conditionMessage)
    exp <- "No seqlevels matching between regions and bam file"
    checkIdentical(obs, exp)
}

## One Seqnames not in bam
test.bam_handler_get_coverage_invalid_regions_one_seqnames_not_in_bam <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    seqlevels(region) <- c(seqlevels(region), "invalid_level")
    seqnames(region)[1] <- "invalid_level"
    obs <- tryCatch(bam_handler$get_coverage(bam_file = bam_file,
                                             regions = region),
                    error = conditionMessage)
    exp <- "Some seqlevels of regions are absent in bam_file"
    checkIdentical(obs, exp)
}

## Invalid regions seqlevels not in bam
test.bam_handler_get_coverage_invalid_regions_seqlevels_not_in_bam <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    seqlevels(region) <- c(seqlevels(region), "invalid_level")
    obs <- tryCatch(bam_handler$get_coverage(bam_file = bam_file,
                                             regions = region),
                    error = conditionMessage)
    exp <- "Some seqlevels of regions are absent in bam_file"
    checkIdentical(obs, exp)
}

## Seqnames not in bam force seqlevels
test.bam_handler_get_coverage_seqnames_not_in_bam_force <- function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    seqlevels(region) <- c(seqlevels(region), "invalid_level")
    seqnames(region)[1] <- "invalid_level"
    obs <- tryCatch(bam_handler$get_coverage(bam_file = bam_file,
                                             regions = region,
                                             force_seqlevels = TRUE),
                    error = conditionMessage)
    checkTrue(class(obs) == "SimpleRleList")
}

###################################################
## Test the bam_handler$get_normalized_coverage()
###################################################

## Valid use
test.bam_handler_get_normalized_coverage_valid_use <- function() {
    bam_handler <- demo_bh$clone()
    region <- regions[[1]][1]
    bam_file <- bam_files[1]
    obs <- bam_handler$get_normalized_coverage(bam_file, region)
    weight <- 1 / (bam_handler$get_aligned_count(bam_file) / 1000000)
    param <- Rsamtools::ScanBamParam(which = reduce(region))
    exp <- GenomicAlignments::readGAlignments(bam_file, param = param)
    exp <- GenomicAlignments::coverage(exp) * weight
    checkIdentical(obs, exp)
}

## Multicore
test.bam_handler_get_normalized_coverage_multicore <- function() {
    bam_file <- bam_files[1]
    region <- regions[[1]]
    bam_handler <- metagene:::demo_bh_multicore$clone()
    coverages <- bam_handler$get_normalized_coverage(bam_file, region)
    checkEquals(length(coverages), 22)
}

## Multiple chromosomes
test.bam_handler_get_normalized_coverage_multiple_chromosomes <- function() {
    bam_file <- metagene:::get_coverage_bam_file()
    region <- rtracklayer::import(metagene:::get_coverage_region())
    param <- Rsamtools::ScanBamParam(which = GenomicRanges::reduce(region))
    exp <- GenomicAlignments::readGAlignments(bam_file, param = param)
    count <- Rsamtools::countBam(bam_file)$records
    weight <- weight <- 1 / (count / 1000000)
    exp <- GenomicAlignments::coverage(exp) * weight
    bam_handler <- metagene:::Bam_Handler$new(bam_file)
    obs <- bam_handler$get_normalized_coverage(bam_file, region)
    checkTrue(all(sum(exp - obs) == 0))
    checkTrue(identical(names(exp), names(obs)))
}

## Duplicated regions
test.bam_handler_get_normalized_coverage_duplicated_regions <- function() {
    bam_handler <- demo_bh$clone()
    region <- regions[[1]]
    bam_file <- bam_files[1]
    obs <- bam_handler$get_normalized_coverage(bam_file, region)
    weight <- 1 / (bam_handler$get_aligned_count(bam_file) / 1000000)
    param <- Rsamtools::ScanBamParam(which = reduce(region))
    exp <- GenomicAlignments::readGAlignments(bam_file, param = param)
    exp <- GenomicAlignments::coverage(exp) * weight
    checkIdentical(obs, exp)
    # Sanity check
    param <- Rsamtools::ScanBamParam(which = region)
    sane <- GenomicAlignments::readGAlignments(bam_file, param = param)
    sane <- GenomicAlignments::coverage(sane) * weight
    checkTrue(!identical(reduce(region), region))
    checkTrue(!identical(obs, sane))
}

## Negative coverage
test.bam_handler_get_normalized_coverage_negative_coverage <- function() {
    bam_handler <- demo_bh$clone()
    region <- regions[[1]][1]
    bam_file <- bam_files[4]
    obs <- bam_handler$get_normalized_coverage(bam_file, region)
    checkTrue(all(sapply(obs, function(x) all(runValue(x) >= 0))))
}

## Invalid bam file
test.bam_handler_get_normalized_coverage_invalid_bam_file <- function() {
    region <- regions[1]
    bam_handler <- demo_bh$clone()
    bam_file <- "not_a_valid_bam_file"
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region),
                    error = conditionMessage)
    exp <- "Bam file not_a_valid_bam_file not found."
    checkIdentical(obs, exp)
}

## Invalid regions class
test.bam_handler_get_normalized_coverage_invalid_regions_class <- function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- "not_a_valid_region"
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region),
                    error = conditionMessage)
    exp <- "Parameter regions must be a GRanges object."
    checkIdentical(obs, exp)
}

## Invalid regions length
test.bam_handler_get_normalized_coverage_invalid_regions_length <- function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- GenomicRanges::GRanges()
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region),
                    error = conditionMessage)
    exp <- "Parameter regions must not be an empty GRanges object"
    checkIdentical(obs, exp)
}

## All seqnames not in bam
test.bam_handler_get_normalized_coverage_invalid_all_seqnames_not_in_bam <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    # All seqlevels
    seqlevels(region) <- "invalid_level"
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region),
                    error = conditionMessage)
    exp <- "Some seqlevels of regions are absent in bam_file"
    checkIdentical(obs, exp)
}

## All seqnames not in bam force seqlevels
test.bam_handler_get_normalized_coverage_seqnames_not_in_bam_force_seqlevels <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    # All seqlevels
    seqlevels(region) <- "invalid_level"
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region,
                                                        force_seqlevels = TRUE),
                    error = conditionMessage)
    exp <- "No seqlevels matching between regions and bam file"
    checkIdentical(obs, exp)
}

## One Seqnames not in bam
test.bam_handler_get_normalized_coverage_invalid_regions_seqnames_not_in_bam <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    seqlevels(region) <- c(seqlevels(region), "invalid_level")
    seqnames(region)[1] <- "invalid_level"
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region),
                    error = conditionMessage)
    exp <- "Some seqlevels of regions are absent in bam_file"
    checkIdentical(obs, exp)
}

## Valid regions seqlevels not in bam
test.bam_handler_get_normalized_coverage_regions_seqlevels_not_in_bam_no_regions <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    seqlevels(region) <- c(seqlevels(region), "invalid_level")
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region),
                    error = conditionMessage)
    checkTrue(obs == "Some seqlevels of regions are absent in bam_file")
}

## Valid regions seqlevels not in bam force
test.bam_handler_get_normalized_coverage_regions_seqlevels_not_in_bam__no_regions_force <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    seqlevels(region) <- c(seqlevels(region), "invalid_level")
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region,
							force_seqlevels = TRUE),
                    error = conditionMessage)
    checkTrue(class(obs) == "SimpleRleList")
}

## Seqnames not in bam force seqlevels
test.bam_handler_get_normalized_coverage_seqnames_not_in_bam_force <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    seqlevels(region) <- c(seqlevels(region), "invalid_level")
    seqnames(region)[1] <- "invalid_level"
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region,
                                                        force_seqlevels = TRUE),
                    error = conditionMessage)
    checkTrue(class(obs) == "SimpleRleList")
}

## No matching seqnames force
test.bam_handler_get_normalized_coverage_no_matching_seqnames_force <-
    function() {
    bam_handler <- demo_bh$clone()
    bam_file <- bam_files[1]
    region <- regions[[1]]
    seqlevels(region) <- "invalid_level"
    obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file,
                                                        regions = region,
                                                        force_seqlevels = TRUE),
                    error = conditionMessage)
    checkTrue(obs == "No seqlevels matching between regions and bam file")
}

###################################################
## Test the bam_handler$get_noise_ratio()
###################################################

chip.bed <- system.file("extdata/align1_rep1.bed", package="metagene")
input.bed <- system.file("extdata/ctrl.bed", package="metagene")
chip.bam <- system.file("extdata/align1_rep1.bam", package="metagene")
chip.bam <- basename(tools::file_path_sans_ext(chip.bam))
input.bam <- system.file("extdata/ctrl.bam", package="metagene")
input.bam <- basename(tools::file_path_sans_ext(input.bam))

## Valid use
test.bam_handler_get_noise_ratio_valid_use <- function() {
    exp <- DBChIP:::NCIS(chip.bed, input.bed, "BED")$est
    bam_handler <- demo_bh$clone()
    obs <- bam_handler$get_noise_ratio(chip.bam, input.bam)
    checkIdentical(obs, exp)
}

## Invalid chip_bam_file
test.bam_handler_get_noise_ratio_invalid_chip_bam_file <- function() {
    bam_handler <- demo_bh$clone()
    chip_bam_name <- "not_a_valid_bam_file"
    obs <- tryCatch(bam_handler$get_noise_ratio(chip_bam_names = chip_bam_name,
                                                input_bam_names = input.bam),
                    error=conditionMessage)
    exp <- "Bam file not_a_valid_bam_file not found."
    checkEquals(obs, exp)
}

## Invalid input_bam_file
test.bam_handler_get_noise_ratio_invalid_input_bam_file <- function() {
    bam_handler <- demo_bh$clone()
    input_b_name <- "not_a_valid_bam_file"
    obs <- tryCatch(bam_handler$get_noise_ratio(chip_bam_names = chip.bam,
                                                input_bam_names = input_b_name),
                    error=conditionMessage)
    exp <- "Bam file not_a_valid_bam_file not found."
    checkEquals(obs, exp)
}

###################################################
## Test the bam_handler$get_bam_name()
###################################################

## Valid bam file
test.bam_handler_get_bam_name_valid_bam_file <- function() {
    bam_handler <- demo_bh$clone()
    obs <- bam_handler$get_bam_name(bam_files[1])
    exp <- tools::file_path_sans_ext(basename(bam_files[1]))
    checkIdentical(obs, exp)
}

## Valid bam name
test.bam_handler_get_bam_name_valid_bam_name <- function() {
    bam_handler <- demo_bh$clone()
    bam_name <- tools::file_path_sans_ext(basename(bam_files[1]))
    obs <- bam_handler$get_bam_name(bam_name)
    exp <- bam_name
    checkIdentical(obs, exp)
}

## Invalid name
test.bam_handler_get_bam_name_invalid_bam_name <- function() {
    bam_handler <- demo_bh$clone()
    obs <- bam_handler$get_bam_name("invalid_bam_name")
    exp <- NULL
    checkIdentical(obs, exp)
}
