## Test functions present in the bam_handler.R file

### {{{ --- Test setup ---

if(FALSE) {
  library( "RUnit" )
  library( "metagene" )
  library( "rtracklayer" )
}

### }}}

bam_files <- get_demo_bam_files()
named_bam_files <- bam_files
names(named_bam_files) <- paste("file", seq(1, length(bam_files)), sep = "_")
not_indexed_bam_file <- metagene:::get_not_indexed_bam_file()
regions <- lapply(metagene:::get_demo_regions(), rtracklayer::import)

###################################################
## Test the Bam_Handler$new() function (initialize)
###################################################

base_msg <- "Bam_Handler initialize - "

## Single valid bam file
test.bam_handler_single_valid_file <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files[1])
  msg <- paste0(base_msg,  
        "Valid initialize call with one bam file did not return correct class")
  checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")), msg = msg)
}

## Invalid bam file - not indexed
test.bam_handler_not_indexed_single_bam_file <- function() {
  obs <- tryCatch(metagene:::Bam_Handler$new(not_indexed_bam_file), 
                    error = conditionMessage)
  exp <- "All BAM files must be indexed"
  msg <- paste0(base_msg, 
        "Single not indexed base file did not return the correct error message")
  checkIdentical(obs, exp, msg)
}

## Multiple bam files, one not indexed
test.bam_handler_multiple_bam_file_one_not_indexed <- function() {
  one_bam_file_not_indexed <- c(bam_files, not_indexed_bam_file)
  obs <- tryCatch(metagene:::Bam_Handler$new(one_bam_file_not_indexed), 
                    error = conditionMessage)
  exp <- "All BAM files must be indexed"
  msg <- paste0(base_msg, 
        "Single not indexed base file did not return the correct error message")
  checkIdentical(obs, exp, msg)
}

## Multiple valid bam files, no cores
test.bam_handler_valid_files_no_cores <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  msg <- paste0(base_msg, "Valid initialize call with multiple bam files ", 
                    "did not return correct class")
  checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")), msg = msg)
}

## Unamed bam files
test.bam_handler_unamed_bam_files <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  obs <- rownames(bam_handler$get_bam_files())
  exp <- tools::file_path_sans_ext(basename(bam_files))
  msg <- paste(base_msg, "Valid initialize call with unnamed bam files ",
                    "did not return correct values.")
  checkIdentical(obs, exp, msg)
}

## Named bam files
test.bam_handler_named_bam_files <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = named_bam_files)
  obs <- rownames(bam_handler$get_bam_files())
  exp <- paste("file", seq(1, length(bam_files)), sep = "_")
  msg <- paste(base_msg, "Valid initialize call with named bam files did ",
                    "not return correct values.")
  checkIdentical(obs, exp, msg)
}

## Valid bam files, numeric cores
test.bam_handler_valid_files_numeric_cores <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files, cores = 2)
  msg <- paste0(base_msg, " Valid initialize call with multiple bam files and ", 
               "numeric cores did not return correct class")
  checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")), msg = msg)
}

## Valid bam files, bpparam cores
test.bam_handler_valid_files_bpparam_cores <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files, 
                                    cores = BiocParallel::MulticoreParam(workers = 2))
  checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")),
            msg = paste0(base_msg, "Valid initialize call with multiple bam ",
                "files and bpparam cores did not return correct class"))
}

## Zero core should not be accepted as an argument
test.bam_handler_initialize_zero_core_number<- function() {
  obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files, cores = 0), 
                  error=conditionMessage)
  exp <- "cores must be a positive numeric or BiocParallelParam instance"
  msg <- paste0(base_msg, "A negative core number argument did not ",
                "generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Negative integer core should not be accepted as an argument
test.bam_handler_initialize_negative_core_number<- function() {
    obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files, 
                    cores = -1), error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    msg <- paste0(base_msg, " A negative core number argument did ",
                  "not generate an exception with expected message.")
    checkIdentical(obs, exp, msg)
}

## Something other than an integer number should not be accepted 
## as an core number argument
test.bam_handler_initialize_not_integer_core_number<- function() {
  obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files,  
            cores = 2.22), error=conditionMessage)
  exp <- "cores must be a positive numeric or BiocParallelParam instance"
  msg <- paste0(base_msg,  "A decimal core number argument did not generate ",
                "an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Something other than an integer number should not be accepted as 
## an core number argument
test.bam_handler_initialize_string_core_number<- function() {
  obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = bam_files,  
            cores ="NotAInteger"), error=conditionMessage)
  exp <- "cores must be a positive numeric or BiocParallelParam instance"
  msg <- paste0(base_msg, "A generic text used as a core number did not ", 
            "generate an exception with expected message." )
  checkIdentical(obs, exp, msg)
}

## All BAM files must be in string format
test.bam_handler_initialize_file_name_not_in_string_format<- function() {
  obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = c(1,2)), 
            error = conditionMessage)
  exp <- "bam_files must be a vector of BAM filenames"
  msg <- paste0(base_msg, "Integers used as files argument did not generate ", 
            "an exception with expected message.")
  checkEquals(obs, exp, msg)
}

## All bam files must exist
test.bam_handler_initialize_with_not_existing_files<- function() {
  obs <- tryCatch(metagene:::Bam_Handler$new(bam_files = 
    c("NotExistingFile", "NotExistingFile2")), error = conditionMessage)
  exp <- "At least one BAM file does not exist"
  msg <- paste0(base_msg, "Not existing BAM file used as ", 
         "argument did not generate an exception with expected message.")
  checkEquals(obs, exp, msg)
}

###################################################
## Test the bam_handler$get_aligned_count()
###################################################

# Note, count were obtained with "samtools view -c -F0x4 ${file}"
exp <- list(4635, 1896, 956, 1999)
names(exp) <- bam_files

## Valid use case
test.bam_handler_get_aligned_count_valid_case <- function() {
  # Note, count were obtained with "samtools view -c -F0x4 ${file}"
  exp <- list(4635, 1896, 956, 1999)
  names(exp) <- bam_files
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  obs <- lapply(bam_files, bam_handler$get_aligned_count)
  names(obs) <- bam_files
  msg <- paste0(base_msg, "Bam_Handler get_aligned_count valid case does ", 
                "not give the expected aligned count")
  checkTrue(all(mapply("==", obs, exp)), msg = msg)
}

## Valid use case, multicore
test.bam_handler_get_aligned_count_valid_case_multicore <- function() {
  # Note, count were obtained with "samtools view -c -F0x4 ${file}"
  exp <- list(4635, 1896, 956, 1999)
  names(exp) <- bam_files
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files, cores = 2)
  obs <- lapply(bam_files, bam_handler$get_aligned_count)
  names(obs) <- bam_files
  msg <- paste0(base_msg, "Bam_Handler get_aligned_count valid case does ", 
                "not give the expected aligned count")
  checkTrue(all(mapply("==", obs, exp)), msg = msg)
}

## Invalid bam file
test.bam_handler_get_aligned_count_invalid_bam_file <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  obs <- tryCatch(bam_handler$get_aligned_count(bam_file = 
                            "not_a_valid_bam_file"), error=conditionMessage)
  exp <- "Bam file not_a_valid_bam_file not found."
  msg <- paste0(base_msg, "Bam_Handler get_aligned_count invalid bam file ", 
                "does not give the expected error message")
  checkEquals(obs, exp, msg = msg)
}

###################################################
## Test the bam_handler$get_rpm_coefficient()
###################################################

## Valid use case
test.bam_handler_get_rpm_coefficient_valid_case <- function() {
  # Note, count were obtained with "samtools view -c -F0x4 ${file}"
  exp <- list(4635/1000000, 1896/1000000, 956/1000000, 1999/1000000)
  names(exp) <- bam_files
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  obs <- lapply(bam_files, bam_handler$get_rpm_coefficient)
  names(obs) <- bam_files
  checkTrue(all(mapply("==", obs, exp)),
            msg = "Bam_Handler get_aligned_count valid case does not give the expected aligned count")
}

## Valid use case, multicore
test.bam_handler_get_rpm_coefficient_valid_case_multicore <- function() {
  # Note, count were obtained with "samtools view -c -F0x4 ${file}"
  exp <- list(4635/1000000, 1896/1000000, 956/1000000, 1999/1000000)
  names(exp) <- bam_files
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files, cores = 2)
  obs <- lapply(bam_files, bam_handler$get_rpm_coefficient)
  names(obs) <- bam_files
  checkTrue(all(mapply("==", obs, exp)),
            msg = "Bam_Handler get_aligned_count valid case does not give the expected aligned count")
}

## Invalid bam file
test.bam_handler_get_rpm_coefficient_invalid_bam_file <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  obs <- tryCatch(bam_handler$get_rpm_coefficient(bam_file = "not_a_valid_bam_file"), error=conditionMessage)
  exp <- "Bam file not_a_valid_bam_file not found."
  checkEquals(obs, exp,
              msg = "Bam_Handler get_aligned_count invalid bam file does not give the expected error message")
}

###################################################
## Test the bam_handler$get_normalized_coverage()
###################################################

base_msg <- "Bam_Handler get_normalized_coverage -"

## Valid use
test.bam_handler_get_normalized_coverage_valid_use <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  region <- regions[[1]][1]
  bam_file <- bam_files[1]
  obs <- bam_handler$get_normalized_coverage(bam_file, region)
  weight <- 1 / (bam_handler$get_aligned_count(bam_file) / 1000000)
  exp <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(bam_file, param = Rsamtools::ScanBamParam(which = reduce(region))))[region] * weight
  msg <- paste(base_msg, "Valid use do not give the expected results")
  checkIdentical(obs, exp, msg)
}

## Duplicated regions
test.bam_handler_get_normalized_coverage_duplicated_regions <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  region <- regions[[1]]
  bam_file <- bam_files[1]
  obs <- bam_handler$get_normalized_coverage(bam_file, region)
  msg <- paste(base_msg, "Duplicated regions do not give expected results")
  weight <- 1 / (bam_handler$get_aligned_count(bam_file) / 1000000)
  exp <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(bam_file, param = Rsamtools::ScanBamParam(which = reduce(region))))[region] * weight
  checkIdentical(obs, exp, msg)
  # Sanity check
  sane <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(bam_file, param = Rsamtools::ScanBamParam(which = region)))[region] * weight
  msg <- paste(base_msg, "Duplicated regions did not pass sanity test - reduced regions")
  checkTrue(!identical(reduce(region), region), msg = msg)
  msg <- paste(base_msg, "Duplicated regions did not pass sanity test - identical")
  checkTrue(!identical(obs, sane), msg = msg)
}

## Negative coverage
test.bam_handler_get_normalized_coverage_negative_coverage <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  region <- regions[[1]][1]
  bam_file <- bam_files[4]
  obs <- bam_handler$get_normalized_coverage(bam_file, region)
  msg <- paste(base_msg, "Negative coverage returns some negative values")
  checkTrue(all(sapply(obs, function(x) all(as.numeric(x) >= 0))), msg = msg)
  # Sanity check
  weight <- 1 / (bam_handler$get_aligned_count(bam_file) / 1000000)
  sane <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(bam_file, param = Rsamtools::ScanBamParam(which = reduce(region))), weight = weight)[region]
  msg <- paste(base_msg, "Negative coverage did not pass sanity test")
  checkTrue(!all(sapply(sane, function(x) all(as.numeric(x) >= 0))), msg = msg)
}

## Invalid bam file
test.bam_handler_get_normalized_coverage_invalid_bam_file <- function() {
  region <- regions[1]
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = "not_a_valid_bam_file", regions = region), error = conditionMessage)
  exp <- "Bam file not_a_valid_bam_file not found."
  msg <- paste(base_msg, "Invalid bam file did not give the expected error message.")
  checkIdentical(obs, exp, msg)
}

## Invalid regions class
test.bam_handler_get_normalized_coverage_invalid_regions_class <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  bam_file <- bam_files[1]
  obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file, regions = "not_a_valid_region"), error = conditionMessage)
  exp <- "Parameter regions must be a GRanges object."
  msg <- paste(base_msg, "Invalid regions class file did not give the expected error message.")
  checkIdentical(obs, exp, msg)
}

## Invalid regions length
test.bam_handler_get_normalized_coverage_invalid_regions_length <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  bam_file <- bam_files[1]
  obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file, regions = GenomicRanges::GRanges()), error = conditionMessage)
  exp <- "Parameter regions must not be an empty GRanges object"
  msg <- paste(base_msg, "Invalid regions length did not give the expected error message.")
  checkIdentical(obs, exp, msg)
}

## Invalid regions levels
test.bam_handler_get_normalized_coverage_invalid_regions_length <- function() {
  bam_handler <- metagene:::Bam_Handler$new(bam_files = bam_files)
  bam_file <- bam_files[1]
  region <- regions[[1]]
  seqlevels(region) <- c(seqlevels(region), "invalid_level")
  obs <- tryCatch(bam_handler$get_normalized_coverage(bam_file = bam_file, regions = region), error = conditionMessage)
  exp <- "Some seqlevels of regions are absent in bam_file header"
  msg <- paste(base_msg, "Invalid regions levels did not give the expected error message.")
  checkIdentical(obs, exp, msg)
}
