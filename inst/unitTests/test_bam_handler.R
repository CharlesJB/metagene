## Test functions present in the bam_handler.R file

### {{{ --- Test setup ---

if(FALSE) {
  library( "RUnit" )
  library( "metagene" )
}

### }}}

bam_files <- get_demo_bam_files()
not_indexed_bam_file <- metagene:::get_not_indexed_bam_file()

###################################################
## Test the Bam_Handler$new() function (initialize)
###################################################

base_msg <- "Bam_Handler initialize -"

## Single valid bam file
test.bam_handler_single_valid_file <- function() {
  bam_handler <- Bam_Handler$new(bam_files[1])
  checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")),
            msg = "Bam_Handler initialize - Valid initialize call with one bam file did not return correct class")
}

## Invalid bam file - not indexed
test.bam_handler_not_indexed_single_bam_file <- function() {
  obs <- tryCatch(Bam_Handler$new(not_indexed_bam_file), error = conditionMessage)
  exp <- "All bam files must be indexed."
  msg <- paste(base_msg, "Single not indexed base file did not return the correct error message")
  checkIdentical(obs, exp, msg)
}

## Multiple bam files, one not indexed
test.bam_handler_multiple_bam_file_one_not_indexed <- function() {
  one_bam_file_not_indexed <- c(bam_files, not_indexed_bam_file)
  obs <- tryCatch(Bam_Handler$new(one_bam_file_not_indexed), error = conditionMessage)
  exp <- "All bam files must be indexed."
  msg <- paste(base_msg, "Single not indexed base file did not return the correct error message")
  checkIdentical(obs, exp, msg)
}

## Multiple valid bam files, no cores
test.bam_handler_valid_files_no_cores <- function() {
  bam_handler <- Bam_Handler$new(bam_files = bam_files)
  checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")),
    msg = "Bam_Handler initialize - Valid initialize call with multiple bam files did not return correct class")
}

## Valid bam files, numeric cores
test.bam_handler_valid_files_numeric_cores <- function() {
  bam_handler <- Bam_Handler$new(bam_files = bam_files, cores = 2)
  checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")),
            msg = "Bam_Handler initialize - Valid initialize call with multiple bam files and numeric cores did not return correct class")
}

## Valid bam files, bpparam cores
test.bam_handler_valid_files_numeric_cores <- function() {
  bam_handler <- Bam_Handler$new(bam_files = bam_files, cores = MulticoreParam(workers = 2))
  checkTrue(all(class(bam_handler) == c("Bam_Handler", "R6")),
            msg = "Bam_Handler initialize - Valid initialize call with multiple bam files and bpparam cores did not return correct class")
}

## Zero core should not be accepted as an argument
test.bam_handler_initialize_zero_core_number<- function() {
  obs <- tryCatch(Bam_Handler$new(bam_files = bam_files, cores = 0), error=conditionMessage)
  exp <- "cores must be positive numeric or BiocParallelParam instance."
  obs <- tryCatch(Bam_Handler$new(bam_files = bam_files, cores = -1), error=conditionMessage)
  exp <- "cores must be positive numeric or BiocParallelParam instance."
  checkIdentical(obs, exp,
    msg = "Bam_Handler initialize - A negative core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument
test.bam_handler_initialize_not_integer_core_number<- function() {
  obs <- tryCatch(Bam_Handler$new(bam_files = bam_files,  cores = 2.22), error=conditionMessage)
  exp <- "cores must be positive numeric or BiocParallelParam instance."
  checkIdentical(obs, exp,
    msg = "Bam_Handler initialize - A decimal core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument
test.bam_handler_initialize_string_core_number<- function() {
  obs <- tryCatch(Bam_Handler$new(bam_files = bam_files,  cores ="NotAInteger"), error=conditionMessage)
  exp <- "cores must be positive numeric or BiocParallelParam instance."
  checkIdentical(obs, exp,
    msg = "Bam_Handler initialize - A generic text used as a core number did not generate an exception with expected message.")
}

## All bam files must be in string format
test.bam_handler_initialize_file_name_not_in_string_format<- function() {
  obs <- tryCatch(Bam_Handler$new(bam_files = c(1,2)), error = conditionMessage)
  exp <- "At least one BAM file name is not a character string."
  checkEquals(obs, exp,
    msg = "Bam_Handler initialize - Integers used as files argument did not generate an exception with expected message.")
}

## All bam files must exist
test.bam_handler_initialize_not_existing_files<- function() {
  obs <- tryCatch(Bam_Handler$new(bam_files = c("NotExistingFile", "NotExistingFile2")), error = conditionMessage)
  exp <- "At least one BAM file does not exist."
  checkEquals(obs, exp,
    msg = "Bam_Handler initialize - Not existing BAM file used as argument did not generate an exception with expected message.")
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
  bam_handler <- Bam_Handler$new(bam_files = bam_files)
  obs <- lapply(bam_files, bam_handler$get_aligned_count)
  names(obs) <- bam_files
  checkTrue(all(mapply("==", obs, exp)),     
    msg = "Bam_Handler get_aligned_count valid case does not give the expected aligned count")
}

## Valid use case, multicore
test.bam_handler_get_aligned_count_valid_case_multicore <- function() {
  # Note, count were obtained with "samtools view -c -F0x4 ${file}"
  exp <- list(4635, 1896, 956, 1999)
  names(exp) <- bam_files
  bam_handler <- Bam_Handler$new(bam_files = bam_files, cores = 2)
  obs <- lapply(bam_files, bam_handler$get_aligned_count)
  names(obs) <- bam_files
  checkTrue(all(mapply("==", obs, exp)),
    msg = "Bam_Handler get_aligned_count valid case does not give the expected aligned count")
}

## Invalid bam file
test.bam_handler_get_aligned_count_invalid_bam_file <- function() {
  bam_handler <- Bam_Handler$new(bam_files = bam_files)
  obs <- tryCatch(bam_handler$get_aligned_count(bam_file = "not_a_valid_bam_file"), error=conditionMessage)
  exp <- "Bam file not_a_valid_bam_file not found."
  checkEquals(obs, exp,
    msg = "Bam_Handler get_aligned_count invalid bam file does not give the expected error message")
}

###################################################
## Test the bam_handler$get_rpm_coefficient()
###################################################

## Valid use case
test.bam_handler_get_rpm_coefficient_valid_case <- function() {
  # Note, count were obtained with "samtools view -c -F0x4 ${file}"
  exp <- list(4635/1000000, 1896/1000000, 956/1000000, 1999/1000000)
  names(exp) <- bam_files
  bam_handler <- Bam_Handler$new(bam_files = bam_files)
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
  bam_handler <- Bam_Handler$new(bam_files = bam_files, cores = 2)
  obs <- lapply(bam_files, bam_handler$get_rpm_coefficient)
  names(obs) <- bam_files
  checkTrue(all(mapply("==", obs, exp)),
            msg = "Bam_Handler get_aligned_count valid case does not give the expected aligned count")
}

## Invalid bam file
test.bam_handler_get_rpm_coefficient_invalid_bam_file <- function() {
  bam_handler <- Bam_Handler$new(bam_files = bam_files)
  obs <- tryCatch(bam_handler$get_rpm_coefficient(bam_file = "not_a_valid_bam_file"), error=conditionMessage)
  exp <- "Bam file not_a_valid_bam_file not found."
  checkEquals(obs, exp,
              msg = "Bam_Handler get_aligned_count invalid bam file does not give the expected error message")
}
