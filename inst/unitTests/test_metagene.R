## Test functions present in the metagene.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
}

### }}}

bam_files <- get_demo_bam_files()
named_bam_files <- bam_files
not_indexed_bam_file <- metagene:::get_not_indexed_bam_file()

###################################################
## Test the metagene$new() function (initialize)
###################################################

base_msg <- "metagene initialize - "

## Invalid verbose value
test.metagene_initialize_invalid_verbose_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(verbose="ZOMBIES"), 
                    error=conditionMessage)
    exp <- "verbose must be a logicial value (TRUE or FALSE)"
    msg <- paste0(base_msg, "An invalid verbose value did not generate an ", 
                "exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Negative padding_size value
test.metagene_initialize_negative_padding_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(padding_size=-1), 
                                            error=conditionMessage)
    exp <- "padding_size must be a non-negative integer"
    msg <- paste0(base_msg, "A negative padding_size value did not generate ",
                "an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Non-integer padding_size value
test.metagene_initialize_invalid_string_padding_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(padding_size="NEW_ZOMBIE"), 
                                            error=conditionMessage)
    exp <- "padding_size must be a non-negative integer"
    msg <- paste0(base_msg, "A character padding_size value did not ", 
                    "generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Numerical padding_size value
test.metagene_initialize_invalid_numerical_padding_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(padding_size=1.2), 
                    error=conditionMessage)
    exp <- "padding_size must be a non-negative integer"
    msg <- paste0(base_msg, "A non-integer padding_size value did not ", 
                  "generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Negative padding_size value
test.metagene_initialize_negative_padding_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(core=-1), 
                    error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    msg <- paste0(base_msg, "A negative core value did not generate ",
                  "an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Non-integer core value
test.metagene_initialize_invalid_string_core_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(core="ZOMBIE2"), 
                    error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    msg <- paste0(base_msg, "A character core value did not ", 
                  "generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Numerical core value
test.metagene_initialize_invalid_numerical_core_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(core=1.2), 
                    error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    msg <- paste0(base_msg, "A non-integer core value did not ", 
                  "generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Zero core value
test.metagene_initialize_invalid_zero_core_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(core=0), 
                    error=conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    msg <- paste0(base_msg, "A zero core value did not ", 
                  "generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Non-character vector bam_files value
test.metagene_initialize_invalid_num_vector_bam_files_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(bam_files=c(2,4,3)), 
                    error=conditionMessage)
    exp <- "bam_files must be a vector of BAM filenames"
    msg <- paste0(base_msg, "A non-character vector bam_files value did not ", 
                  "generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

## Non-vector bam_files value
test.metagene_initialize_invalid_list_bam_files_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(bam_files=list(a="ZOMBIE_01.txt",
            b="ZOMBIE_02.txt")), error=conditionMessage)
    exp <- "bam_files must be a vector of BAM filenames"
    msg <- paste0(base_msg, "A non-vector bam_files value did not ", 
                  "generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Not indexed bam in bam_files value
test.metagene_initialize_invalid_no_index_bam_files_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(bam_files=not_indexed_bam_file), 
            error=conditionMessage)
    exp <- "All BAM files must be indexed"
    msg <- paste0(base_msg, "A not indexed BAM in bam_files value value ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Multiple bam files, only one not indexed in bam_files value
test.metagene_initialize_multiple_bam_file_one_not_indexed <- function() {
    obs <- tryCatch(metagene:::metagene$new(bam_files = c(bam_files, 
            not_indexed_bam_file)), error = conditionMessage)
    exp <- "All BAM files must be indexed"
    msg <- paste0(base_msg, "Only one not indexed BAM in bam_files value value",
                  " did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Not valid object in region value 
test.metagene_initialize_invalid_array_region_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(bam_files=named_bam_files, 
            region=array(data = NA, dim = c(2,2,2))), 
            error=conditionMessage)
    exp <- paste0("regions must be either a vector of BED filenames, a ",
            "GRanges object or a GrangesList object")
    msg <- paste0(base_msg, "A not indexed bam in bam_files value value ",
            "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 