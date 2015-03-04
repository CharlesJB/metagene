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
regions <- metagene:::get_demo_regions()

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


###################################################
## Test the metagene$plot() function 
###################################################

base_msg <- "metagene plot - "

mg <- metagene:::metagene$new(bam_files=named_bam_files, regions=regions)

# Not valid design object
test.metagene_plot_invalid_design <- function() {
    obs <- tryCatch(mg$plot(design=c(1,2)), 
                    error=conditionMessage)
    exp <- "design must be a data.frame object"
    msg <- paste0(base_msg, "A vector design object ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 

# Design data.frame with not enough columns
test.metagene_plot_invalid_design_data_frame <- function() {
    obs <- tryCatch(mg$plot(design=data.frame(a=c("ZOMBIE_ONE", "ZOMBIE_TWO"))), 
                    error=conditionMessage)
    exp <- "design must have at least 2 columns"
    msg <- paste0(base_msg, "A design data.frame with only one column ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 

# Design data.frame with invalid first column
test.metagene_plot_invalid_design_first_column <- function() {
    obs <- tryCatch(mg$plot(design=data.frame(a=c(1,3), 
                zombies=c("ZOMBIE_ONE", "ZOMBIE_TWO"))), 
                error=conditionMessage)
    exp <- "The first column of design must be BAM filenames"
    msg <- paste0(base_msg, "A design data.frame with numbers in first column ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 

# Design data.frame with invalid second column
test.metagene_plot_invalid_design_second_column <- function() {
    designTemp<-data.frame(a=named_bam_files, 
                           zombies=rep("ZOMBIE_ONE", length(named_bam_files)))
    obs <- tryCatch(mg$plot(design=designTemp), error=conditionMessage)
    exp <- paste0("All design column, except the first one, must be in ", 
                    "numeric format")
    msg <- paste0(base_msg, "A design data.frame with characters in second column ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 

# Design data.frame with invalid second column
test.metagene_plot_invalid_design_not_defined_file <- function() {
    designNew<-data.frame(a=c(named_bam_files, "I am not a file"), 
                          b=rep(1, length(named_bam_files) + 1))
    obs <- tryCatch(mg$plot(design=designNew), 
                        error=conditionMessage)
    exp <- "At least one BAM file does not exist"
    msg <- paste0(base_msg, "A design data.frame with not existing file in ", 
            "first column did not generate an exception with expected message.")
    checkIdentical(obs, exp, msg)
} 

# Design using zero file (0 in all rows of the design object)
test.metagene_plot_design_using_no_file <- function() {
    designNew<-data.frame(a=named_bam_files, 
                          b=rep(0, length(named_bam_files)))
    obs <- tryCatch(mg$plot(design=designNew), 
                    error=conditionMessage)
    exp <- "At least one BAM file must be used in the design"
    msg <- paste0(base_msg, "A design data.frame which does not use BAM file ", 
                  "did not generate an exception with expected message.")
    checkIdentical(obs, exp, msg)
} 

# Invalid regions_group object
test.metagene_plot_invalid_region_group <- function() {
    obs <- tryCatch(mg$plot(regions_group=array(NA, dim = c(2,2,2))), 
                    error=conditionMessage)
    exp <- "regions_group should be a list or a vector"
    msg <- paste0(base_msg, "A invalid regions_group object ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 

# Invalid regions_group object
test.metagene_plot_region <- function() {
    obs <- tryCatch(mg$plot(regions_group=c(regions, "Hello Word!")), 
                    error=conditionMessage)
    exp <- paste0("All elements in regions_group should be regions ",
                  "defined during the creation of metagene object")
    msg <- paste0(base_msg, "A invalid regions_group object ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 

###################################################
## Test the metagene$heatmap() function 
###################################################

base_msg <- "metagene heatmap - "

# Invalid negative bin_size
test.metagene_heatmap_negative_bin_size <- function() {
    obs <- tryCatch(mg$heatmap(bin_size=-2), 
                    error=conditionMessage)
    exp <- "bin_size must be a positive integer"
    msg <- paste0(base_msg, "A negative bin_size ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 

# Invalid zero bin_size
test.metagene_heatmap_negative_bin_size <- function() {
    obs <- tryCatch(mg$heatmap(bin_size=0), 
                    error=conditionMessage)
    exp <- "bin_size must be a positive integer"
    msg <- paste0(base_msg, "A negative bin_size ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 

# Invalid numerical bin_size
test.metagene_heatmap_negative_bin_size <- function() {
    obs <- tryCatch(mg$heatmap(bin_size=2.3), 
                    error=conditionMessage)
    exp <- "bin_size must be a positive integer"
    msg <- paste0(base_msg, "A negative bin_size ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
} 


