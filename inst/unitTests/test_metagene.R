## Test functions present in the metagene.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
}

### }}}

bam_files <- get_demo_bam_files()
named_bam_files <- bam_files
names(named_bam_files) <- letters[1:(length(named_bam_files))]
not_indexed_bam_file <- metagene:::get_not_indexed_bam_file()
regions <- metagene:::get_demo_regions()
design <- data.frame(Samples = c("align1_rep1.bam", "align1_rep2.bam",
                     "align2_rep1.bam", "align2_rep2.bam", "ctrl.bam"),
                     align1 = c(1,1,0,0,2), align2 = c(0,0,1,1,2))
design$Samples <- paste0(system.file("extdata", package="metagene"), "/",
                         design$Samples)
regions_strand <- lapply(regions, rtracklayer::import)
stopifnot(length(unique(vapply(regions, length, numeric(1)))) == 1)
set.seed(1)
index_strand <- sample(1:length(regions_strand[[1]]),
            round(length(regions_strand[[1]])/2))
regions_strand <- lapply(regions_strand,
                         function(x) { strand(x[index_strand]) <- "-"; x })
demo_mg <- metagene$new(regions = get_demo_regions(), bam_files = get_demo_bam_files())
region <- regions[1]
bam_file <- bam_files[1]
demo_mg_min <- metagene$new(regions = region, bam_files = bam_file)

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

## Invalid force_seqlevels value
test.metagene_initialize_invalid_force_seqlevels_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(force_seqlevels="ZOMBIES"),
                    error=conditionMessage)
    exp <- "force_seqlevels must be a logicial value (TRUE or FALSE)"
    msg <- paste0(base_msg, "An invalid force_seqlevels value did not generate an ",
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
    obs <- tryCatch(metagene:::metagene$new(bam_files=bam_files,
            region=array(data = NA, dim = c(2,2,2))),
            error=conditionMessage)
    exp <- paste0("regions must be either a vector of BED filenames, a ",
            "GRanges object or a GrangesList object")
    msg <- paste0(base_msg, "A not indexed bam in bam_files value value ",
            "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Valid regions with extra seqlevels
test.metagene_initialize_valid_regions_supplementary_seqlevels <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- c(GenomeInfoDb::seqlevels(region),
				      "extra_seqlevels")
    mg <- tryCatch(metagene$new(regions = region, bam_files = bam_files[1]),
		   error = conditionMessage)
    msg <- paste0(base_msg, "Valid regions with extra seqlevels did not ")
    msg <- paste0(msg, "return a valid metagene object.")
    checkIdentical(class(mg), c("metagene", "R6"), msg = msg)
}

# Invalid Extra seqnames
test.metagene_initialize_invalid_extra_seqnames <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- "extra_seqlevels"
    obs <- tryCatch(metagene$new(regions = region, bam_files = bam_files[1]),
		   error = conditionMessage)
    exp <- "Some seqnames of regions are absent in bam_file header"
    msg <- paste(base_msg, "Invalid regions seqnames did not give the expected")
    msg <- paste(msg, "error message.")
    checkIdentical(obs, exp, msg)
}

# Extra seqnames with force
test.metagene_initialize_one_extra_seqnames_force_seqlevels <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- c(GenomeInfoDb::seqlevels(region),
					 "extra_seqlevels")
    GenomeInfoDb::seqnames(region)[1] <- "extra_seqlevels"
    mg <- tryCatch(metagene$new(regions = region, bam_files = bam_files[1],
				 force_seqlevels = TRUE),
		   error = conditionMessage)
    msg <- paste(base_msg, "Supplementary seqnames should not have raised an")
    msg <- paste(msg, "error with force_seqlevels = TRUE.")
    checkIdentical(class(mg), c("metagene", "R6"), msg = msg)
}

# Invalid all extra seqnames with force
test.metagene_initialize_all_extra_seqnames_force_seqlevels <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- "extra_seqlevels"
    obs <- tryCatch(metagene$new(regions = region, bam_files = bam_files[1],
				 force_seqlevels = TRUE),
		   error = conditionMessage)
    exp <- "Parameter regions must not be an empty GRanges object"
    msg <- paste(base_msg, "Invalid all extra seqnames did not generate the")
    msg <- paste(msg, "expected error with force_seqlevels = TRUE.")
    checkIdentical(obs, exp, msg)
}

# Valid regions narrowPeak
test.metagene_initialize_valid_narrowpeak <- function() {
    region <- metagene:::get_narrowpeak_region()
    mg <- metagene$new(regions = region, bam_files = bam_files[1])
    obs <- mg$get_regions()$list1
    extraCols <- c(signalValue = "numeric", pValue = "numeric",
		   qValue = "numeric", peak = "integer")
    exp <- rtracklayer::import(region, format = "BED", extraCols = extraCols)
    msg <- paste(base_msg, "Valid narrowPeak file did not generate the")
    mgs <- paste(msg, "expected regions.")
    checkIdentical(obs, exp, msg)
}

# Valid regions broadPeak
test.metagene_initialize_valid_broadpeak <- function() {
    region <- metagene:::get_broadpeak_region()
    mg <- metagene$new(regions = region, bam_files = bam_files[1])
    obs <- mg$get_regions()$list1
    extraCols <- c(signalValue = "numeric", pValue = "numeric",
		   qValue = "numeric")
    exp <- rtracklayer::import(region, format = "BED", extraCols = extraCols)
    msg <- paste(base_msg, "Valid broadPeak file did not generate the")
    mgs <- paste(msg, "expected regions.")
    checkIdentical(obs, exp, msg)
}

# Valid named bam files
test.metagene_initialize_valid_named_bam_files <- function() {
    mg <- metagene$new(regions = regions[1], bam_files = named_bam_files[1])
    obs <- mg$get_params()[["bam_files"]]
    exp <- named_bam_files[1]
    msg <- paste(base_msg, "Valid named bam files did not generate the")
    msg <- paste(msg, "expected bam_files")
    checkIdentical(obs, exp, msg)
    obs <- names(mg$get_raw_coverages())
    exp <- names(named_bam_files)[1]
    msg <- paste(base_msg, "Valid named bam files did not generate the")
    msg <- paste(msg, "expected coverages")
    checkIdentical(obs, exp, msg)
}

# Valid unnamed bam files
test.metagene_initialize_valid_unnamed_bam_files <- function() {
    mg <- metagene$new(regions = regions[1], bam_files = bam_files[1])
    obs <- mg$get_params()[["bam_files"]]
    exp <- bam_files[1]
    names(exp) <- tools::file_path_sans_ext(basename(bam_files[1]))
    msg <- paste(base_msg, "Valid unnamed bam files did not generate the")
    mgs <- paste(msg, "expected results")
    checkIdentical(obs, exp, msg)
    obs <- names(mg$get_raw_coverages())
    exp <- tools::file_path_sans_ext(basename(bam_files[1]))
    msg <- paste(base_msg, "Valid named bam files did not generate the")
    mgs <- paste(msg, "expected coverages")
    checkIdentical(obs, exp, msg)
}

###################################################
## Test the metagene$plot() function 
###################################################

## Valid default
test.metagene_plot_default <- function() {
  mg <- demo_mg$clone()
  mg$produce_data_frame(sample_count = 10)
  pdf(NULL)
  mg$plot()
  dev.off()
  plot <- mg$get_plot()
  checkTrue(all(class(plot) == c("gg", "ggplot")))
}

## Valid show_friedman false
test.metagene_plot_valid_show_friedman_false <- function() {
  mg <- demo_mg_min$clone()
  mg$produce_data_frame(sample_count = 10)
  pdf(NULL)
  mg$plot(show_friedman = FALSE)
  dev.off()
  plot <- mg$get_plot()
  checkTrue(all(class(plot) == c("gg", "ggplot")))
}

## Valid show_friedman true
test.metagene_plot_valid_show_friedman_true <- function() {
  mg <- demo_mg_min$clone()
  mg$produce_data_frame(sample_count = 10)
  pdf(NULL)
  mg$plot(show_friedman = TRUE)
  dev.off()
  plot <- mg$get_plot()
  checkTrue(all(class(plot) == c("gg", "ggplot")))
}

## Invalid show_friedman class
test.metagene_plot_invalid_show_friedman_class <- function() {
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$plot(show_friedman = 1),
                  error=conditionMessage)
  exp <- "is.logical(show_friedman) is not TRUE"
  checkIdentical(obs, exp)
}

## Invalid show_friedman length
test.metagene_plot_invalid_show_friedman_length <- function() {
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$plot(show_friedman = c(TRUE, FALSE)),
                  error=conditionMessage)
  exp <- "length(show_friedman) == 1 is not TRUE"
  checkIdentical(obs, exp)
}

##################################################
# Test the metagene$get_params() function
##################################################

## Valid usage
test.metagene_get_params_valid_usage <- function() {
    mg <- demo_mg$clone()
    params <- mg$get_params()
    checkIdentical(unname(params[["bam_files"]]), get_demo_bam_files())
    checkIdentical(params[["padding_size"]], 0)
    checkIdentical(params[["verbose"]], FALSE)
    checkIdentical(params[["force_seqlevels"]], FALSE)
    checkIdentical(params[["flip_regions"]], FALSE)
}

##################################################
# Test the metagene$get_design() function
##################################################

## Valid usage
test.metagene_get_design_valid_usage <- function() {
    mg <- demo_mg$clone()
    mg$add_design(get_demo_design())
    design <- mg$get_design()
    checkIdentical(design, get_demo_design())
}

##################################################
# Test the metagene$get_regions() function
##################################################

## Valid usage default
test.metagene_get_regions_valid_usage_default <- function() {
    mg <- demo_mg$clone()
    regions <- mg$get_regions()
    exp <- get_demo_regions()
    exp <- tools::file_path_sans_ext(basename(exp))
    checkIdentical(names(regions), exp)
}

## Valid usage subset
test.metagene_get_regions_valid_usage_subset <- function() {
    mg <- demo_mg$clone()
    regions <- mg$get_regions(get_demo_regions()[1])
    exp <- get_demo_regions()[1]
    exp <- tools::file_path_sans_ext(basename(exp))
    checkIdentical(names(regions), exp)
}

## Invalid usage region_names class
test.metagene_get_regions_invalid_usage_region_names_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_regions(1), error = conditionMessage)
    exp <- "a character vector argument expected"
    checkIdentical(obs, exp)
}

## Invalid usage region_names empty
test.metagene_get_regions_invalid_usage_region_names_empty <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_regions(""), error = conditionMessage)
    exp <- "all(region_names %in% names(private$regions)) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage region_names absent
test.metagene_get_regions_invalid_usage_region_names_absent <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_regions("not_valid_name"), error = conditionMessage)
    exp <- "all(region_names %in% names(private$regions)) is not TRUE"
    checkIdentical(obs, exp)
}


##################################################
# Test the metagene$get_matrices() function
##################################################

## Valid usage default
test.metagene_get_matrices_valid_usage_default <- function() {
    mg <- demo_mg$clone()
    mg$produce_matrices()
    matrices <- mg$get_matrices()
    exp_regions <- tools::file_path_sans_ext(basename(get_demo_regions()))
    checkIdentical(names(matrices), exp_regions)
    exp_bam <- tools::file_path_sans_ext(basename(get_demo_bam_files()))
    checkIdentical(names(matrices[[1]]), exp_bam)
    checkTrue(is.matrix(matrices[[1]][[1]][["input"]]))
    checkTrue(length(matrices) == 2)
    checkTrue(length(matrices[[1]]) == 5)
}

## Valid usage subset
test.metagene_get_matrices_valid_usage_subset <- function() {
    mg <- demo_mg$clone()
    mg$produce_matrices()
    matrices <- mg$get_matrices(get_demo_regions()[1],
                                get_demo_bam_files()[1:2])
    exp_regions <- tools::file_path_sans_ext(basename(get_demo_regions()[1]))
    checkIdentical(names(matrices), exp_regions)
    exp_bam <- tools::file_path_sans_ext(basename(get_demo_bam_files()[1:2]))
    checkIdentical(names(matrices[[1]]), exp_bam)
    checkTrue(is.matrix(matrices[[1]][[1]][["input"]]))
    checkTrue(length(matrices) == 1)
    checkTrue(length(matrices[[1]]) == 2)
}

## Valid usage no matrices
test.metagene_get_matrices_valid_usage_no_matrices <- function() {
    mg <- demo_mg$clone()
    matrices <- mg$get_matrices()
    checkTrue(is.null(matrices))
    matrices_subset <- mg$get_matrices(get_demo_regions()[1],
                                get_demo_bam_files()[1:2])
    checkTrue(is.null(matrices_subset))
}

## Invalid usage region_names class
test.metagene_get_matrices_invalid_usage_region_names_class <- function() {
    mg <- demo_mg$clone()$produce_matrices()
    obs <- tryCatch(mg$get_matrices(region_names = 1),
                    error = conditionMessage)
    exp <- "is.character(region_names) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage region_names empty
test.metagene_get_matrices_invalid_usage_region_names_empty <- function() {
    mg <- demo_mg$clone()$produce_matrices()
    obs <- tryCatch(mg$get_matrices(region_names = ""),
                    error = conditionMessage)
    exp <- "all(region_names %in% names(private$matrices)) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage region_names absent
test.metagene_get_matrices_invalid_usage_region_names_absent <- function() {
    mg <- demo_mg$clone()$produce_matrices()
    obs <- tryCatch(mg$get_matrices(region_names = "not_valid_name"),
                    error = conditionMessage)
    exp <- "all(region_names %in% names(private$matrices)) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage exp_names class
test.metagene_get_matrices_invalid_usage_exp_names_class <- function() {
    mg <- demo_mg$clone()$produce_matrices()
    obs <- tryCatch(mg$get_matrices(exp_names = 1),
                    error = conditionMessage)
    exp <- "is.character(exp_names) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage exp_names empty
test.metagene_get_matrices_invalid_usage_exp_names_empty <- function() {
    mg <- demo_mg$clone()$produce_matrices()
    obs <- tryCatch(mg$get_matrices(exp_names = ""),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage exp_names absent
test.metagene_get_matrices_invalid_usage_exp_names_absent <- function() {
    mg <- demo_mg$clone()$produce_matrices()
    obs <- tryCatch(mg$get_matrices(exp_names = "not_valid_name"),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

##################################################
# Test the metagene$get_data_frame() function
##################################################

## Valid usage default
test.metagene_get_data_frame_valid_usage_default <- function() {
    mg <- demo_mg$clone()
    mg$produce_matrices()$produce_data_frame(sample_count = 10)
    df <- mg$get_data_frame()
    regions <- get_demo_regions()
    bam_files <- get_demo_bam_files()
    checkTrue(is.data.frame(df))
    checkTrue(ncol(df) == 5)
    checkTrue(nrow(df) == length(regions) * length(bam_files) * 100)
}

## Valid usage subset
test.metagene_get_data_frame_valid_usage_subset <- function() {
    regions <- get_demo_regions()[1]
    bam_files <- get_demo_bam_files()[1:2]
    mg <- demo_mg$clone()
    mg$produce_matrices()$produce_data_frame(sample_count = 10)
    df <- mg$get_data_frame(region_names = regions, exp_names = bam_files)
    checkTrue(is.data.frame(df))
    checkTrue(ncol(df) == 5)
    checkTrue(nrow(df) == length(regions) * length(bam_files) * 100)
}

## Valid usage no matrices
test.metagene_get_data_frame_valid_usage_no_matrices <- function() {
    mg <- demo_mg$clone()
    df <- mg$get_data_frame()
    checkTrue(is.null(df))
    df_subset <- mg$get_data_frame(get_demo_regions()[1],
                                get_demo_bam_files()[1:2])
    checkTrue(is.null(df_subset))
}

## Invalid usage region_names class
test.metagene_get_data_frame_invalid_usage_region_names_class <- function() {
    mg <- demo_mg$clone()
    mg <- mg$produce_matrices()$produce_data_frame(sample_count = 10)
    obs <- tryCatch(mg$get_data_frame(region_names = 1),
                    error = conditionMessage)
    exp <- "is.character(region_names) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage region_names empty
test.metagene_get_data_frame_invalid_usage_region_names_empty <- function() {
    mg <- demo_mg$clone()
    mg <- mg$produce_matrices()$produce_data_frame(sample_count = 10)
    obs <- tryCatch(mg$get_data_frame(region_names = ""),
                    error = conditionMessage)
    exp <- "all(region_names %in% names(private$matrices)) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage region_names absent
test.metagene_get_data_frame_invalid_usage_region_names_absent <- function() {
    mg <- demo_mg$clone()
    mg <- mg$produce_matrices()$produce_data_frame(sample_count = 10)
    obs <- tryCatch(mg$get_data_frame(region_names = "not_valid_name"),
                    error = conditionMessage)
    exp <- "all(region_names %in% names(private$matrices)) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage exp_names class
test.metagene_get_data_frame_invalid_usage_exp_names_class <- function() {
    mg <- demo_mg$clone()
    mg <- mg$produce_matrices()$produce_data_frame(sample_count = 10)
    obs <- tryCatch(mg$get_data_frame(exp_names = 1),
                    error = conditionMessage)
    exp <- "is.character(exp_names) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage exp_names empty
test.metagene_get_data_frame_invalid_usage_exp_names_empty <- function() {
    mg <- demo_mg$clone()
    mg <- mg$produce_matrices()$produce_data_frame(sample_count = 10)
    obs <- tryCatch(mg$get_data_frame(exp_names = ""),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid usage exp_names absent
test.metagene_get_data_frame_invalid_usage_exp_names_absent <- function() {
    mg <- demo_mg$clone()
    mg <- mg$produce_matrices()$produce_data_frame(sample_count = 10)
    obs <- tryCatch(mg$get_data_frame(exp_names = "not_valid_name"),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

##################################################
# Test the metagene$get_plot() function
##################################################

## Valid case no graph
test.metagene_get_plot_valid_case_no_graph <- function() {
    mg <- demo_mg$clone()
    plot <- mg$get_plot()
    checkTrue(is.null(plot))
}

## Valid case graph
test.metagene_get_plot_valid_case_graph <- function() {
    pdf(NULL)
    mg <- demo_mg$clone()
    mg$produce_data_frame(sample_count = 10)$plot()
    plot <- mg$get_plot()
    dev.off()
    checkTrue(all(class(plot) == c("gg", "ggplot")))
}

##################################################
# Test the metagene$get_raw_coverages() function
##################################################

base_msg <- "metagene get_raw_coverages - "

exp_raw <- GenomicAlignments::readGAlignments(bam_files[1])
exp_raw <- GenomicAlignments::coverage(exp_raw)

## Default filenames
test.metagene_get_raw_coverages_default_filenames <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_raw_coverages()[[1]]
    checkTrue(all(vapply(1:length(obs), function(i)
			 identical(obs[[i]], exp_raw[[i]]), logical(1))))
}

## NULL filenames
test.metagene_get_raw_coverages_null_filenames <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_raw_coverages(filenames = NULL)[[1]]
    checkTrue(all(vapply(1:length(obs), function(i)
			 identical(obs[[i]], exp_raw[[i]]), logical(1))))
}

## One filename
test.metagene_get_raw_coverages_one_filename <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_raw_coverages(filenames = bam_files[1])[[1]]
    checkTrue(all(vapply(1:length(obs), function(i)
			 identical(obs[[i]], exp_raw[[i]]), logical(1))))
}

## All filenames
test.metagene_get_raw_coverages_all_filename <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_raw_coverages(filenames = bam_files)[[1]]
    checkTrue(all(vapply(1:length(obs), function(i)
			 identical(obs[[i]], exp_raw[[i]]), logical(1))))
}

## Invalid filenames class
test.metagene_get_raw_coverages_invalid_filenames_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_raw_coverages(filenames = 1), error = conditionMessage)
    exp <- "is.character(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid empty filename
test.metagene_get_raw_coverages_invalid_empty_filename <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_raw_coverages(filenames = ""), error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename alone
test.metagene_get_raw_coverages_invalid_filename_alone <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_raw_coverages(filenames = "asdf"), error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename among valid
test.metagene_get_raw_coverages_invalid_filename_among_valid <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_raw_coverages(filenames = c("asdf", bam_files)), error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

##################################################
# Test the metagene$get_normalized_coverages() function
##################################################

base_msg <- "metagene get_normalized_coverages - "

count <- Rsamtools::countBam(bam_files[1])$records
weight <- 1 / (count / 1000000)
exp_norm <- exp_raw * weight

## Default filenames
test.metagene_get_normalized_coverages_default_filenames <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_normalized_coverages()[[1]]
    checkTrue(all(vapply(1:length(obs), function(i)
			 identical(obs[[i]], exp_norm[[i]]), logical(1))))
}

## NULL filenames
test.metagene_get_normalized_coverages_null_filenames <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_normalized_coverages(filenames = NULL)[[1]]
    checkTrue(all(vapply(1:length(obs), function(i)
			 identical(obs[[i]], exp_norm[[i]]), logical(1))))
}

## One filename
test.metagene_get_normalized_coverages_one_filename <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_normalized_coverages(filenames = bam_files[1])[[1]]
    checkTrue(all(vapply(1:length(obs), function(i)
			 identical(obs[[i]], exp_norm[[i]]), logical(1))))
}

## All filenames
test.metagene_get_normalized_coverages_all_filename <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_normalized_coverages(filenames = bam_files)[[1]]
    checkTrue(all(vapply(1:length(obs), function(i)
			 identical(obs[[i]], exp_norm[[i]]), logical(1))))
}

## Invalid filenames class
test.metagene_get_normalized_coverages_invalid_filenames_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_normalized_coverages(filenames = 1), error = conditionMessage)
    exp <- "is.character(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid empty filename
test.metagene_get_normalized_coverages_invalid_empty_filename <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_normalized_coverages(filenames = ""), error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename alone
test.metagene_get_normalized_coverages_invalid_filename_alone <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_normalized_coverages(filenames = "asdf"), error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename among valid
test.metagene_get_normalized_coverages_invalid_filename_among_valid <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_normalized_coverages(filenames = c("asdf", bam_files)), error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

##################################################
# Test the metagene$add_design() function
##################################################

## Valid design data frame
test.metagene_add_design_valid_design_data_frame <- function() {
    mg <- demo_mg$clone()
    mg$add_design(get_demo_design())
    checkIdentical(mg$get_design(), get_demo_design())
}

## Valid design NULL
test.metagene_add_design_valid_design_null <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NULL)
    checkIdentical(colnames(mg$get_design())[-1], names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, NA first
test.metagene_add_design_valid_design_na_na_first <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NA)
    checkIdentical(colnames(mg$get_design())[-1], names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, NULL first
test.metagene_add_design_valid_design_na_null_first <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NULL)
    mg$add_design(design = NA)
    checkIdentical(colnames(mg$get_design())[-1], names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, design first
test.metagene_add_design_valid_design_na_design_first <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = get_demo_design())
    checkIdentical(mg$get_design(), get_demo_design())
}

## Valid check_bam_files TRUE NA design
test.metagene_add_design_valid_check_bam_files_true_na_design <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NA, check_bam_files = TRUE)
    checkIdentical(colnames(mg$get_design())[-1], names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid check_bam_files TRUE NULL design
test.metagene_add_design_valid_check_bam_files_true_null_design <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NA, check_bam_files = TRUE)
    checkIdentical(colnames(mg$get_design())[-1], names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid check_bam_files TRUE design design
test.metagene_add_design_valid_check_bam_files_true_design_design <- function()
{
    mg <- demo_mg$clone()
    mg$add_design(design = get_demo_design(), check_bam_files = TRUE)
    checkIdentical(mg$get_design(), get_demo_design())
}

## Invalid design class
test.metagene_add_design_invalid_design_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$add_design(design = 1), error=conditionMessage)
    exp <- "design must be a data.frame object, NULL or NA"
    checkIdentical(obs, exp)
}

## Invalid design column
test.metagene_add_design_invalid_design_column <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design <- design[, 1, drop = FALSE]
    obs <- tryCatch(mg$add_design(design = design), error=conditionMessage)
    exp <- "design must have at least 2 columns"
    checkIdentical(obs, exp)
}

## Invalid design column one class
test.metagene_add_design_invalid_design_column_one_class <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design[,1] <- seq_along(design[,1])
    obs <- tryCatch(mg$add_design(design = design), error=conditionMessage)
    exp <- "The first column of design must be BAM filenames"
    checkIdentical(obs, exp)
}

## Invalid design column two plus class
test.metagene_add_design_invalid_design_columns_two_plus_class <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design[,2] <- letters[seq_along(design[,2])]
    obs <- tryCatch(mg$add_design(design = design), error=conditionMessage)
    exp <- "All design column, except the first one, must be in numeric format"
    checkIdentical(obs, exp)
}

## Invalid check_bam_files class
test.metagene_add_design_invalid_check_bam_files_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$add_design(check_bam_files = 1), error=conditionMessage)
    exp <- "is.logical(check_bam_files) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid bam file check_bam_files TRUE
test.metagene_add_design_invalid_bam_file_check_bam_files_true <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design[1,1] <- "not_a_valid_bam_file"
    obs <- tryCatch(mg$add_design(design = design, check_bam_files = TRUE),
                    error=conditionMessage)
    exp <- "Design contains bam files absent from metagene."
    checkIdentical(obs, exp)
}

## Invalid bam file check_bam_files FALSE
test.metagene_add_design_invalid_bam_file_check_bam_files_false <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design[1,1] <- "not_a_valid_bam_file"
    obs <- mg$add_design(design = design, check_bam_files = FALSE)
    checkIdentical(design, mg$get_design())
}

##################################################
# Test the metagene$produce_matrices() function 
##################################################
#base_msg <- "metagene produce_matrices - "

test.metagene_produce_matrices_valid_default <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$get_params(), FALSE)
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    mg$produce_matrices()
    checkIdentical(mg$get_params()[["bin_size"]], NULL)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(is.list(mg$get_matrices()), TRUE)
    checkIdentical(length(mg$get_matrices()) == 2, TRUE)
    checkIdentical(all(sapply(mg$get_matrices(), class) == c("list", "list")), TRUE)
    checkIdentical(all(sapply(mg$get_matrices(), length) == c(5,5)), TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[5]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[5]][[1]]) == c(50,100)), TRUE)
}

test.metagene_produce_matrices_valid_design <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$get_params(), FALSE)
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    mg$produce_matrices(design = design)
    checkIdentical(mg$get_params()[["bin_size"]], NULL)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(is.list(mg$get_matrices()), TRUE)
    checkIdentical(length(mg$get_matrices()) == 2, TRUE)
    checkIdentical(all(sapply(mg$get_matrices(), class) == c("list", "list")), TRUE)
    checkIdentical(all(sapply(mg$get_matrices(), length) == c(2,2)), TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[2]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[2]][[1]]), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[2]][[1]]) == c(50,100)), TRUE)
}


test.metagene_produce_matrices_valid_bin_count <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$get_params(), FALSE)
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    mg$produce_matrices(bin_count = 200)
    checkIdentical(mg$get_params()[["bin_size"]], NULL)
    checkIdentical(mg$get_params()[["bin_count"]], 200)
    checkIdentical(length(mg$get_matrices()[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[1]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[2]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[3]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[4]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[5]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[1]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[2]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[3]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[4]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[5]][[1]]) == c(50,200)), TRUE)
}

test.metagene_produce_matrices_valid_bin_size <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$get_params(), FALSE)
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    mg$produce_matrices(bin_size = 10)
    checkIdentical(mg$get_params()[["bin_size"]], 10)
    checkIdentical(mg$get_params()[["bin_count"]], 200)
    checkIdentical(length(mg$get_matrices()[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[1]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[2]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[3]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[4]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[5]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[1]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[2]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[3]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[4]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[5]][[1]]) == c(50,200)), TRUE)
}

# Not valid design object
test.metagene_produce_matrices_invalid_design <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$produce_matrices(design=c(1,2)), 
                    error=conditionMessage)
    exp <- "design must be a data.frame object, NULL or NA"
    msg <- paste0(base_msg, "A vector design object ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Design data.frame with not enough columns
test.metagene_produce_matrices_invalid_design_data_frame <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$produce_matrices(design=data.frame(a=c("ZOMBIE_ONE", "ZOMBIE_TWO"))), 
                    error=conditionMessage)
    exp <- "design must have at least 2 columns"
    msg <- paste0(base_msg, "A design data.frame with only one column ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Design data.frame with invalid first column
test.metagene_produce_matrices_invalid_design_first_column <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$produce_matrices(design=data.frame(a=c(1,3), 
                zombies=c("ZOMBIE_ONE", "ZOMBIE_TWO"))), 
                error=conditionMessage)
    exp <- "The first column of design must be BAM filenames"
    msg <- paste0(base_msg, "A design data.frame with numbers in first column ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Design data.frame with invalid second column
test.metagene_produce_matrices_invalid_design_second_column <- function() {
    mg <- demo_mg$clone()
    designTemp<-data.frame(a=named_bam_files, 
                           zombies=rep("ZOMBIE_ONE", length(named_bam_files)))
    obs <- tryCatch(mg$produce_matrices(design=designTemp), error=conditionMessage)
    exp <- paste0("All design column, except the first one, must be in ", 
                    "numeric format")
    msg <- paste0(base_msg, "A design data.frame with characters in second column ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Design data.frame with invalid second column
test.metagene_produce_matrices_invalid_design_not_defined_file <- function() {
    mg <- demo_mg$clone()
    designNew<-data.frame(a=c(bam_files, "I am not a file"),
                          b=rep(1, length(bam_files) + 1))
    obs <- tryCatch(mg$produce_matrices(design=designNew), 
                        error=conditionMessage)
    exp <- "At least one BAM file does not exist"
    msg <- paste0(base_msg, "A design data.frame with not existing file in ", 
            "first column did not generate an exception with expected message.")
    checkIdentical(obs, exp, msg)
}

# Design using zero file (0 in all rows of the design object)
test.metagene_produce_matrices_design_using_no_file <- function() {
    mg <- demo_mg$clone()
    designNew<-data.frame(a=bam_files,
                          b=rep(0, length(bam_files)))
    obs <- tryCatch(mg$produce_matrices(design=designNew), 
                    error=conditionMessage)
    exp <- "At least one BAM file must be used in the design"
    msg <- paste0(base_msg, "A design data.frame which does not use BAM file ", 
                  "did not generate an exception with expected message.")
    checkIdentical(obs, exp, msg)
}

# Invalid bin_count class
test.metagene_produce_matrices_invalid_bin_count_class <- function() {
    mg <- demo_mg$clone()
   obs <- tryCatch(mg$produce_matrices(bin_count = "a"), error = conditionMessage)
   exp <- "bin_count must be NULL or a positive integer"
   msg <- paste0(base_msg, "Invalid bin_count class did not generate the ")
   msg <- paste0(msg, "expected error message.")
   checkIdentical(obs, exp, msg)
}

# Invalid bin_count negative value
test.metagene_produce_matrices_invalid_bin_count_negative_value <- function() {
    mg <- demo_mg$clone()
   obs <- tryCatch(mg$produce_matrices(bin_count = -1), error = conditionMessage)
   exp <- "bin_count must be NULL or a positive integer"
   msg <- paste0(base_msg, "Invalid bin_count negative value did not generate ")
   msg <- paste0(msg, "the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Invalid bin_count decimals
test.metagene_produce_matrices_invalid_bin_count_decimals <- function() {
    mg <- demo_mg$clone()
   obs <- tryCatch(mg$produce_matrices(bin_count = 1.2), error = conditionMessage)
   exp <- "bin_count must be NULL or a positive integer"
   msg <- paste0(base_msg, "Invalid bin_count decimals did not generate ")
   msg <- paste0(msg, "the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Invalid bin_size class
test.metagene_produce_matrices_invalid_bin_size_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$produce_matrices(bin_size = "a"), error = conditionMessage)
    exp <- "bin_size must be NULL or a positive integer"
    msg <- paste0(base_msg, "Invalid bin_size class did not generate the ")
    msg <- paste0(msg, "expected error message.")
    checkIdentical(obs, exp ,msg)
}

# Invalid bin_size negative value
test.metagene_produce_matrices_invalid_bin_size_negative_value <- function() {
   mg <- demo_mg$clone()
   obs <- tryCatch(mg$produce_matrices(bin_size = -1), error = conditionMessage)
   exp <- "bin_size must be NULL or a positive integer"
   msg <- paste0(base_msg, "Invalid bin_size negative value did not generate ")
   msg <- paste0(msg, "the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Invalid bin_size decimals
test.metagene_produce_matrices_invalid_bin_size_decimals <- function() {
   mg <- demo_mg$clone()
   obs <- tryCatch(mg$produce_matrices(bin_size = 1.2), error = conditionMessage)
   exp <- "bin_size must be NULL or a positive integer"
   msg <- paste0(base_msg, "Invalid bin_size decimals did not generate ")
   msg <- paste0(msg, "the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Invalid bin_size regions widths
test.metagene_produce_matrices_invalid_bin_size_regions_width <- function() {
   region <- lapply(regions[1:2], rtracklayer::import)
   width(region[[1]]) <- 1000
   mg <- metagene$new(bam_files=bam_files[1], regions=region)
   obs <- tryCatch(mg$produce_matrices(bin_size = 100), error = conditionMessage)
   exp <- "bin_size can only be used if all selected regions have"
   exp <- paste(exp, "same width")
   msg <- paste0(base_msg, "Invalid bin_size regions width did not generate ")
   msg <- paste0(msg, "the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Warning width not multiple of bin_size
test.metagene_produce_matrices_invalid_bin_size_regions_width_not_multiple <- function() {
   mg <- demo_mg$clone()
   bin_size <- 1234
   width <- 2000
   obs <- tryCatch(mg$produce_matrices(bin_size = 1234), warning = conditionMessage)
   exp <- paste0("width (", width, ") is not a multiple of ")
   exp <- paste0(exp, "bin_size (", bin_size, "), last bin ")
   exp <- paste0(exp, "will be removed.")
   msg <- paste0(base_msg, "Invalid bin_size decimals did not generate ")
   msg <- paste0(msg, "the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Invalid noise_rate class
test.metagene_produce_matrices_invalid_noise_removal_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$produce_matrices(noise_removal = 1234),
                    error = conditionMessage)
    exp <- "noise_removal must be NA, NULL, \"NCIS\" or \"RPM\"."
    msg <- paste0(base_msg, "Invalid noise_removal class did not generate ")
    msg <- paste0("the expected error message.")
    checkIdentical(obs, exp, msg)
}

# Invalid noise_rate value
test.metagene_produce_matrices_invalid_noise_removal_value <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$produce_matrices(noise_removal = "CSI"),
                    error = conditionMessage)
    exp <- "noise_removal must be NA, NULL, \"NCIS\" or \"RPM\"."
    msg <- paste0(base_msg, "Invalid noise_removal class did not generate ")
    msg <- paste0("the expected error message.")
    checkIdentical(obs, exp, msg)
}

# Valid noise_removal NCIS
test.metagene_produce_matrices_valid_noise_removal_ncis <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()[,1:2]
    design[,2][2] <- 0
    mg$produce_matrices(noise_removal = "NCIS", design = design)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["noise_removal"]], "NCIS")
    checkIdentical(length(mg$get_matrices()[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[1]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[1]][[1]]), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[1]][[1]]) == c(50,100)), TRUE)
}

# Invalid normalization class
test.metagene_produce_matrices_invalid_normalization_class <- function() {
   mg <- demo_mg$clone()
   obs <- tryCatch(mg$produce_matrices(normalization = 1234),
                   error = conditionMessage)
   exp <- "normalization must be NA, NULL or \"RPM\"."
   msg <- paste0(base_msg, "Invalid normalization class did not generate ")
   msg <- paste0("the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Invalid normalization value
test.metagene_produce_matrices_invalid_normalization_value <- function() {
   mg <- demo_mg$clone()
   obs <- tryCatch(mg$produce_matrices(normalization = "CSI"),
                   error = conditionMessage)
   exp <- "normalization must be NA, NULL or \"RPM\"."
   msg <- paste0(base_msg, "Invalid normalization class did not generate ")
   msg <- paste0("the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Valid normalization RPM
test.metagene_produce_matrices_valid_normalization_rpm <- function() {
    mg <- demo_mg$clone()
    mg$produce_matrices(normalization = "RPM")
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["normalization"]], "RPM")
    checkIdentical(length(mg$get_matrices()[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[5]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[5]][[1]]) == c(50,100)), TRUE)
}

# Invalid flip_regions class
test.metagene_produce_matrices_invalid_flip_regions_class <- function() {
   mg <- demo_mg$clone()
   obs <- tryCatch(mg$produce_matrices(flip_regions = 1234),
                   error = conditionMessage)
   exp <- "flip_regions must be a logical."
   msg <- paste0(base_msg, "Invalid flip_regions class did not generate ")
   msg <- paste0("the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Valid flip_regions true
test.metagene_produce_matrices_valid_flip_regions_true <- function() {
    mg <- demo_mg$clone()
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_matrices(flip_regions = TRUE)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[5]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[5]][[1]]) == c(50,100)), TRUE)
}

# Valid flip_regions false
test.metagene_produce_matrices_valid_flip_regions_false <- function() {
    mg <- demo_mg$clone()
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_matrices(flip_regions = FALSE)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    checkIdentical(length(mg$get_matrices()[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$get_matrices()[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$get_matrices()[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[1]][[5]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$get_matrices()[[2]][[5]][[1]]) == c(50,100)), TRUE)
}

##################################################
# Test the metagene$produce_data_frame() function
##################################################

## Valid stat bootstrap
test.metagene_produce_data_frame_valid_stat_bootstrap <- function() {
  mg <- demo_mg$clone()
  mg$produce_data_frame(sample_count = 10)
  mg$produce_data_frame(stat = "bootstrap")
  df <- mg$get_data_frame()
  regions <- get_demo_regions()
  bam_files <- get_demo_bam_files()
  checkTrue(is.data.frame(df))
  checkTrue(ncol(df) == 5)
  checkTrue(nrow(df) == length(regions) * length(bam_files) * 100)
}

## Valid stat basic
test.metagene_produce_data_frame_valid_stat_basic <- function() {
  mg <- demo_mg$clone()
  mg$produce_data_frame(sample_count = 10)
  mg$produce_data_frame(stat = "bootstrap")
  df <- mg$get_data_frame()
  regions <- get_demo_regions()
  bam_files <- get_demo_bam_files()
  checkTrue(is.data.frame(df))
  checkTrue(ncol(df) == 5)
  checkTrue(nrow(df) == length(regions) * length(bam_files) * 100)
}

## Invalid stat class
test.metagene_produce_data_frame_invalid_stat_class <- function() {
  mg <- demo_mg$clone()
  obs <- tryCatch(mg$produce_data_frame(stat = 1),
                  error=conditionMessage)
  exp <- "is.character(stat) is not TRUE"
  checkIdentical(obs, exp)
}

## Invalid stat empty
test.metagene_produce_data_frame_invalid_stat_empty <- function() {
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$produce_data_frame(stat = ""),
                  error=conditionMessage)
  exp <- "stat %in% c(\"bootstrap\", \"basic\") is not TRUE"
  checkIdentical(obs, exp)
}

## Invalid stat length greater than one
test.metagene_produce_data_frame_invalid_stat_length_greater_than_one <-
    function() {
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$produce_data_frame(stat = c("bootstrap", "basic")),
                  error=conditionMessage)
  exp <- "length(stat) == 1 is not TRUE"
  checkIdentical(obs, exp)
}

## Invalid stat value
test.metagene_produce_data_frame_invalid_stat_value <- function() {
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$produce_data_frame(stat = "invalid_stat_value"),
                  error=conditionMessage)
  exp <- "stat %in% c(\"bootstrap\", \"basic\") is not TRUE"
  checkIdentical(obs, exp)
}

###################################################
## Test the metagene$flip_regions() function
###################################################

# Valid case not previously flipped
test.metagene_flip_regions_not_previously_flipped <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files[1], regions=regions_strand)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_matrices()
    m1 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$flip_regions()
    m2 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
    mg$flip_regions()
    m3 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
    # Compare the matrices
    checkTrue(identical(m1, m2) == FALSE)
    checkTrue(identical(m2, m3) == TRUE)
    i <- index_strand
    checkIdentical(m1[!(i),], m2[!(i),])
    checkIdentical(m1[i,ncol(m1):1], m2[i,])
}

# Valid case previously flipped
test.metagene_flip_regions_previously_flipped <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files[1], regions=regions_strand)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_matrices(flip_regions = TRUE)
    m1 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
    mg$flip_regions()
    m2 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
    # Compare the matrices
    checkTrue(identical(m1, m2) == TRUE)
}


###################################################
## Test the metagene$unflip_regions() function
###################################################

# Valid case not previously flipped
test.metagene_unflip_regions_not_previously_flipped <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files[1], regions=regions_strand)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_matrices()
    m1 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$unflip_regions()
    m2 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    # Compare the matrices
    checkTrue(identical(m1, m2) == TRUE)
}

# Valid case previously flipped
test.metagene_unflip_regions_previously_flipped <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files[1], regions=regions_strand)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_matrices(flip_regions = TRUE)
    m1 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
    mg$unflip_regions()
    m2 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$unflip_regions()
    m3 <- mg$get_matrices()[[1]][[1]][[1]]
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    # Compare the matrices
    checkTrue(identical(m1, m2) == FALSE)
    checkTrue(identical(m2, m3) == TRUE)
    i <- index_strand
    checkIdentical(m1[!(i),], m2[!(i),])
    checkIdentical(m1[i,ncol(m1):1], m2[i,])
}
