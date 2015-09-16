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
    mg <- tryCatch(metagene$new(regions = region, bam_files = bam_files),
		   error = conditionMessage)
    msg <- paste0(base_msg, "Valid regions with extra seqlevels did not ")
    msg <- paste0(msg, "return a valid metagene object.")
    checkIdentical(class(mg), c("metagene", "R6"), msg = msg)
}

# Invalid Extra seqnames
test.metagene_initialize_invalid_extra_seqnames <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- "extra_seqlevels"
    obs <- tryCatch(metagene$new(regions = region, bam_files = bam_files),
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
    mg <- tryCatch(metagene$new(regions = region, bam_files = bam_files,
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
    obs <- tryCatch(metagene$new(regions = region, bam_files = bam_files,
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
    mg <- metagene$new(regions = region, bam_files = get_demo_bam_files()[1])
    obs <- mg$regions$list1
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
    mg <- metagene$new(regions = region, bam_files = get_demo_bam_files()[1])
    obs <- mg$regions$list1
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
    obs <- mg$params[["bam_files"]]
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
    obs <- mg$params[["bam_files"]]
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
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  pdf(NULL)
  res <- mg$plot()
  dev.off()
  msg <- paste0(base_msg, "Valid default did not return the expected class.")
  checkTrue(class(res) == "list", msg)
  msg <- paste0(base_msg, "Valid default did not return the expected content.")
  checkIdentical(names(res), c("DF", "graph"))
}

## Valid show_friedman false
test.metagene_plot_valid_show_friedman_false <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  pdf(NULL)
  res <- mg$plot(show_friedman = FALSE)
  dev.off()
  msg <- paste0(base_msg, "Valid show_friedman false did not return the")
  msg <- paste(msg, "expected class.")
  checkTrue(class(res) == "list", msg)
  msg <- paste0(base_msg, "Valid show_friedman false did not return the")
  msg <- paste(msg, "expected content.")
  checkIdentical(names(res), c("DF", "graph"))
}

## Valid show_friedman true
test.metagene_plot_valid_show_friedman_true <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  pdf(NULL)
  res <- mg$plot(show_friedman = TRUE)
  dev.off()
  msg <- paste0(base_msg, "Valid show_friedman true did not return the")
  msg <- paste(msg, "expected class.")
  checkTrue(class(res) == "list", msg)
  msg <- paste0(base_msg, "Valid show_friedman true did not return the")
  msg <- paste(msg, "expected content.")
  checkIdentical(names(res), c("DF", "graph"))
}

## Invalid show_friedman class
test.metagene_plot_invalid_show_friedman_class <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$plot(show_friedman = 1),
                  error=conditionMessage)
  exp <- "is.logical(show_friedman) is not TRUE"
  msg <- paste0(base_msg, "Invalid show_friedman class did not return the")
  msg <- paste(msg, "expected error")
  checkIdentical(obs, exp, msg)
}

## Invalid show_friedman length
test.metagene_plot_invalid_show_friedman_length <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$plot(show_friedman = c(TRUE, FALSE)),
                  error=conditionMessage)
  exp <- "length(show_friedman) == 1 is not TRUE"
  msg <- paste0(base_msg, "Invalid show_friedman length did not return the")
  msg <- paste(msg, "expected error")
  checkIdentical(obs, exp, msg)
}

## Valid stat bootstrap
test.metagene_plot_valid_stat_bootstrap <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  pdf(NULL)
  res <- mg$plot(stat = "bootstrap")
  dev.off()
  msg <- paste0(base_msg, "Valid stat bootstrap did not return the expected class.")
  checkTrue(class(res) == "list", msg)
  msg <- paste0(base_msg, "Valid stat bootstrap did not return the expected content.")
  checkIdentical(names(res), c("DF", "graph"))
}

## Valid stat basic
test.metagene_plot_valid_stat_basic <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  pdf(NULL)
  res <- mg$plot(stat = "basic")
  dev.off()
  msg <- paste0(base_msg, "Valid stat basic did not return the expected class.")
  checkTrue(class(res) == "list", msg)
  msg <- paste0(base_msg, "Valid stat basic did not return the expected content.")
  checkIdentical(names(res), c("DF", "graph"))
}

## Invalid stat class
test.metagene_plot_invalid_stat_class <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$plot(stat = 1),
                  error=conditionMessage)
  exp <- "stat %in% c(\"bootstrap\", \"basic\") is not TRUE"
  msg <- paste0(base_msg, "Invalid stat class did not return the expected error")
  checkIdentical(obs, exp, msg)
}

## Invalid stat empty
test.metagene_plot_invalid_stat_empty <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$plot(stat = ""),
                  error=conditionMessage)
  exp <- "stat %in% c(\"bootstrap\", \"basic\") is not TRUE"
  msg <- paste0(base_msg, "Invalid stat empty did not return the expected error")
  checkIdentical(obs, exp, msg)
}

## Invalid stat length greater than one
test.metagene_plot_invalid_stat_empty <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$plot(stat = c("bootstrap", "basic")),
                  error=conditionMessage)
  exp <- "length(stat) == 1 is not TRUE"
  msg <- paste0(base_msg, "Invalid stat empty did not return the expected error")
  checkIdentical(obs, exp, msg)
}


## Invalid stat value
test.metagene_plot_invalid_stat_value <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  obs <- tryCatch(mg$plot(stat = "invalid_stat_value"),
                  error=conditionMessage)
  exp <- "stat %in% c(\"bootstrap\", \"basic\") is not TRUE"
  msg <- paste0(base_msg, "Invalid stat value did not return the expected error")
  checkIdentical(obs, exp, msg)
}

## Valid bin_size
test.metagene_plot_valid_bin_size <- function() {
  base_msg <- "metagene plot - "
  mg <- demo_mg_min$clone()
  pdf(NULL)
  res <- mg$plot(bin_size = 100)
  dev.off()
  msg <- paste0(base_msg, "Valid bin_size did not return the expected class.")
  checkTrue(class(res) == "list", msg)
  msg <- paste0(base_msg, "Valid bin_size did not return the expected content.")
  checkIdentical(names(res), c("DF", "graph"))
}

###################################################
## Test the metagene$heatmap() function 
###################################################

base_msg <- "metagene heatmap - "

# Invalid negative bin_size
test.metagene_heatmap_negative_bin_size <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$heatmap(bin_size=-2), 
                    error=conditionMessage)
    exp <- "bin_size must be a positive integer"
    msg <- paste0(base_msg, "A negative bin_size ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Invalid zero bin_size
test.metagene_heatmap_zero_bin_size <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$heatmap(bin_size=0), 
                    error=conditionMessage)
    exp <- "bin_size must be a positive integer"
    msg <- paste0(base_msg, "A zero bin_size ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Invalid numerical bin_size
test.metagene_heatmap_decimal_bin_size <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$heatmap(bin_size=2.3), 
                    error=conditionMessage)
    exp <- "bin_size must be a positive integer"
    msg <- paste0(base_msg, "A decimal bin_size ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
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
# Test the metagene$produce_matrices() function 
##################################################
#base_msg <- "metagene produce_matrices - "

test.metagene_produce_matrices_valid_default <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$params, FALSE)
    checkIdentical("bin_count" %in% mg$params, FALSE)
    mg$produce_matrices()
    checkIdentical(mg$params[["bin_size"]], NULL)
    checkIdentical(mg$params[["bin_count"]], 100)
    checkIdentical(is.list(mg$matrices), TRUE)
    checkIdentical(length(mg$matrices) == 2, TRUE)
    checkIdentical(all(sapply(mg$matrices, class) == c("list", "list")), TRUE)
    checkIdentical(all(sapply(mg$matrices, length) == c(5,5)), TRUE)
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[5]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[5]][[1]]) == c(50,100)), TRUE)
}

test.metagene_produce_matrices_valid_design <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$params, FALSE)
    checkIdentical("bin_count" %in% mg$params, FALSE)
    mg$produce_matrices(design = design)
    checkIdentical(mg$params[["bin_size"]], NULL)
    checkIdentical(mg$params[["bin_count"]], 100)
    checkIdentical(is.list(mg$matrices), TRUE)
    checkIdentical(length(mg$matrices) == 2, TRUE)
    checkIdentical(all(sapply(mg$matrices, class) == c("list", "list")), TRUE)
    checkIdentical(all(sapply(mg$matrices, length) == c(2,2)), TRUE)
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[2]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[2]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[2]][[1]]) == c(50,100)), TRUE)
}

test.metagene_produce_matrices_valid_select_region <- function() {
    select_regions <- "list1"
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$params, FALSE)
    checkIdentical("bin_count" %in% mg$params, FALSE)
    mg$produce_matrices(select_regions = select_regions)
    checkIdentical(mg$params[["bin_size"]], NULL)
    checkIdentical(mg$params[["bin_count"]], 100)
    checkIdentical(is.list(mg$matrices), TRUE)
    checkIdentical(length(mg$matrices) == 1, TRUE)
    checkIdentical(all(sapply(mg$matrices, class) == "list"), TRUE)
    checkIdentical(all(sapply(mg$matrices, length) == 5), TRUE)
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[5]][[1]]) == c(50,100)), TRUE)
}

test.metagene_produce_matrices_valid_bin_count <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$params, FALSE)
    checkIdentical("bin_count" %in% mg$params, FALSE)
    mg$produce_matrices(bin_count = 200)
    checkIdentical(mg$params[["bin_size"]], NULL)
    checkIdentical(mg$params[["bin_count"]], 200)
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[3]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[4]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[5]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[1]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[2]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[3]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[4]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[5]][[1]]) == c(50,200)), TRUE)
}

test.metagene_produce_matrices_valid_bin_size <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_size" %in% mg$params, FALSE)
    checkIdentical("bin_count" %in% mg$params, FALSE)
    mg$produce_matrices(bin_size = 10)
    checkIdentical(mg$params[["bin_size"]], 10)
    checkIdentical(mg$params[["bin_count"]], 200)
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[3]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[4]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[5]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[1]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[2]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[3]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[4]][[1]]) == c(50,200)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[5]][[1]]) == c(50,200)), TRUE)
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

# Invalid select_regions class
test.metagene_produce_matrices_invalid_select_regions_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$produce_matrices(select_regions=array(NA, dim = c(2,2,2))), 
                    error=conditionMessage)
    exp <- "select_regions must be a character vector."

    msg <- paste0(base_msg, "A invalid select_regions object class ",
                  "did not generate an exception with expected message." )
    checkIdentical(obs, exp, msg)
}

# Invalid select_regions content
test.metagene_produce_matrices_invalid_select_regions_content <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$produce_matrices(select_regions=c(regions, "Hello Word!")), 
                    error=conditionMessage)
    exp <- paste0("All elements in select_regions should be regions ",
                  "defined during the creation of metagene object")
    msg <- paste0(base_msg, "A invalid select_regions object content ",
                  "did not generate an exception with expected message." )
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
   mg <- demo_mg$clone()
   region <- lapply(regions[1:2], rtracklayer::import)
   width(region[[1]]) <- 1000
   mg <- metagene$new(bam_files=bam_files, regions=region)
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
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
    obs <- tryCatch(mg$produce_matrices(noise_removal = 1234),
                    error = conditionMessage)
    exp <- "noise_removal must be NA, NULL, \"NCIS\" or \"RPM\"."
    msg <- paste0(base_msg, "Invalid noise_removal class did not generate ")
    msg <- paste0("the expected error message.")
    checkIdentical(obs, exp, msg)
}

# Invalid noise_rate value
test.metagene_produce_matrices_invalid_noise_removal_value <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
    obs <- tryCatch(mg$produce_matrices(noise_removal = "CSI"),
                    error = conditionMessage)
    exp <- "noise_removal must be NA, NULL, \"NCIS\" or \"RPM\"."
    msg <- paste0(base_msg, "Invalid noise_removal class did not generate ")
    msg <- paste0("the expected error message.")
    checkIdentical(obs, exp, msg)
}

# Valid noise_removal NCIS
test.metagene_produce_matrices_valid_noise_removal_ncis <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
    mg$produce_matrices(noise_removal = "NCIS", design = design)
    checkIdentical(mg$params[["bin_count"]], 100)
    checkIdentical(mg$params[["noise_removal"]], "NCIS")
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[2]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[2]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[2]][[1]]) == c(50,100)), TRUE)
}

# Invalid normalization class
test.metagene_produce_matrices_invalid_normalization_class <- function() {
   mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
   obs <- tryCatch(mg$produce_matrices(normalization = 1234),
                   error = conditionMessage)
   exp <- "normalization must be NA, NULL or \"RPM\"."
   msg <- paste0(base_msg, "Invalid normalization class did not generate ")
   msg <- paste0("the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Invalid normalization value
test.metagene_produce_matrices_invalid_normalization_value <- function() {
   mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
   obs <- tryCatch(mg$produce_matrices(normalization = "CSI"),
                   error = conditionMessage)
   exp <- "normalization must be NA, NULL or \"RPM\"."
   msg <- paste0(base_msg, "Invalid normalization class did not generate ")
   msg <- paste0("the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Valid normalization RPM
test.metagene_produce_matrices_valid_normalization_rpm <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
    mg$produce_matrices(normalization = "RPM")
    checkIdentical(mg$params[["bin_count"]], 100)
    checkIdentical(mg$params[["normalization"]], "RPM")
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[5]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[5]][[1]]) == c(50,100)), TRUE)
}

# Invalid flip_regions class
test.metagene_produce_matrices_invalid_flip_regions_class <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
   obs <- tryCatch(mg$produce_matrices(flip_regions = 1234),
                   error = conditionMessage)
   exp <- "flip_regions must be a logical."
   msg <- paste0(base_msg, "Invalid flip_regions class did not generate ")
   msg <- paste0("the expected error message.")
   checkIdentical(obs, exp, msg)
}

# Valid flip_regions true
test.metagene_produce_matrices_valid_flip_regions_true <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$produce_matrices(flip_regions = TRUE)
    checkIdentical(mg$params[["bin_count"]], 100)
    checkIdentical(mg$params[["flip_regions"]], TRUE)
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[5]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[5]][[1]]) == c(50,100)), TRUE)
}

# Valid flip_regions false
test.metagene_produce_matrices_valid_flip_regions_false <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions)
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$produce_matrices(flip_regions = FALSE)
    checkIdentical(mg$params[["bin_count"]], 100)
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    checkIdentical(length(mg$matrices[[1]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[1]][[5]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[1]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[2]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[3]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[4]]) == 1, TRUE)
    checkIdentical(length(mg$matrices[[2]][[5]]) == 1, TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[1]][[5]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[1]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[2]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[3]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[4]][[1]]), TRUE)
    checkIdentical(is.matrix(mg$matrices[[2]][[5]][[1]]), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[1]][[5]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[1]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[2]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[3]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[4]][[1]]) == c(50,100)), TRUE)
    checkIdentical(all(dim(mg$matrices[[2]][[5]][[1]]) == c(50,100)), TRUE)
}


###################################################
## Test the metagene$flip_regions() function
###################################################

# Valid case not previously flipped
test.metagene_flip_regions_not_previously_flipped <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions_strand)
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$produce_matrices()
    m1 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$flip_regions()
    m2 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], TRUE)
    mg$flip_regions()
    m3 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], TRUE)
    # Compare the matrices
    checkTrue(identical(m1, m2) == FALSE)
    checkTrue(identical(m2, m3) == TRUE)
    i <- index_strand
    checkIdentical(m1[!(i),], m2[!(i),])
    checkIdentical(m1[i,ncol(m1):1], m2[i,])
}

# Valid case previously flipped
test.metagene_flip_regions_previously_flipped <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions_strand)
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$produce_matrices(flip_regions = TRUE)
    m1 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], TRUE)
    mg$flip_regions()
    m2 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], TRUE)
    # Compare the matrices
    checkTrue(identical(m1, m2) == TRUE)
}


###################################################
## Test the metagene$unflip_regions() function
###################################################

# Valid case not previously flipped
test.metagene_unflip_regions_not_previously_flipped <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions_strand)
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$produce_matrices()
    m1 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$unflip_regions()
    m2 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    # Compare the matrices
    checkTrue(identical(m1, m2) == TRUE)
}

# Valid case previously flipped
test.metagene_unflip_regions_previously_flipped <- function() {
    mg <- metagene:::metagene$new(bam_files=bam_files, regions=regions_strand)
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$produce_matrices(flip_regions = TRUE)
    m1 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], TRUE)
    mg$unflip_regions()
    m2 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    mg$unflip_regions()
    m3 <- mg$matrices[[1]][[1]][[1]]
    checkIdentical(mg$params[["flip_regions"]], FALSE)
    # Compare the matrices
    checkTrue(identical(m1, m2) == FALSE)
    checkTrue(identical(m2, m3) == TRUE)
    i <- index_strand
    checkIdentical(m1[!(i),], m2[!(i),])
    checkIdentical(m1[i,ncol(m1):1], m2[i,])
}
