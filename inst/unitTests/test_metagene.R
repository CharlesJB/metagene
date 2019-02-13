## Test functions present in the metagene.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
    library( "data.table" )
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
design$Samples <- paste0(system.file("extdata", package = "metagene"), "/",
                        design$Samples)
regions_strand <- lapply(regions, rtracklayer::import)
stopifnot(length(unique(vapply(regions, length, numeric(1)))) == 1)
set.seed(1)
index_strand <- sample(1:length(regions_strand[[1]]),
                    round(length(regions_strand[[1]])/2))
regions_strand <- lapply(regions_strand,
                        function(x) { strand(x[index_strand]) <- "-"; x })
demo_mg <- metagene$new(regions = get_demo_regions(),
                        bam_files = get_demo_bam_files())
region <- regions[1]
bam_file <- bam_files[1]
demo_mg_min <- metagene$new(regions = region, bam_files = bam_file)

###################################################
## Test the metagene$new() function (initialize)
###################################################

## Invalid verbose value
test.metagene_initialize_invalid_verbose_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(verbose = "ZOMBIES"),
                    error = conditionMessage)
    exp <- "verbose must be a logicial value (TRUE or FALSE)"
    checkIdentical(obs, exp)
}

## Invalid force_seqlevels value
test.metagene_initialize_invalid_force_seqlevels_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(force_seqlevels = "ZOMBIES"),
                    error = conditionMessage)
    exp <- "force_seqlevels must be a logicial value (TRUE or FALSE)"
    checkIdentical(obs, exp)
}

## Negative padding_size value
test.metagene_initialize_negative_padding_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(padding_size = -1),
                    error = conditionMessage)
    exp <- "padding_size must be a non-negative integer"
    checkIdentical(obs, exp)
}

## Non-integer padding_size value
test.metagene_initialize_invalid_string_padding_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(padding_size = "NEW_ZOMBIE"),
                    error = conditionMessage)
    exp <- "padding_size must be a non-negative integer"
    checkIdentical(obs, exp)
}

## Numerical padding_size value
test.metagene_initialize_invalid_numerical_padding_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(padding_size = 1.2),
                    error = conditionMessage)
    exp <- "padding_size must be a non-negative integer"
    checkIdentical(obs, exp)
}

## Negative padding_size value
test.metagene_initialize_negative_padding_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(core = -1),
                    error = conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    checkIdentical(obs, exp)
}

## Non-integer core value
test.metagene_initialize_invalid_string_core_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(core = "ZOMBIE2"),
                    error = conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    checkIdentical(obs, exp)
}

## Numerical core value
test.metagene_initialize_invalid_numerical_core_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(core = 1.2),
                    error = conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    checkIdentical(obs, exp)
}

## Zero core value
test.metagene_initialize_invalid_zero_core_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(core = 0),
                    error = conditionMessage)
    exp <- "cores must be a positive numeric or BiocParallelParam instance"
    checkIdentical(obs, exp)
}

## Non-character vector bam_files value
test.metagene_initialize_invalid_num_vector_bam_files_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(bam_files = c(2,4,3)),
                    error = conditionMessage)
    exp <- "bam_files must be a vector of BAM filenames"
    checkIdentical(obs, exp)
}

## Non-vector bam_files value
test.metagene_initialize_invalid_list_bam_files_value <- function() {
    bam_files <- list(a = "ZOMBIE_01.txt", b = "ZOMBIE_02.txt")
    obs <- tryCatch(metagene:::metagene$new(bam_files = bam_files),
                    error = conditionMessage)
    exp <- "bam_files must be a vector of BAM filenames"
    checkIdentical(obs, exp)
}
# Not indexed bam in bam_files value
test.metagene_initialize_invalid_no_index_bam_files_value <- function() {
    obs <- tryCatch(metagene:::metagene$new(regions = regions, bam_files = not_indexed_bam_file),
                    error = conditionMessage)
    exp <- "All BAM files must be indexed"
    checkIdentical(obs, exp)
}

# Multiple bam files, only one not indexed in bam_files value
test.metagene_initialize_multiple_bam_file_one_not_indexed <- function() {
    bam_files <- c(bam_files, not_indexed_bam_file)
    obs <- tryCatch(metagene:::metagene$new(regions = regions, bam_files = bam_files),
                    error = conditionMessage)
    exp <- "All BAM files must be indexed"
    checkIdentical(obs, exp)
}

# not value for argument region
test.metagene_invalid_initialize_without_region_argument <- function() {
    obs <- tryCatch(metagene:::metagene$new(bam_files = bam_files),
                    error = conditionMessage)
    exp <- 'argument "regions" is missing, with no default'
    checkIdentical(obs, exp)
}

# not value for argument region
test.metagene_invalid_initialize_without_bam_files_argument <- function() {
    obs <- tryCatch(metagene:::metagene$new(regions = regions),
                    error = conditionMessage)
    exp <- 'argument "bam_files" is missing, with no default'
    checkIdentical(obs, exp)
}

# Not valid object in region value
test.metagene_initialize_invalid_array_region_value <- function() {
    region <- array(data = NA, dim = c(2,2,2))
    obs <- tryCatch(metagene:::metagene$new(bam_files = bam_files,
                                            region = region),
                    error = conditionMessage)
    exp <- paste0("regions must be either a vector of BED filenames, a ",
                "GRanges object or a GrangesList object")
    checkIdentical(obs, exp)
}
# Valid regions with extra seqlevels
test.metagene_initialize_valid_regions_supplementary_seqlevels <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- c(GenomeInfoDb::seqlevels(region),
                                        "extra_seqlevels")
    obs <- tryCatch(metagene$new(regions = region, bam_files = bam_files[1]),
                    error = conditionMessage)
    exp <- "Some seqlevels of regions are absent in bam_file"
    checkIdentical(obs, exp)
}

# Valid regions with extra seqlevels force
test.metagene_initialize_valid_regions_supplementary_seqlevels_force <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- c(GenomeInfoDb::seqlevels(region),
                                        "extra_seqlevels")
    obs <- tryCatch(mg <- metagene$new(regions = region, bam_files = bam_files[1],
                force_seqlevels = TRUE),
                    error = conditionMessage)
    checkIdentical(class(mg), c("metagene", "R6"))
}

# Invalid Extra seqnames
test.metagene_initialize_invalid_extra_seqnames <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- "extra_seqlevels"
    obs <- tryCatch(metagene$new(regions = region, bam_files = bam_files[1]),
                    error = conditionMessage)
    exp <- "Some seqlevels of regions are absent in bam_file"
    checkIdentical(obs, exp)
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
    checkIdentical(class(mg), c("metagene", "R6"))
}

# Invalid all extra seqnames with force
test.metagene_initialize_all_extra_seqnames_force_seqlevels <- function() {
    region <- rtracklayer::import(regions[1])
    GenomeInfoDb::seqlevels(region) <- "extra_seqlevels"
    obs <- tryCatch(metagene$new(regions = region, bam_files = bam_files[1],
                                force_seqlevels = TRUE),
                    error = conditionMessage)
    exp <- "No seqlevels matching between regions and bam file"
    checkIdentical(obs, exp)
}

# Valid regions narrowPeak
test.metagene_initialize_valid_narrowpeak <- function() {
    region <- metagene:::get_narrowpeak_region()
    mg <- metagene$new(regions = region, bam_files = bam_files[1])
    obs <- mg$get_regions()$list1
    extraCols <- c(signalValue = "numeric", pValue = "numeric",
                qValue = "numeric", peak = "integer")
    exp <- rtracklayer::import(region, format = "BED", extraCols = extraCols)
    checkIdentical(obs, exp)
}

# Valid regions broadPeak
test.metagene_initialize_valid_broadpeak <- function() {
    region <- metagene:::get_broadpeak_region()
    mg <- metagene$new(regions = region, bam_files = bam_files[1])
    obs <- mg$get_regions()$list1
    extraCols <- c(signalValue = "numeric", pValue = "numeric",
                qValue = "numeric")
    exp <- rtracklayer::import(region, format = "BED", extraCols = extraCols)
    checkIdentical(obs, exp)
}

# Valid named bam files
test.metagene_initialize_valid_named_bam_files <- function() {
    mg <- metagene$new(regions = regions[1], bam_files = named_bam_files[1])
    obs <- mg$get_params()[["bam_files"]]
    exp <- named_bam_files[1]
    checkIdentical(obs, exp)
    obs <- names(mg$get_raw_coverages())
    exp <- names(named_bam_files)[1]
    checkIdentical(obs, exp)
}

# Valid unnamed bam files
test.metagene_initialize_valid_unnamed_bam_files <- function() {
    mg <- metagene$new(regions = regions[1], bam_files = bam_files[1])
    obs <- mg$get_params()[["bam_files"]]
    exp <- bam_files[1]
    names(exp) <- tools::file_path_sans_ext(basename(bam_files[1]))
    checkIdentical(obs, exp)
    obs <- names(mg$get_raw_coverages())
    exp <- tools::file_path_sans_ext(basename(bam_files[1]))
    checkIdentical(obs, exp)
}

###################################################
## Test the metagene$plot() function
###################################################

## Valid default
#test.metagene_plot_default <- function() {
#    mg <- demo_mg$clone()
#    mg$produce_data_frame(sample_count = 10)
#    pdf(NULL)
#    mg$plot()
#    dev.off()
#    plot <- mg$get_plot()
#    checkTrue(all(class(plot) == c("gg", "ggplot")))
#}

## Valid show_friedman false
#test.metagene_plot_valid_show_friedman_false <- function() {
#    mg <- demo_mg_min$clone()
#    mg$produce_data_frame(sample_count = 10)
#    pdf(NULL)
#    mg$plot(show_friedman = FALSE)
#    dev.off()
#    plot <- mg$get_plot()
#    checkTrue(all(class(plot) == c("gg", "ggplot")))
#}

## Valid show_friedman true
#test.metagene_plot_valid_show_friedman_true <- function() {
#    mg <- demo_mg_min$clone()
#    mg$produce_data_frame(sample_count = 10)
#    pdf(NULL)
#    mg$plot(show_friedman = TRUE)
#    dev.off()
#    plot <- mg$get_plot()
#    checkTrue(all(class(plot) == c("gg", "ggplot")))
#}

# ## Invalid show_friedman class
# test.metagene_plot_invalid_show_friedman_class <- function() {
    # mg <- demo_mg_min$clone()
    # obs <- tryCatch(mg$plot(show_friedman = 1),
                    # error = conditionMessage)
    # exp <- "is.logical(show_friedman) is not TRUE"
    # checkIdentical(obs, exp)
# }

# ## Invalid show_friedman length
# test.metagene_plot_invalid_show_friedman_length <- function() {
    # mg <- demo_mg_min$clone()
    # obs <- tryCatch(mg$plot(show_friedman = c(TRUE, FALSE)),
                    # error = conditionMessage)
    # exp <- "length(show_friedman) == 1 is not TRUE"
    # checkIdentical(obs, exp)
# }

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
# Test the metagene$get_table() function
##################################################

## Valid usage default
test.metagene_get_table_valid_usage_default <- function() {
 mg <- demo_mg$clone()
 mg$produce_table()
 tab <- mg$get_table()
 checkTrue(is.data.frame(tab))
 checkTrue(dim(tab)[1] > 0)
 checkTrue(dim(tab)[2] == 5)
 checkIdentical(colnames(tab), c('region', 'design', 'bin', 'value', 'strand'))
}

## Valid usage without producing table before
test.metagene_get_table_without_producing_table_before <- function() {
 mg <- demo_mg$clone()
 tab <- mg$get_table()
 checkIdentical(tab, NULL)
}

## Valid usage get_table return by copy of table
test.metagene_get_table_check_copy_of_table <- function() {
 mg <- demo_mg$clone()
 mg$produce_table()
 tab <- mg$get_table()
 #modification of table by reference
 tab[,c := rep(1:5, length=.N)]
 #Is table copied and unchanged ? 
 tab2 <- mg$get_table()
 checkIdentical(ncol(tab) == ncol(tab2), FALSE)
}

##################################################
# Test the metagene$get_matrice() function
##################################################

test.metagene_get_matrices_valid_usage_default = function(){
    mg <- demo_mg$clone()
    mg$produce_table()
    m <- mg$get_matrices()
    checkIdentical(dim(m$list1$align1_rep1$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list1$align1_rep2$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list1$align2_rep1$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list1$align2_rep2$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list1$ctrl$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list2$align1_rep1$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list2$align1_rep2$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list2$align2_rep1$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list2$align2_rep2$input) == c(50,100), c(TRUE,TRUE))
    checkIdentical(dim(m$list2$ctrl$input) == c(50,100), c(TRUE,TRUE))
}

test.metagene_get_matrices_without_producing_table_before <- function() {
 mg <- demo_mg$clone()
 m <- mg$get_matrices()
 checkIdentical(m, NULL)
}

##################################################
# Test the metagene$get_data_frame() function
##################################################

## Valid usage default
test.metagene_get_data_frame_valid_usage_default <- function() {
 mg <- demo_mg$clone()
 mg$produce_table()$produce_data_frame(sample_count = 10)
 df <- mg$get_data_frame()
 regions <- get_demo_regions()
 bam_files <- get_demo_bam_files()
 checkTrue(is.data.frame(df))
 checkTrue(ncol(df) == 8)
 checkTrue(nrow(df) == length(regions) * length(bam_files) * mg$get_params()$bin_count)
}
#
### Valid usage subset
test.metagene_get_data_frame_valid_usage_subset <- function() {
 regions <- tools::file_path_sans_ext(basename(get_demo_regions()[1]))
 bam_files <- tools::file_path_sans_ext(basename(get_demo_bam_files()[1:2]))
 mg <- demo_mg$clone()
 mg$produce_table()$produce_data_frame(sample_count = 10)
 df <- mg$get_data_frame(region_names = regions, design_names = bam_files)
 checkTrue(is.data.frame(df))
 checkTrue(ncol(df) == 8)
 checkTrue(nrow(df) == length(regions) * length(bam_files) * mg$get_params()$bin_count)
}

#
## Valid usage get_data_frame return by copy of data_frame
test.metagene_get_data_frame_check_copy_of_data_frame <- function() {
 mg <- demo_mg$clone()
 mg$produce_table()
 mg$produce_data_frame()
 df1 <- mg$get_data_frame()
 #modification of table by reference
 df1$c <- 1:1000
 #Is table copied and unchanged ? 
 df2 <- mg$get_data_frame()
 checkIdentical(ncol(df1) == ncol(df2), FALSE)
}

#
## Valid usage no data_frame produced
test.metagene_get_data_frame_valid_usage_no_data_frame <- function() {
 mg <- demo_mg$clone()
 df <- mg$get_data_frame()
 checkTrue(is.null(df))
 df_subset <- mg$get_data_frame(get_demo_regions()[1],
                                get_demo_bam_files()[1:2])
 checkTrue(is.null(df_subset))
}
#
### Invalid usage region_names class
test.metagene_get_data_frame_invalid_usage_region_names_class <- function() {
 mg <- demo_mg$clone()
 mg <- mg$produce_table()$produce_data_frame(sample_count = 10)
 obs <- tryCatch(mg$get_data_frame(region_names = 1),
                error = conditionMessage)
 exp <- "is.character(region_names) is not TRUE"
 checkIdentical(obs, exp)
}
#
### Invalid usage region_names empty
test.metagene_get_data_frame_invalid_usage_region_names_empty <- function() {
 mg <- demo_mg$clone()
 mg <- mg$produce_table()$produce_data_frame(sample_count = 10)
 obs <- tryCatch(mg$get_data_frame(region_names = ""),
                error = conditionMessage)
 exp <- "all(region_names %in% unique(private$table$region)) is not TRUE"
 checkIdentical(obs, exp)
}
#
### Invalid usage region_names absent
test.metagene_get_data_frame_invalid_usage_region_names_absent <- function() {
 mg <- demo_mg$clone()
 mg <- mg$produce_table()$produce_data_frame(sample_count = 10)
 obs <- tryCatch(mg$get_data_frame(region_names = "not_valid_name"),
                error = conditionMessage)
 exp <- "all(region_names %in% unique(private$table$region)) is not TRUE"
 checkIdentical(obs, exp)
}
#
### Valid usage exp_name no design
test.metagene_get_data_frame_valid_usage_design_names_no_design <- function() {
 mg <- demo_mg$clone()$produce_table()
 mg$produce_data_frame(sample_count = 10)
 exp_name <- tools::file_path_sans_ext(basename(get_demo_bam_files()[1]))
 nodesign <- unique(mg$get_data_frame(design_names = exp_name)$design)
 checkIdentical(nodesign, exp_name)
}
#
### Valid usage exp_name design
test.metagene_get_data_frame_valid_usage_design_names_exist_design <- function() {
 mg <- demo_mg$clone()$produce_table(design = get_demo_design())
 mg$produce_data_frame(sample_count = 10)
 exp_name <- unique(mg$get_table()$design)
 yesdesign <- unique(mg$get_data_frame(design_names = exp_name)$design)
 checkIdentical(yesdesign, exp_name)
}
#
### Invalid usage exp_name bam_file design
test.metagene_get_data_frame_invalid_usage_design_names_bam_file_design <-
    function() {
 mg <- demo_mg$clone()$produce_table(design = get_demo_design())
 mg$produce_data_frame(sample_count = 10)
 exp_name <- tools::file_path_sans_ext(basename(get_demo_bam_files()[1]))
 obs <- tryCatch(mg$get_data_frame(design_names = exp_name),
                error = conditionMessage)
 exp <- "all(design_names %in% unique(private$table$design)) is not TRUE"
 checkIdentical(obs, exp)
}
#
### Invalid usage design_names class
test.metagene_get_data_frame_invalid_usage_design_names_class <- function() {
 mg <- demo_mg$clone()
 mg <- mg$produce_table()$produce_data_frame(sample_count = 10)
 obs <- tryCatch(mg$get_data_frame(design_names = 1),
                error = conditionMessage)
 exp <- "is.character(design_names) is not TRUE"
 checkIdentical(obs, exp)
}
#
### Invalid usage design_names empty
test.metagene_get_data_frame_invalid_usage_design_names_empty <- function() {
 mg <- demo_mg$clone()
 mg <- mg$produce_table()$produce_data_frame(sample_count = 10)
 obs <- tryCatch(mg$get_data_frame(design_names = ""),
                error = conditionMessage)
 exp <- "all(design_names %in% unique(private$table$design)) is not TRUE"
 checkIdentical(obs, exp)
}
#
### Invalid usage design_names absent
test.metagene_get_data_frame_invalid_usage_design_names_absent <- function() {
 mg <- demo_mg$clone()
 mg <- mg$produce_table()$produce_data_frame(sample_count = 10)
 obs <- tryCatch(mg$get_data_frame(design_names = "not_valid_name"),
                error = conditionMessage)
 exp <- "all(design_names %in% unique(private$table$design)) is not TRUE"
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
#test.metagene_get_plot_valid_case_graph <- function() {
#    pdf(NULL)
#    mg <- demo_mg$clone()
#    mg$produce_data_frame(sample_count = 10)$plot()
#    plot <- mg$get_plot()
#    dev.off()
#    checkTrue(all(class(plot) == c("gg", "ggplot")))
#}

##################################################
# Test the metagene$get_raw_coverages() function
##################################################

exp_raw <- GenomicAlignments::readGAlignments(bam_files[1])
exp_raw <- GenomicAlignments::coverage(exp_raw)

## Default filenames
test.metagene_get_raw_coverages_default_filenames <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_raw_coverages()[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) identical(obs[[i]], exp_raw[[i]]),
                        logical(1))))
}

## NULL filenames
test.metagene_get_raw_coverages_null_filenames <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_raw_coverages(filenames = NULL)[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) identical(obs[[i]], exp_raw[[i]]),
                        logical(1))))
}

## One filename
test.metagene_get_raw_coverages_one_filename <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_raw_coverages(filenames = bam_files[1])[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) identical(obs[[i]], exp_raw[[i]]),
                        logical(1))))
}

## All filenames
test.metagene_get_raw_coverages_all_filename <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_raw_coverages(filenames = bam_files)[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) identical(obs[[i]], exp_raw[[i]]),
                        logical(1))))
}

## Invalid filenames class
test.metagene_get_raw_coverages_invalid_filenames_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_raw_coverages(filenames = 1),
                    error = conditionMessage)
    exp <- "is.character(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid empty filename
test.metagene_get_raw_coverages_invalid_empty_filename <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_raw_coverages(filenames = ""),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename alone
test.metagene_get_raw_coverages_invalid_filename_alone <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_raw_coverages(filenames = "asdf"),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename among valid
test.metagene_get_raw_coverages_invalid_filename_among_valid <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_raw_coverages(filenames = c("asdf", bam_files)),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

##################################################
# Test the metagene$get_normalized_coverages() function
##################################################

count <- Rsamtools::countBam(bam_files[1])$records
weight <- 1 / (count / 1000000)
exp_norm <- exp_raw * weight

## Default filenames
test.metagene_get_normalized_coverages_default_filenames <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_normalized_coverages()[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) identical(obs[[i]], exp_norm[[i]]),
                        logical(1))))
}

## NULL filenames
test.metagene_get_normalized_coverages_null_filenames <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_normalized_coverages(filenames = NULL)[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) identical(obs[[i]], exp_norm[[i]]),
                        logical(1))))
}

## One filename
test.metagene_get_normalized_coverages_one_filename <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_normalized_coverages(filenames = bam_files[1])[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) identical(obs[[i]], exp_norm[[i]]),
                        logical(1))))
}

## All filenames
test.metagene_get_normalized_coverages_all_filename <- function() {
    mg <- demo_mg$clone()
    obs <- mg$get_normalized_coverages(filenames = bam_files)[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) identical(obs[[i]], exp_norm[[i]]),
                        logical(1))))
}

## Invalid filenames class
test.metagene_get_normalized_coverages_invalid_filenames_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_normalized_coverages(filenames = 1),
                    error = conditionMessage)
    exp <- "is.character(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid empty filename
test.metagene_get_normalized_coverages_invalid_empty_filename <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_normalized_coverages(filenames = ""),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename alone
test.metagene_get_normalized_coverages_invalid_filename_alone <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$get_normalized_coverages(filenames = "asdf"),
                    error = conditionMessage)
    exp <- "private$check_bam_files(filenames) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid filename among valid
test.metagene_get_normalized_coverages_invalid_filename_among_valid <- 
    function() {
    mg <- demo_mg$clone()
    filenames <- c("asdf", bam_files)
    obs <- tryCatch(mg$get_normalized_coverages(filenames = filenames),
                    error = conditionMessage)
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
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, NA first
test.metagene_add_design_valid_design_na_na_first <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NA)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, NULL first
test.metagene_add_design_valid_design_na_null_first <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NULL)
    mg$add_design(design = NA)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid design NA, design first
test.metagene_add_design_valid_design_na_design_first <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = get_demo_design())
    checkIdentical(mg$get_design(), get_demo_design())
}

## Valid design, factor sample names
test.metagene_add_design_valid_design_factor_sample_names <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design[,1] <- factor(design[,1])
    mg$add_design(design)
    checkTrue(is.factor(design[,1]))
    checkTrue(is.character(mg$get_design()[,1]))
    checkIdentical(design[,-1], mg$get_design()[,-1])
    checkIdentical(as.character(design[,1]), mg$get_design()[,1])
}

## Valid check_bam_files TRUE NA design
test.metagene_add_design_valid_check_bam_files_true_na_design <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NA, check_bam_files = TRUE)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
    checkTrue(all(apply(mg$get_design()[,-1], 2, sum) == 1))
}

## Valid check_bam_files TRUE NULL design
test.metagene_add_design_valid_check_bam_files_true_null_design <- function() {
    mg <- demo_mg$clone()
    mg$add_design(design = NA, check_bam_files = TRUE)
    checkIdentical(colnames(mg$get_design())[-1],
                names(mg$get_params()[["bam_files"]]))
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
    obs <- tryCatch(mg$add_design(design = 1), error = conditionMessage)
    exp <- "design must be a data.frame object, NULL or NA"
    checkIdentical(obs, exp)
}

## Invalid design column
test.metagene_add_design_invalid_design_column <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design <- design[, 1, drop = FALSE]
    obs <- tryCatch(mg$add_design(design = design), error = conditionMessage)
    exp <- "design must have at least 2 columns"
    checkIdentical(obs, exp)
}

## Invalid design column one class
test.metagene_add_design_invalid_design_column_one_class <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design[,1] <- seq_along(design[,1])
    obs <- tryCatch(mg$add_design(design = design), error = conditionMessage)
    exp <- "The first column of design must be BAM filenames"
    checkIdentical(obs, exp)
}

## Invalid design column two plus class
test.metagene_add_design_invalid_design_columns_two_plus_class <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design[,2] <- letters[seq_along(design[,2])]
    obs <- tryCatch(mg$add_design(design = design), error = conditionMessage)
    exp <- "All design column, except the first one, must be in numeric format"
    checkIdentical(obs, exp)
}

## Invalid check_bam_files class
test.metagene_add_design_invalid_check_bam_files_class <- function() {
    mg <- demo_mg$clone()
    obs <- tryCatch(mg$add_design(check_bam_files = 1),
                    error = conditionMessage)
    exp <- "is.logical(check_bam_files) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid bam file check_bam_files TRUE
test.metagene_add_design_invalid_bam_file_check_bam_files_true <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()
    design[1,1] <- "not_a_valid_bam_file"
    obs <- tryCatch(mg$add_design(design = design, check_bam_files = TRUE),
                    error = conditionMessage)
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
# Test the metagene$produce_table() function
##################################################

test.metagene_produce_table_valid_without_design <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    mg$produce_table()
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(is.data.frame(mg$get_table()), TRUE)
    #length of table : number of region * number of design * number of bin * number of range by region (demo = 50,000 lines)
    tablength <- length(mg$get_regions())*length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(mg$get_table())[1] == tablength, TRUE)
    checkIdentical(dim(mg$get_table())[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(tools::file_path_sans_ext(basename(mg$get_params()$bam_files)), unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in tools::file_path_sans_ext(basename(mg$get_params()$bam_files))){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}

test.metagene_produce_table_valid_with_design <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    demo_design <- get_demo_design()
    mg$produce_table(design = demo_design)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(is.data.frame(mg$get_table()), TRUE)
    #length of table : number of region * number of design * number of bin * number of range by region (demo = 50,000 lines)
    tablength <- length(mg$get_regions())*(dim(demo_design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(mg$get_table())[1] == tablength, TRUE)
    checkIdentical(dim(mg$get_table())[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(names(demo_design)[-1], unique(tab$design))
    
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- (dim(demo_design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in names(demo_design)[-1]){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}

test.metagene_produce_table_valid_without_design_bin_count_50 <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    mg$produce_table(bin_count = 50)
    checkIdentical(mg$get_params()[["bin_count"]], 50)
    checkIdentical(is.data.frame(mg$get_table()), TRUE)
    #length of table : number of region * number of design * number of bin * number of range by region (demo = 50,000 lines)
    tablength <- length(mg$get_regions())*length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(mg$get_table())[1] == tablength, TRUE)
    checkIdentical(dim(mg$get_table())[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(tools::file_path_sans_ext(basename(mg$get_params()$bam_files)), unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in tools::file_path_sans_ext(basename(mg$get_params()$bam_files))){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}

test.metagene_produce_table_valid_with_design_bin_count_50 <- function() {
    mg <- demo_mg$clone()
    checkIdentical("bin_count" %in% mg$get_params(), FALSE)
    demo_design <- get_demo_design()
    mg$produce_table(design = demo_design, bin_count = 50)
    checkIdentical(mg$get_params()[["bin_count"]], 50)
    checkIdentical(is.data.frame(mg$get_table()), TRUE)
    #length of table : number of region * number of design * number of bin * number of range by region (demo = 50,000 lines)
    tablength <- length(mg$get_regions())*(dim(demo_design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(mg$get_table())[1] == tablength, TRUE)
    checkIdentical(dim(mg$get_table())[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(names(demo_design)[-1], unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- (dim(demo_design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in names(demo_design)[-1]){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}
#
# Not valid design object
test.metagene_produce_table_invalid_design <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(design = c(1,2)),
                error = conditionMessage)
 exp <- "design must be a data.frame object, NULL or NA"
 checkIdentical(obs, exp)
}
#
# Design data.frame with not enough columns
test.metagene_produce_table_invalid_design_data_frame <- function() {
 mg <- demo_mg$clone()
 design <- data.frame(a = c("ZOMBIE_ONE", "ZOMBIE_TWO"))
 obs <- tryCatch(mg$produce_table(design = design),
                error = conditionMessage)
 exp <- "design must have at least 2 columns"
 checkIdentical(obs, exp)
}

# Design data.frame with invalid first column
test.metagene_produce_table_invalid_design_first_column <- function() {
 mg <- demo_mg$clone()
 design <- data.frame(a = c(1,3), zombies = c("ZOMBIE_ONE", "ZOMBIE_TWO"))
 obs <- tryCatch(mg$produce_table(design = design),
                error = conditionMessage)
 exp <- "The first column of design must be BAM filenames"
 checkIdentical(obs, exp)
}

# Design data.frame with invalid second column
test.metagene_produce_table_invalid_design_second_column <- function() {
 mg <- demo_mg$clone()
 designTemp<-data.frame(a = named_bam_files,
                        zombies = rep("ZOMBIE_ONE", length(named_bam_files)))
 obs <- tryCatch(mg$produce_table(design = designTemp),
                error = conditionMessage)
 exp <- paste0("All design column, except the first one, must be in ",
                "numeric format")
 checkIdentical(obs, exp)
}

# Design data.frame with invalid second column
test.metagene_produce_table_invalid_design_not_defined_file <- function() {
 mg <- demo_mg$clone()
 designNew<-data.frame(a = c(bam_files, "I am not a file"),
                        b = rep(1, length(bam_files) + 1))
 obs <- tryCatch(mg$produce_table(design = designNew),
                error = conditionMessage)
 exp <- "At least one BAM file does not exist"
 checkIdentical(obs, exp)
}

# Design using zero file (0 in all rows of the design object)
test.metagene_produce_table_design_using_no_file <- function() {
 mg <- demo_mg$clone()
 designNew<-data.frame(a = bam_files,
                        b = rep(0, length(bam_files)))
 obs <- tryCatch(mg$produce_table(design = designNew),
                error = conditionMessage)
 exp <- "At least one BAM file must be used in the design"
 checkIdentical(obs, exp)
}

# Invalid bin_count class
test.metagene_produce_table_invalid_bin_count_class <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(bin_count = "a"),
                error = conditionMessage)
 exp <- "bin_count must be NULL or a positive integer"
 checkIdentical(obs, exp)
}

# Invalid bin_count negative value
test.metagene_produce_table_invalid_bin_count_negative_value <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(bin_count = -1),
                error = conditionMessage)
 exp <- "bin_count must be NULL or a positive integer"
 checkIdentical(obs, exp)
}

# Invalid bin_count decimals
test.metagene_produce_table_invalid_bin_count_decimals <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(bin_count = 1.2),
                error = conditionMessage)
 exp <- "bin_count must be NULL or a positive integer"
 checkIdentical(obs, exp)
}

# Invalid noise_rate class
test.metagene_produce_table_invalid_noise_removal_class <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(noise_removal = 1234),
                error = conditionMessage)
 exp <- "noise_removal must be NA, NULL, \"NCIS\" or \"RPM\"."
 checkIdentical(obs, exp)
}

# Invalid noise_rate value
test.metagene_produce_table_invalid_noise_removal_value <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(noise_removal = "CSI"),
                error = conditionMessage)
 exp <- "noise_removal must be NA, NULL, \"NCIS\" or \"RPM\"."
 checkIdentical(obs, exp)
}

# Valid noise_removal NCIS
test.metagene_produce_table_valid_noise_removal_ncis <- function() {
    mg <- demo_mg$clone()
    design <- get_demo_design()[,1:2]
    design[,2][2] <- 0
    mg$produce_table(noise_removal = "NCIS", design = design)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["noise_removal"]], "NCIS")
    tab <- mg$get_table()
    tablength <- length(mg$get_regions())*(dim(design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(tab)[1] == tablength, TRUE)
    checkIdentical(dim(tab)[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(names(design)[-1], unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- (dim(design)[2]-1)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in names(design)[-1]){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}

# Invalid normalization class
test.metagene_produce_table_invalid_normalization_class <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(normalization = 1234),
                error = conditionMessage)
 exp <- "normalization must be NA, NULL or \"RPM\"."
 checkIdentical(obs, exp)
}

# Invalid normalization value
test.metagene_produce_table_invalid_normalization_value <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(normalization = "CSI"),
                error = conditionMessage)
 exp <- "normalization must be NA, NULL or \"RPM\"."
 checkIdentical(obs, exp)
}
#
## Valid normalization RPM
test.metagene_produce_table_valid_normalization_rpm <- function() {
    mg <- demo_mg$clone()
    mg$produce_table(normalization = "RPM")
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["normalization"]], "RPM")
    tab <- mg$get_table()
    tablength <- length(mg$get_regions())*length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    checkIdentical(dim(tab)[1] == tablength, TRUE)
    checkIdentical(dim(tab)[2] == length(c('region', 'design', 'bin', 'value', 'strand')), TRUE)
    
    tab <- mg$get_table()
    #check for presence of levels of factors (region, design, strand)
    checkIdentical(names(mg$get_regions()), unique(tab$region))
    checkIdentical(tools::file_path_sans_ext(basename(mg$get_params()$bam_files)), unique(tab$design))
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(tab$strand[which(tab$region == region_names)]))
    }
    #check for number of line by factor region
    reglength <- length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(tab$region[which(tab$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)*length(mg$get_regions()[[1]])
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in tools::file_path_sans_ext(basename(mg$get_params()$bam_files))){
            #print(design_names)
            checkIdentical(length(tab$design[which(tab$region == region_names & tab$design == design_names)]) == designlength , TRUE)
        }
    }
    print(TRUE)
}
#
## Invalid flip_regions class
test.metagene_produce_table_invalid_flip_regions_class <- function() {
 mg <- demo_mg$clone()
 obs <- tryCatch(mg$produce_table(flip_regions = 1234),
                error = conditionMessage)
 exp <- "flip_regions must be a logical."
 checkIdentical(obs, exp)
}
#
## Valid flip_regions true
test.metagene_produce_table_valid_flip_regions_true <- function() {
    mg <- demo_mg$clone()
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_table(flip_regions = TRUE)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
    
    #modifier strand char for regions
    
    
    #test expected == observed
    tab <- mg$get_table()
    # print(tab)
    expect <- c()
    for (region_names in names(mg$get_regions())){
        if(as.vector(strand(mg$get_regions()[[region_names]]))[1] == "-") {
            expect <- c(expect,rep(100:1,length(names(mg$get_design())[-1])*length(as.vector(strand(mg$get_regions()[[region_names]])))))
        } else {
            expect <- c(expect,rep(1:100,length(names(mg$get_design())[-1])*length(as.vector(strand(mg$get_regions()[[region_names]])))))
        }
    }
    # print(class(tab$bin))
    # print(class(expect))
    checkIdentical(as.numeric(tab$bin), as.numeric(expect))
}
#
## Valid flip_regions false
test.metagene_produce_table_valid_flip_regions_false <- function() {
    mg <- demo_mg$clone()
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    mg$produce_table(flip_regions = FALSE)
    checkIdentical(mg$get_params()[["bin_count"]], 100)
    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
    
    #modifier strand char dans regions
    
    
    #test expected == observed
    tab <- mg$get_table()
    # print(tab)
    expect <- c()
    for (region_names in names(mg$get_regions())){
        # if(as.vector(strand(mg$get_regions()[[region_names]]))[1] == "-") {
            # expect <- c(expect,rep(100:1,length(names(mg$get_design())[-1])*length(as.vector(strand(mg$get_regions()[[region_names]])))))
        # } else {
            expect <- c(expect,rep(1:100,length(names(mg$get_design())[-1])*length(as.vector(strand(mg$get_regions()[[region_names]])))))
        # }
    }
    # print(class(tab$bin))
    # print(class(expect))
    checkIdentical(as.numeric(tab$bin), as.numeric(expect))
}

##################################################
# Test the metagene$produce_data_frame() function
##################################################

## Valid default usage
test.metagene_produce_data_frame_default_arguments <- function(){
    mg <- demo_mg$clone()
    mg$produce_table()
    mg$produce_data_frame()
    df <- mg$get_data_frame()
    #check the nrow & ncol
    checkIdentical(ncol(df), length(c('region', 'design', 'bin', 'value', 'strand', 'qinf', 'qsup', 'group')))
    expectedNbRow <- length(mg$get_regions())*length(mg$get_params()$bam_files)*(mg$get_params()$bin_count)
    checkIdentical(nrow(df) == expectedNbRow, TRUE)
    #check colnames
    checkIdentical(colnames(df), c('region', 'design', 'bin', 'value', 'strand', 'qinf', 'qsup','group'))
    #check region, design repartition in data_frame
    checkIdentical(names(mg$get_regions()), unique(df$region))
    checkIdentical(names(mg$get_design()[-1]), unique(df$design))
    #check for number of line by factor region
    reglength <- length(names(mg$get_design())[-1])*(mg$get_params()$bin_count)
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(length(df$region[which(df$region == region_names)]) == reglength , TRUE)
    }
    #check for number of line by factor design
    designlength <- (mg$get_params()$bin_count)
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        for (design_names in names(mg$get_design()[-1])){
            #print(design_names)
            checkIdentical(length(df$design[which(df$region == region_names & df$design == design_names)]) == designlength , TRUE)
        }
    }
    #check for bin repartition
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(sum(df$bin), sum(1:100)*length(names(mg$get_design())[-1])*length(mg$get_regions()))
    }
    #strand matches
    for (region_names in names(mg$get_regions())){
        #print(region_names)
        checkIdentical(unique(as.vector(strand(mg$get_regions())[[region_names]])), unique(df$strand[which(df$region == region_names)]))
    }
    print(TRUE)
}

### Invalid alpha class
test.metagene_produce_data_frame_invalid_alpha_class <- function(){
    mg <- demo_mg$clone()
    mg$produce_table()
    obs <- tryCatch(mg$produce_data_frame(alpha='test'),
                error = conditionMessage)
    exp <- "is.numeric(alpha) is not TRUE"
    checkIdentical(obs, exp)
}
### Invalid alpha value
test.metagene_produce_data_frame_invalid_alpha_value <- function(){
    mg <- demo_mg$clone()
    mg$produce_table()
    obs <- tryCatch(mg$produce_data_frame(alpha=-0.8),
                error = conditionMessage)
    exp <- "alpha >= 0 & alpha <= 1 is not TRUE"
    checkIdentical(obs, exp)
}
### Invalid sample_count class
test.metagene_produce_data_frame_invalid_sample_count_class <- function(){
    mg <- demo_mg$clone()
    mg$produce_table()
    obs <- tryCatch(mg$produce_data_frame(sample_count='test'),
                error = conditionMessage)
    exp <- "is.numeric(sample_count) is not TRUE"
    checkIdentical(obs, exp)
}
### Invalid sample_count value
test.metagene_produce_data_frame_invalid_sample_count_value <- function(){
    mg <- demo_mg$clone()
    mg$produce_table()
    obs <- tryCatch(mg$produce_data_frame(sample_count=0),
                error = conditionMessage)
    exp <- "sample_count > 0 is not TRUE"
    checkIdentical(obs, exp)
    obs <- tryCatch(mg$produce_data_frame(sample_count=-10),
                error = conditionMessage)
    exp <- "sample_count > 0 is not TRUE"
    checkIdentical(obs, exp)
}

###################################################
## Test the metagene$flip_regions() function
###################################################
# TODO: Re-code later
## Valid case not previously flipped
#test.metagene_flip_regions_not_previously_flipped <- function() {
#    mg <- metagene:::metagene$new(bam_files = bam_files[1],
#                                regions = regions_strand)
#    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
#    mg$produce_table()
#    m1 <- mg$get_table()[[1]][[1]][[1]]
#    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
#    mg$flip_regions()
#    m2 <- mg$get_table()[[1]][[1]][[1]]
#    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
#    mg$flip_regions()
#    m3 <- mg$get_table()[[1]][[1]][[1]]
#    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
#    # Compare the matrices
#    checkTrue(identical(m1, m2) == FALSE)
#    checkTrue(identical(m2, m3) == TRUE)
#    i <- index_strand
#    checkIdentical(m1[!(i),], m2[!(i),])
#    checkIdentical(m1[i,ncol(m1):1], m2[i,])
#}
#
## Valid case previously flipped
#test.metagene_flip_regions_previously_flipped <- function() {
#    mg <- metagene:::metagene$new(bam_files = bam_files[1],
#                                regions = regions_strand)
#    checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
#    mg$produce_table(flip_regions = TRUE)
#    m1 <- mg$get_table()[[1]][[1]][[1]]
#    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
#    mg$flip_regions()
#    m2 <- mg$get_table()[[1]][[1]][[1]]
#    checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
#    # Compare the matrices
#    checkTrue(identical(m1, m2) == TRUE)
#}
#
#
####################################################
### Test the metagene$unflip_regions() function
####################################################
#
## Valid case not previously flipped
test.metagene_unflip_regions_not_previously_flipped <- function() {
 mg <- metagene:::metagene$new(bam_files = bam_files[1],
                                regions = regions_strand)
 checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
 mg$produce_table()
 tab1 <- mg$get_table()
 checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
 mg$unflip_regions()
 tab2 <- mg$get_table()
 checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
 # Compare the table
 checkTrue(identical(tab1, tab2) == TRUE)
}
#
## Valid case previously flipped
test.metagene_unflip_regions_previously_flipped <- function() {
 mg <- metagene:::metagene$new(bam_files = bam_files[1],
                                regions = regions_strand)
 checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
 mg$produce_table(flip_regions = TRUE)
 tab1 <- mg$get_table()
 checkIdentical(mg$get_params()[["flip_regions"]], TRUE)
 mg$unflip_regions()
 tab2 <- mg$get_table()
 checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
 mg$unflip_regions()
 tab3 <- mg$get_table()
 checkIdentical(mg$get_params()[["flip_regions"]], FALSE)
 # Compare the table
 # print(tab1)
 # print(tab2)
 # print(tab3)
 checkTrue(identical(tab1, tab2) == FALSE)
 checkTrue(identical(tab2, tab3) == TRUE)
}
