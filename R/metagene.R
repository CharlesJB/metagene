#' A class to manage metagene analysis.
#'
#' This class will allow to load, convert and normalize alignments and regions
#' files/data. Once the data is ready, the user can then chose to produce
#' metagene plots on the data (or a subset of the data).
#'
#' @section Constructor:
#' \describe{
#'     \item{}{\code{mg <- metagene$new(regions, bam_files, padding_size = 0,
#'                              cores = SerialParam(), verbose = FALSE,
#'                              force_seqlevels = FALSE)}}
#'     \item{regions}{Either a \code{vector} of BED, narrowPeak or broadPeak
#'                    filenames, a \code{GRanges} object or a \code{GRangesList}
#'                    object.}
#'     \item{bam_files}{A \code{vector} of BAM filenames. The BAM files must be
#'                      indexed. i.e.: if a file is named file.bam, there must
#'                      be a file named file.bam.bai in the same directory.}
#'     \item{padding_size}{The regions will be extended on each side by the
#'                         value of this parameter. The padding_size must be a
#'                         non-negative integer. Default = 0.}
#'     \item{cores}{The number of cores available to parallelize the analysis.
#'                  Either a positive integer or a \code{BiocParallelParam}.
#'                  Default: \code{SerialParam()}.}
#'     \item{verbose}{Print progression of the analysis. A logical constant.
#'                    Default: \code{FALSE}.}
#'     \item{force_seqlevels}{If \code{TRUE}, Remove regions that are not found
#'                            in bam file header. Default: \code{FALSE}.}
#' }
#'
#'     \code{metagene$new} returns a \code{metagene} object that contains the
#'         coverages for every BAM files in the regions from the \code{regions}
#'         param.
#'
#' @return
#' \code{metagene$new} returns a \code{metagene} object which contains the
#' normalized coverage values for every regions and for every BAM files.
#'
#' @section Methods:
#' \describe{
#'     \item{}{\code{mg$plot(region_names = NULL, exp_names = NULL,
#'                   title = NULL, x_label = NULL)}}
#'     \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                         all the regions are returned. Default: \code{NULL}.}
#'     \item{exp_names}{The names of the experiments to extract. If a design was
#'                      added to the \code{metagene} object, \code{exp_names}
#'                      correspond to the column names in the design, otherwise
#'                      \code{exp_names} corresponds to the BAM name or the BAM
#'                      filename. If \code{NULL}, all the experiments are
#'                      returned. Default: \code{NULL}.}
#'     \item{title}{A title to add to the graph. If \code{NULL}, will be
#'                         automatically created. Default: NULL}
#'     \item{x_label}{X-axis label to add to the metagene plot. If \code{NULL},
#'                    metagene will use generic label. Default: \code{NULL}.}
#' }
#' \describe{
#'     \item{}{\code{mg$produce_matrices(design, bin_count, noise_removal,
#'                   normalization, flip_regions, bin_size = NULL}}
#'     \item{design}{A \code{data.frame} that describe to experiment to plot.
#'                   see \code{plot} function for more details. \code{NA} can be
#'                   used keep previous design value. Default: \code{NA}.}
#'     \item{bin_count}{The number of bin to create. \code{NA} can be used to
#'                      keep previous bin_count value. A bin_count value of 100
#'                      will be used if no value is specified. Default:
#'                      \code{NA}.}
#'     \item{noise_removal}{The algorithm to use to remove control(s). Possible
#'                          values are \code{NA}, \code{NULL} or "NCIS". By
#'                          default, value is \code{NULL}. Use \code{NA} keep
#'                          previous \code{noise_removal} value (i.e. if
#'                          \code{produce_matrices} was called before). See
#'                          Liand and Keles 2012 for the NCIS algorithm.}
#'     \item{normalization}{The algorithm to use to normalize samples. Possible
#;                          values are \code{NA}, \code{NULL} or "RPM". By
#'                          default, value is \code{NULL} and no normalization
#'                          will be performed. Use \code{NA} keep
#'                          previous \code{normalization} value (i.e. if
#'                          \code{produce_matrices} was called before).}
#'     \item{flip_regions}{Should regions on negative strand be flip_regions?
#'                         Default: \code{FALSE}.}
#'     \item{bin_size}{Deprecated.}
#' }
#' \describe{
#'     \item{}{\code{mg$produce_data_frame(stat = "bootstrap", range = NULL
#'                   ...)}}
#'     \item{stat}{The stat to use to calculate the values of the ribbon in the
#'                 metagene plot. Must be "bootstrap" or "basic". "bootstrap"
#'                 will estimate the range of average ("mean" or "median") that
#'                 could have produce the observed distribution of values for
#'                 each bin. With the "basic" approach, the ribbon will
#'                 represent the range of the values between
#'                 \code{1 - alpha / 2} and \code{alpha / 2} (see \code{...}
#'                 param.}
#'     \item{range}{Deprecated.}
#'     \item{...}{Extra params for the calculation of the ribbon values. See
#'                following param descriptions.}
#'     \item{alpha}{The range of the estimation to be shown with the ribbon.
#'                  \code{1 - alpha / 2} and \code{alpha / 2} will be used.
#'                  Default: 0.05.}
#'     \item{average}{The function to use to summarize the values of each
#'                    bins. "mean" or "median". Default: "mean".}
#'     \item{sample_count}{With "bootstrap" only. The number of draw to do
#'                         in the bootstrap calculation. Default: 1000.}
#' }
#' \describe{
#'     \item{}{mg$get_params()}
#' }
#' \describe{
#'     \item{}{mg$get_design()}
#' }
#' \describe{
#'     \item{}{mg$get_regions(region_names = NULL)}
#'     \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                         all the regions are returned. Default: \code{NULL}.}
#' }
#' \describe{
#'     \item{}{mg$get_matrices(region_names = NULL, exp_name = NULL)}
#'     \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                         all the regions are returned. Default: \code{NULL}.}
#'     \item{exp_names}{The names of the experiments to extract. If a design was
#'                      added to the \code{metagene} object, \code{exp_names}
#'                      correspond to the column names in the design, otherwise
#'                      \code{exp_names} corresponds to the BAM name or the BAM
#'                      filename. If \code{NULL}, all the experiments are
#'                      returned. Default: \code{NULL}.}
#' }
#' \describe{
#'     \item{}{mg$get_data_frame(region_names = NULL, exp_name = NULL)}
#'     \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                         all the regions are returned. Default: \code{NULL}.}
#'     \item{exp_names}{The names of the experiments to extract. If a design was
#'                      added to the \code{metagene} object, \code{exp_names}
#'                      correspond to the column names in the design, otherwise
#'                      \code{exp_names} corresponds to the BAM name or the BAM
#'                      filename. If \code{NULL}, all the experiments are
#'                      returned. Default: \code{NULL}.}
#' }
#' \describe{
#'     \item{}{get_plot = function()}
#' }
#' \describe{
#'     \item{}{get_raw_coverages = function(filenames)}
#'     \item{filenames}{The name of the file to extract raw coverages. Can be
#'                      the filename with the extension of the name of the bam
#'                      file (if a named bam files was used during the creation
#'                      of the metagene object). If \code{NULL}, returns the
#'                      coverage of every bam files. Default: \code{NULL}.}
#' }
#' \describe{
#'     \item{}{get_normalized_coverages = function(filenames)}
#'     \item{filenames}{The name of the file to extract normalized coverages (in
#'                      RPM). Can be the filename with the extension of the name
#'                      of the bam file (if a named bam files was used during
#'                      the creation of the metagene object). If \code{NULL},
#'                      returns the coverage every bam files. Default:
#'                      \code{NULL}.}
#' }
#' \describe{
#'     \item{}{\code{mg$export(bam_file, region, file)}}
#'     \item{bam_file}{The name of the bam file to export.}
#'     \item{region}{The name of the region to export.}
#'     \item{file}{The name of the ouput file.}
#' }
#' \describe{
#'     \item{}{\code{mg$add_design(design = NULL, check_bam_files = FALSE)}}
#'     \item{design}{A \code{data.frame} that describe to experiment to plot.
#'                   See \code{plot} function for more details. \code{NA} can be
#'                   used keep previous design value. Default: \code{NA}.}
#'     \item{check_bam_files}{Force check that all the bam files from the first
#'                            columns of the design are present in current
#'                            metagene object. Default: \code{FALSE}}
#' }
#' \describe{
#'     \item{}{\code{mg$unflip_regions()}}
#' }
#'
#' \describe{
#'     \item{}{\code{mg$flip_regions()}}
#' }
#' \describe{
#'     \item{}{\code{mg$unflip_regions()}}
#' }
#'
#' @examples
#' region <- get_demo_regions()[1]
#' bam_file <- get_demo_bam_files()[1]
#' mg <- metagene$new(regions = region, bam_files = bam_file)
#' \dontrun{
#'     df <- metagene$plot()
#' }
#'
#' @importFrom R6 R6Class
#' @export
#' @format A metagene experiment manager

metagene <- R6Class("metagene",
    public = list(
    # Methods
        initialize = function(regions, bam_files, padding_size = 0,
                              cores = SerialParam(), verbose = FALSE,
                              force_seqlevels = FALSE) {
            # Check params...
            private$check_param(regions = regions, bam_files = bam_files,
                                padding_size = padding_size,
                                cores = cores, verbose = verbose,
                                force_seqlevels = force_seqlevels)
            # Save params
            private$parallel_job <- Parallel_Job$new(cores)
            private$params[["padding_size"]] <- padding_size
            private$params[["verbose"]] <- verbose
            if (is.null(names(bam_files))) {
                new_names <- tools::file_path_sans_ext(basename(bam_files))
                names(bam_files) <- new_names
            }
            private$params[["bin_count"]] <- 100
            private$params[["bam_files"]] <- bam_files
            private$params[["force_seqlevels"]] <- force_seqlevels
            private$params[["flip_regions"]] <- FALSE

            # Prepare bam files
            private$print_verbose("Prepare bam files...")
            private$bam_handler <- Bam_Handler$new(bam_files, cores = cores)

            # Prepare regions
            private$print_verbose("Prepare regions...")
            private$regions <- private$prepare_regions(regions)

            # Parse bam files
            private$print_verbose("Parse bam files...\n")
            private$print_verbose("  coverages...\n")
            private$coverages <- private$produce_coverages()
        },
        get_bam_count = function(filename) {
            # Parameters validation are done by Bam_Handler object
            private$bam_handler$get_aligned_count(filename)
        },
        get_params = function() {
            private$params
        },
        get_design = function() {
            private$design
        },
        get_regions = function(region_names = NULL) {
            if (is.null(region_names)) {
                private$regions
            } else {
                new_names <- tools::file_path_sans_ext(basename(region_names))
                region_names <- new_names
                stopifnot(all(region_names %in% names(private$regions)))
                private$regions[region_names]
            }
        },
        get_matrices = function(region_names = NULL, exp_names = NULL) {
            if (length(private$matrices) == 0) {
                NULL
            } else if (is.null(region_names) & is.null(exp_names)) {
                private$matrices
            } else {
                if (!is.null(region_names)) {
                    stopifnot(is.character(region_names))
                    region_names <-
                        tools::file_path_sans_ext(basename(region_names))
                    stopifnot(all(region_names %in% names(private$matrices)))
                } else {
                    region_names <- names(private$regions)
                }
                if (!is.null(exp_names)) {
                    stopifnot(is.character(exp_names))
                    exp_names <- private$get_bam_names(exp_names)
                    stopifnot(all(exp_names %in% names(private$matrices[[1]])))
                } else {
                    exp_names <- colnames(private$design)[-1]
                }
                matrices <- private$matrices[region_names]
                lapply(matrices, `[`, exp_names)
            }
        },
        get_data_frame = function(region_names = NULL, exp_names = NULL) {
            if (nrow(private$df) == 0) {
                NULL
            } else if (is.null(region_names) & is.null(exp_names)) {
                private$df
            } else {
                if (!is.null(region_names)) {
                    stopifnot(is.character(region_names))
                    region_names <- basename(region_names)
                    region_names <- tools::file_path_sans_ext(region_names)
                    stopifnot(all(region_names %in% names(private$matrices)))
                } else {
                    region_names <- names(private$regions)
                }
                if (!is.null(exp_names)) {
                    stopifnot(is.character(exp_names))
                    exp_names <- private$get_bam_names(exp_names)
                    stopifnot(all(exp_names %in% names(private$matrices[[1]])))
                } else {
                    exp_names <- colnames(private$design)[-1]
                }
                eg <- expand.grid(exp_names, region_names)
                groups <- do.call(paste, c(eg, sep = "_"))
                i <- private$df$group %in% groups
                private$df[i,]
            }
        },
        get_plot = function() {
            if (is.character(private$graph)) {
                NULL
            } else {
                private$graph
            }
        },
        get_raw_coverages = function(filenames = NULL) {
            if (is.null(filenames)) {
                private$coverages
            } else {
                stopifnot(is.character(filenames))
                stopifnot(length(filenames) > 0)
                bam_names <- private$get_bam_names(filenames)
                stopifnot(length(bam_names) == length(filenames))
                private$coverages[bam_names]
            }
        },
        get_normalized_coverages = function(filenames = NULL) {
            normalize_coverage <- function(filename) {
            count <- private$bam_handler$get_aligned_count(filename)
                weight <- 1 / (count / 1000000)
                coverages[[filename]] <- coverages[[filename]] * weight
        }
        coverages <- self$get_raw_coverages(filenames)
            coverage_names <- names(coverages)
            coverages <-
                private$parallel_job$launch_job(data = coverage_names,
                                                FUN = normalize_coverage)
            names(coverages) <- coverage_names
            coverages
        },
        add_design = function(design, check_bam_files = FALSE) {
            private$design = private$fetch_design(design, check_bam_files)
        },
		produce_rna_seq_matrices = function() {
			stitch_gene <- function(cov, gr) {
			    chr <- unique(as.character(GenomeInfoDb::seqnames(gr)))
			    view <- Views(cov[[chr]], start(gr), end(gr))
			    as.numeric(do.call("c", viewApply(view, function(x) x)))
			}
			coverages <- self$get_normalized_coverages()
			for (region in names(self$get_regions())) {
				gr <- self$get_regions()[[region]]
				matrices[[region]] <- list()
				for (bam_name in names(coverages)) {
					matrices[[region]][[bam_name]] <- list()

					matrices[[region]][[design_name]][["input"]] <-
						stitch_gene(coverages[[bam_name]], gr)
				}
			}
		},
        produce_matrices = function(design = NA, bin_count = NA,
                                    noise_removal = NA, normalization = NA,
                                    flip_regions = FALSE, bin_size = NULL) {
            if (!is.null(bin_size)) {
                warning("bin_size is now deprecated. Please use bin_count.")
                bin_size <- NULL
            }
            design = private$fetch_design(design)
            private$check_produce_matrices_params(bin_count = bin_count,
                                                  bin_size = bin_size,
                                                  design = design,
                                                  noise_removal = noise_removal,
                                                  normalization = normalization,
                                                  flip_regions = flip_regions)
            bin_count <- private$get_param_value(bin_count, "bin_count")
            noise_removal <- private$get_param_value(noise_removal,
                                                     "noise_removal")
            normalization <- private$get_param_value(normalization,
                                                     "normalization")
            if (is.null(bin_count)) {
                bin_count = 100
            }
            if (private$matrices_need_update(design = design,
                                             bin_count = bin_count,
                                             bin_size = bin_size,
                                             noise_removal = noise_removal,
                                             normalization = normalization)) {
                coverages <- private$coverages
                if (!is.null(noise_removal)) {
                  coverages <- private$remove_controls(coverages, design)
                } else {
                  coverages <- private$merge_chip(coverages, design)
                }
                if (!is.null(normalization)) {
                  coverages <- private$normalize_coverages(coverages, design)
                }
                matrices <- list()
                for (region in names(self$get_regions())) {
                    matrices[[region]] <- list()
                    for (design_name in colnames(design)[-1]) {
                        matrices[[region]][[design_name]] <- list()

                        matrices[[region]][[design_name]][["input"]] <-
                        private$get_matrix(coverages[[design_name]], region,
                                           bin_count)
                    }
                }
                private$matrices <- matrices
                private$params[["bin_count"]] <- bin_count
                private$params[["noise_removal"]] <- noise_removal
                private$params[["normalization"]] <- normalization
                private$design <- design
            }
            if (flip_regions == TRUE) {
                self$flip_regions()
            } else {
                self$unflip_regions()
            }
            invisible(self)
        },
		produce_data_frame_rna_seq = function(stat = "bootstrap", gene_name, ...) {
			stopifnot(is.character(gene_name))
			stopifnot(length(gene_name) == 1)
			stopifnot(gene_name %in% names(self$get_regions()))
            stopifnot(is.character(stat))
            stopifnot(length(stat) == 1)
            stopifnot(stat %in% c("bootstrap", "basic"))

			gr <- self$get_regions()[[gene_name]]
			range <- c(1, sum(width(gr)))

            # 1. Get the correctly formatted matrices
            if (length(private$matrices) == 0) {
                self$produce_rna_seq_matrices()
            }

            # 2. Produce the data.frame
            if (private$data_frame_need_update(stat, ...) == TRUE) {
                sample_size = NULL
                if (stat == "bootstrap") {
                    stat <- Bootstrap_Stat
                    sample_size <- sapply(private$matrices, sapply, sapply,
                                          nrow)
                    sample_size <- as.integer(min(unlist(sample_size)))
                } else {
                    stat <- Basic_Stat
                }
                df <- data.frame(group = character(), position = numeric(),
                                 value = numeric(), qinf = numeric(),
                                 qsup = numeric())

                m <- private$matrices
                regions <- names(m)
                stopifnot(length(unique(lapply(m, names))) == 1)
                exp_names <- unique(lapply(m, names))[[1]]

                combinations <- expand.grid(regions, exp_names)
                combinations <- split(combinations, 1:nrow(combinations))

                get_statistics <- function(combination) {
                    region <- as.character(combination[1,1])
                    exp_name <- as.character(combination[1,2])
                    group_name <- paste(exp_name, region, sep = "_")
                    message(group_name)
                    data <- m[[region]][[exp_name]][["input"]]
                    ctrl <- m[[region]][[exp_name]][["ctrl"]]
                    current_stat <- ""
                    if (! is.null(sample_size)) {
                        current_stat <- stat$new(data = data, ctrl = ctrl,
                                                 sample_size = sample_size,
                                                 range = range, ...)
                    } else {
                        current_stat <- stat$new(data = data, ctrl = ctrl,
                                                 range = range, ...)
                    }
                    current_df <- current_stat$get_statistics()
                    current_df <- cbind(rep(group_name, nrow(current_df)),
                                        current_df)
                    colnames(current_df)[1] <- "group"
                    current_df
                }
                df <- private$parallel_job$launch_job(data = combinations,
                                                      FUN = get_statistics)
                df <- do.call("rbind", df)
                private$df <- df
            }
            invisible(self)
		},
        produce_data_frame = function(stat = "bootstrap", range = NULL,
                                      ...) {
            # Prepare range
            if (!is.null(range)) {
                warning("range param in produce_data_frame is now deprecated.")
            }
            bin_count <- private$params[["bin_count"]]
            dist <- floor(bin_count / 2)
            range <- c(-dist, dist)

            stopifnot(is.character(stat))
            stopifnot(length(stat) == 1)
            stopifnot(stat %in% c("bootstrap", "basic"))
            # 1. Get the correctly formatted matrices
            if (length(private$matrices) == 0) {
                self$produce_matrices()
            }

            # 2. Produce the data.frame
            if (private$data_frame_need_update(stat, ...) == TRUE) {
                sample_size = NULL
                if (stat == "bootstrap") {
                    stat <- Bootstrap_Stat
                    sample_size <- sapply(private$matrices, sapply, sapply,
                                          nrow)
                    sample_size <- as.integer(min(unlist(sample_size)))
                } else {
                    stat <- Basic_Stat
                }
                df <- data.frame(group = character(), position = numeric(),
                                 value = numeric(), qinf = numeric(),
                                 qsup = numeric())

                m <- private$matrices
                regions <- names(m)
                stopifnot(length(unique(lapply(m, names))) == 1)
                exp_names <- unique(lapply(m, names))[[1]]

                combinations <- expand.grid(regions, exp_names)
                combinations <- split(combinations, 1:nrow(combinations))

                get_statistics <- function(combination) {
                    region <- as.character(combination[1,1])
                    exp_name <- as.character(combination[1,2])
                    group_name <- paste(exp_name, region, sep = "_")
                    message(group_name)
                    data <- m[[region]][[exp_name]][["input"]]
                    ctrl <- m[[region]][[exp_name]][["ctrl"]]
                    current_stat <- ""
                    if (! is.null(sample_size)) {
                        current_stat <- stat$new(data = data, ctrl = ctrl,
                                                 sample_size = sample_size,
                                                 range = range, ...)
                    } else {
                        current_stat <- stat$new(data = data, ctrl = ctrl,
                                                 range = range, ...)
                    }
                    current_df <- current_stat$get_statistics()
                    current_df <- cbind(rep(group_name, nrow(current_df)),
                                        current_df)
                    colnames(current_df)[1] <- "group"
                    current_df
                }
                df <- private$parallel_job$launch_job(data = combinations,
                                                      FUN = get_statistics)
                df <- do.call("rbind", df)
                private$df <- df
            }
            invisible(self)
        },
        plot = function(region_names = NULL, exp_names = NULL, title = NULL, x_label = NULL) {
            # 1. Get the correctly formatted matrices
            if (length(private$matrices) == 0) {
                self$produce_matrices()
            }

            # 2. Produce the data frame
            if (nrow(private$df) == 0) {
                self$produce_data_frame()
            }
            df <- self$get_data_frame(region_names = region_names,
                                      exp_names = exp_names)
            # 3. Produce the graph
            if (is.null(title)) {
                title <- paste(unique(private$df[["group"]]), collapse=" vs ")
            }
            p <- private$plot_graphic(df = df, title = title, x_label = x_label)
            print(p)
            private$graph <- p
            invisible(self)
        },
        export = function(bam_file, region, file) {
            region <- tools::file_path_sans_ext(basename(region))
            region <- private$regions[[region]]
            param <- Rsamtools::ScanBamParam(which = region)
            alignments <- GenomicAlignments::readGAlignments(bam_file,
                                                                param = param)
            weight <- 1 - private$bam_handler$get_rpm_coefficient(bam_file)
            seqlevels(alignments) <- seqlevels(region)
            # TODO: don't use the weight param of coverage
            coverage <- GenomicAlignments::coverage(alignments, weight=weight)
            rtracklayer::export(coverage, file, "BED")
            invisible(coverage)
        },
        flip_regions = function() {
            if (private$params[["flip_regions"]] == FALSE) {
                private$flip_matrices()
                private$params[["flip_regions"]] <- TRUE
            }
            invisible(self)
        },
        unflip_regions = function() {
            if (private$params[["flip_regions"]] == TRUE) {
                private$flip_matrices()
                private$params[["flip_regions"]] <- FALSE
            }
            invisible(self)
        }
    ),
        private = list(
        params = list(),
        regions = GRangesList(),
        matrices = list(),
        design = data.frame(),
        coverages = list(),
        df = data.frame(),
        graph = "",
        bam_handler = "",
        parallel_job = "",
        check_param = function(regions, bam_files, padding_size,
                               cores, verbose, force_seqlevels) {
            # Check parameters validity
            if (!is.logical(verbose)) {
                stop("verbose must be a logicial value (TRUE or FALSE)")
            }
            if (!is.logical(force_seqlevels)) {
                stop("force_seqlevels must be a logicial value (TRUE or FALSE)")
            }
            if (!(is.numeric(padding_size) || is.integer(padding_size)) ||
                padding_size < 0 || as.integer(padding_size) != padding_size) {
                stop("padding_size must be a non-negative integer")
            }
            isBiocParallel = is(cores, "BiocParallelParam")
            isInteger = ((is.numeric(cores) || is.integer(cores)) &&
                             cores > 0 &&  as.integer(cores) == cores)
            if (!isBiocParallel && !isInteger) {
                stop(paste0("cores must be a positive numeric or ",
                            "BiocParallelParam instance"))
            }
            if (!is.vector(bam_files, "character")) {
                stop("bam_files must be a vector of BAM filenames")
            }
            if (!all(sapply(bam_files, file.exists))) {
                stop("At least one BAM file does not exist")
            }
            if (!all(sapply(paste0(bam_files, ".bai"), file.exists))) {
                stop("All BAM files must be indexed")
            }
            if (!is(regions, "GRangesList") && !is.character(regions)
                && !is(regions, "GRanges") && !is.list(regions)) {
                stop(paste0("regions must be either a vector of BED ",
                    "filenames, a GRanges object or a GrangesList object"))
            }
            # Validation specific to regions as a vector
            if (is.character(regions) && !all(sapply(regions, file.exists))) {
                stop("regions must be a list of existing BED files")
            }
        },
        check_design = function(design, check_bam_files = FALSE) {
            stopifnot(is.logical(check_bam_files))
            if(!is.null(design) && !is.data.frame(design) &&
               !identical(design, NA)) {
                stop("design must be a data.frame object, NULL or NA")
            }
            if (is.data.frame(design)) {
                if(ncol(design) < 2) {
                    stop("design must have at least 2 columns")
                }
                if (!(is.character(design[,1]) || is.factor(design[,1]))) {
                    stop("The first column of design must be BAM filenames")
                }
                if (!all(apply(design[, -1, drop=FALSE], MARGIN=2,
                               is.numeric))) {
                    stop(paste0("All design column, except the first one, must",
                                    " be in numeric format"))
                }
                if (check_bam_files == TRUE) {
                    samples <- as.character(design[,1])
                    if (!all(private$check_bam_files(samples))) {
                        stop("Design contains bam files absent from metagene.")
                    }
                }
            }
        },
        check_produce_matrices_params = function(bin_count, bin_size, design,
                                                 noise_removal, normalization,
                                                 flip_regions) {
            # At least one file must be used in the design
            if (!identical(design, NA)) {
                if (!is.null(design)) {
                    if (sum(rowSums(design[ , -1, drop=FALSE]) > 0) == 0) {
                        stop("At least one BAM file must be used in the design")
                    }
                }
            }
            # Test only BAM file used in the design
            if (!identical(design, NA)) {
                if(!is.null(design) &&
                    !all(apply(design[rowSums(design[ , -1, drop=FALSE]) > 0, 1,
                                    drop=FALSE], MARGIN = 2,
                            FUN=private$check_bam_files))) {
                    stop("At least one BAM file does not exist")
                }
            }
            # bin_count should be a numeric value without digits
            if (!identical(bin_count, NA)) {
                if (!is.null(bin_count)) {
                    if (!is.numeric(bin_count) || bin_count < 0 ||
                        as.integer(bin_count) != bin_count) {
                        stop("bin_count must be NULL or a positive integer")
                    }
                }
            }
            if (!identical(noise_removal, NA)) {
                if (!is.null(noise_removal)) {
                    if (!noise_removal %in% c("NCIS", "RPM")) {
                        msg <- 'noise_removal must be NA, NULL, "NCIS" or '
                        msg <- paste0(msg, '"RPM".')
                        stop(msg)
                    }
                }
            }
            if (!identical(normalization, NA)) {
                if (!is.null(normalization)) {
                    if (!normalization == "RPM") {
                        msg <- "normalization must be NA, NULL or \"RPM\"."
                        stop(msg)
                    }
                }
            }
            if (!is.logical(flip_regions)) {
                msg <- "flip_regions must be a logical."
                stop(msg)
            }
        },
        matrices_need_update = function(design, bin_count, bin_size,
                                        noise_removal, normalization) {
            if (!identical(private$design, design)) {
                return(TRUE)
            }
            if (!identical(private$params[["bin_count"]], bin_count)) {
                return(TRUE)
            }
            if (!identical(private$params[["noise_removal"]], noise_removal)) {
                return(TRUE)
            }
            if (!identical(private$params[["normalization"]], normalization)) {
                return(TRUE)
            }
            return(FALSE)
        },
        data_frame_need_update = function(stat = NA, alpha = NA, average = NA,
                                          sample_count = NA) {
            need_update = FALSE
            # Fetch saved values
            stat = private$get_param_value(stat, "stat")
            alpha = private$get_param_value(alpha, "alpha")
            average = private$get_param_value(average, "average")
            sample_count = private$get_param_value(sample_count, "sample_count")

            # Add default, if needed
            if (is.null(stat)) {
                stat <- "bootstrap"
            }
            if (is.null(alpha)) {
                alpha <- 0.05
            }
            if (is.null(average)) {
                average <- "mean"
            }
            if (is.null(sample_count)) {
                sample_count <- 1000
            }
            if (nrow(private$df) == 0) {
                need_update <- TRUE
                private$params[["sample_count"]] <- sample_count
                private$params[["average"]] <- average
                private$params[["alpha"]] <- alpha
                private$params[["stat"]] <- stat
            } else {
                # Check if data frame need update
                if (!identical(private$params[["stat"]], stat)) {
                    need_update <- TRUE
                    private$params[["stat"]] <- stat
                }
                if (!identical(private$params[["alpha"]], alpha)) {
                    need_update <- TRUE
                    private$params[["alpha"]] <- alpha
                }
                if (!identical(private$params[["average"]], average)) {
                    need_update <- TRUE
                    private$params[["average"]] <- average
                }
                if (!identical(private$params[["sample_count"]],
                               sample_count)) {
                    need_update <- TRUE
                    private$params[["sample_count"]] <- sample_count
                }
            }
            need_update
        },
        get_param_value = function(param_value, param_name) {
            param_name <- as.character(param_name)
            if (identical(param_value, NA)) {
                if (! param_name %in% names(private$params)) {
                    param_value <- NULL
                } else {
                    param_value <- private$params[[param_name]]
                }
            }
            return(param_value)
        },
        print_verbose = function(to_print) {
            if (private$params[["verbose"]]) {
                cat(paste0(to_print, "\n"))
            }
        },
        get_matrix = function(coverages, region, bcount) {
            gr <- private$regions[[region]]
            grl <- split(gr, GenomeInfoDb::seqnames(gr))
            i <- vapply(grl, length, numeric(1)) > 0
            m <- do.call("rbind", lapply(grl[i], private$get_view_means,
                                         bcount = bcount, cov = coverages))
        },
        get_view_means = function(gr, bcount, cov) {
            chr <- unique(as.character(GenomeInfoDb::seqnames(gr)))
            gr <- intoNbins(gr, bcount)
            stopifnot(length(chr) == 1)
            views <- Views(cov[[chr]], start(gr), end(gr))
            matrix(viewMeans(views), ncol = bcount, byrow = TRUE)
        },
        prepare_regions = function(regions) {
            if (class(regions) == "character") {
                names <- sapply(regions, function(x)
                    file_path_sans_ext(basename(x)))
                import_file <- function(region) {
            ext <- tolower(tools::file_ext(region))
                    if (ext == "narrowpeak") {
                        extraCols <- c(signalValue = "numeric",
                                       pValue = "numeric", qValue = "numeric",
                                       peak = "integer")
                        rtracklayer::import(region, format = "BED",
                                            extraCols = extraCols)
                } else if (ext == "broadpeak") {
                        extraCols <- c(signalValue = "numeric",
                                       pValue = "numeric", qValue = "numeric")
                        rtracklayer::import(region, format = "BED",
                                            extraCols = extraCols)
                    } else {
                        rtracklayer::import(region)
                    }
                }
                regions <- private$parallel_job$launch_job(data = regions,
                                                           FUN = import_file)
                names(regions) <- names
            } else if (class(regions) == "GRanges") {
                regions <- GRangesList("regions" = regions)
            } else if (class(regions) == "list") {
                regions <- GRangesList(regions)
        }
            if (is.null(names(regions))) {
            names(regions) <- sapply(seq_along(regions), function(x) {
                    paste("region", x, sep = "_")
                              })
        }
            # TODO: Check if there is a id column in the mcols of every ranges.
            #       If not, add one by merging seqnames, start and end.
            GRangesList(lapply(regions, function(x) {
                # Add padding
                start(x) <- start(x) - private$params$padding_size
                start(x)[start(x) < 0] <- 1
                end(x) <- end(x) + private$params$padding_size
                # Clean seqlevels
            x <- sortSeqlevels(x)
    #            seqlevels(x) <- unique(as.character(seqnames(x)))
                x
            }))
        },
        produce_coverages = function() {
        regions <- GenomicRanges::reduce(BiocGenerics::unlist(private$regions))
            res <- private$parallel_job$launch_job(
                        data = private$params[["bam_files"]],
                        FUN = private$bam_handler$get_coverage,
                        regions = regions,
                        force_seqlevels= private$params[["force_seqlevels"]])
            names(res) <- names(private$params[["bam_files"]])
            lapply(res, GenomeInfoDb::sortSeqlevels)
        },
        plot_graphic = function(df, title, x_label) {
            # Prepare x label
            if (is.null(x_label)) {
                x_label <- "Distance in bins from regions center"
            }

            # Prepare y label
            y_label <- "Mean coverage"
            if (is.null(private$params[["normalization"]])) {
                y_label <- paste(y_label, "(raw)")
            } else {
                y_label <- paste(y_label("(RPM)"))
            }

            # Produce plot
            p <- plot_metagene(df) +
                ylab(y_label) +
                xlab(x_label) +
                ggtitle(title)
            p
        },
        fetch_design = function(design, check_bam_files = FALSE) {
            private$check_design(design = design, check_bam_files)
            get_complete_design <- function() {
                bam_files <- names(private$params[["bam_files"]])
                design <- data.frame(bam_files = bam_files)
                for (bam_file in names(private$coverages)) {
                    colname <- file_path_sans_ext(basename(bam_file))
                    design[[colname]] <- rep(0, length(bam_files))
                    i <- bam_files == bam_file
                    design[[colname]][i] <- 1
                }
                design
            }
            if (is.null(design)) {
                return(get_complete_design())
            }
            if (identical(design, NA)) {
                if (all(dim(private$design) == c(0, 0))) {
                    return(get_complete_design())
                } else {
                    return(private$design)
                }
            }
            design[,1] <- as.character(design[,1])
            return(design)
        },
        remove_controls = function(coverages, design) {
            results <- list()
            for (design_name in colnames(design)[-1]) {
                i <- design[[design_name]] == 1
                j <- design[[design_name]] == 2
            chip_bam_files <- as.character(design[,1][i])
            chip_names <- private$get_bam_names(chip_bam_files)
                input_bam_files <- as.character(design[,1][j])
            input_names <- private$get_bam_names(input_bam_files)
                chip_coverages <- coverages[chip_names]
                chip_coverages <- Reduce("+", chip_coverages)
                if (length(input_bam_files) > 0) {
                    noise_ratio <-
                        private$bam_handler$get_noise_ratio(chip_names,
                                                            input_names)
                    input_coverages <- coverages[input_names]
                    input_coverages <- noise_ratio * Reduce("+",
                                                            input_coverages)
                    results[design_name] <- chip_coverages - input_coverages
                    i <- results[[design_name]] < 0
                    results[[design_name]][i] <- 0
                } else {
                    results[design_name] <- chip_coverages
                }
            }
            results
        },
        normalize_coverages = function(coverages, design) {
            for (design_name in colnames(design)[-1]) {
                bam_files <- as.character(design[,1][1])
                counts <- lapply(bam_files,
                                 private$bam_handler$get_aligned_count)
                count <- do.call("+", counts)
                weight <- 1 / (count / 1000000)
                coverages[[design_name]] <- coverages[[design_name]] * weight
            }
            coverages
        },
        merge_chip = function(coverages, design) {
            result <- list()
            for (design_name in colnames(design)[-1]) {
                i <- design[[design_name]] == 1
                bam_files <- as.character(design[,1][i])
            bam_names <- private$get_bam_names(bam_files)
                cov <- coverages[bam_names]
                result[[design_name]] <- Reduce("+", cov)
            }
            result
        },
        flip_matrices = function() {
            for (region_name in names(private$matrices)) {
                region <- private$regions[[region_name]]
                i <- as.logical(strand(region) == "-")
                flip <- function(x) {
                    x[i,] <- x[i,ncol(x):1]
                    x
                }
                m <- private$matrices[[region_name]]
                m <- lapply(m, lapply, flip)
                private$matrices[[region_name]] <- m
            }
        },
        get_bam_names = function(filenames) {
            if (all(filenames %in% colnames(private$design)[-1])) {
                filenames
            } else {
                stopifnot(private$check_bam_files(filenames))
                vapply(filenames,
                   private$bam_handler$get_bam_name,
                   character(1))
            }
        },
        check_bam_files = function(bam_files) {
            all(vapply(bam_files,
            function(x) {
                !is.null((private$bam_handler$get_bam_name(x)))
            },
            logical(1)))
        }
    )
)
