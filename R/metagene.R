#' A class to manage metagene analysis.
#'
#' This class will allow to load, convert and normalize alignments and regions
#' files/data. Once the data is ready, the user can then chose to produce
#' metagene plots on the data (or a subset of the data).
#'
#' @section Constructor:
#' \describe{
#'   \item{}{\code{mg <- metagene$new(regions, bam_files, padding_size = 0,
#'                              cores = SerialParam(), verbose = FALSE,
#'                              force_seqlevels = FALSE)}}
#'   \item{regions}{Either a \code{vector} of BED filenames, a \code{GRanges}
#'                  object or a \code{GRangesList} object.}
#'   \item{bam_files}{A \code{vector} of BAM filenames. The BAM files must be 
#'                    indexed. i.e.: if a file is named file.bam, there must be 
#'                    a file named file.bam.bai in the same directory.}
#'   \item{padding_size}{The regions will be extended on each side by the value
#'                       of this parameter. The padding_size must be a 
#'                       non-negative integer. Default = 0.}
#'   \item{cores}{The number of cores available to parallelize the analysis.
#'                Either a positive integer or a \code{BiocParallelParam}.
#'                Default: \code{SerialParam()}.}
#'   \item{verbose}{Print progression of the analysis. A logical constant. 
#'                  Default: \code{FALSE}.}
#'   \item{force_seqlevels}{If \code{TRUE}, Remove regions that are not found in
#'                          bam file header. Default: \code{FALSE}.}
#' }
#'
#'  \code{metagene$new} returns a \code{metagene} object that contains the
#'  coverages for every BAM files in the regions from the \code{regions} param.
#'
#' @return
#' \code{metagene$new} returns a \code{metagene} object which contains the
#' normalized coverage values for every regions and for every BAM files.
#'
#' @section Methods:
#' \describe{
#'   \item{}{\code{df <- mg$plot(design = NULL, regions_group = NULL,
#'               bin_size = 100, alpha = 0.05, sample_count = 1000,
#'               title = NULL)}}
#'   \item{design}{A \code{data.frame} that describe to experiment to plot. The
#'                 first column must be the existing BAM filenames. The other 
#'                 columns (at least one is required) represent how the files 
#'                 should be grouped. All those colums should be in numeric 
#'                 format. All the files in the same group will be 
#'                 combined before doing the statistical analysis. For each 
#'                 column, there will be one line with its ribon in 
#'                 the metagene plot. At least one file should be selected (not
#'                 zero). If \code{NULL}, all BAM files passed as argument 
#'                 during \code{metagene} object creation will be used to 
#'                 create a default design. 
#'                 Default = \code{NULL}.
#'
#'                 0: Do not use file.
#'                 1: File is input.
#'                 2: File is control.}
#'   \item{regions_group}{A \code{list} or a \code{vector} of region names to
#'                      include in the analysis. If \code{NULL}, all the
#'                      regions used when creating the
#'                      \code{metagene} object will be used.
#'                      Default: \code{NULL}}
#'   \item{bin_size}{The size of bin to use before calculating the statistics.
#'                   larger \code{bin_size} will reduce the calculation time
#'                   and produce smoother curves. It must be a positive 
#'                   integer. Default = 100.}
#'   \item{alpha}{The condifence interval (CI) to represent with a ribbon. Must
#'                be a value between 0 and 1. With a value of 0.05, the ribbon
#'                will represent 95% of the estimated values of the means.
#'                Default = 0.05.}
#'   \item{sample_count}{The number of draw to do during bootstrap analysis.
#'                       Default: 1000.}
#'   \item{title}{A title to add to the graph. If \code{NULL}, will be
#'                       automatically created. Default: NULL}
#' }
#' \describe{
#'   \item{}{\code{mg$export(bam_file, region, file)}}
#'   \item{bam_file}{The name of the bam file to export.}
#'   \item{region}{The name of the region to export.}
#'   \item{file}{The name of the ouput file.}
#' }
#' \describe{
#'   \item{}{\code{mg$heatmap(region, bam_file, bin_size)}}
#'   \item{region}{The name of the region to export.}
#'   \item{bam_file}{The name of the bam file to export.}
#'   \item{bin_size}{The size of the bin to produce before creating heatmap. It
#'                      must be a positive integer. Default: 10.}
#' }
#'
#' @examples
#'  regions <- get_demo_regions()
#'  bam_files <- get_demo_bam_files()
#'  mg <- metagene$new(regions, bam_files)
#'  \dontrun{
#'  df <- metagene$plot()}
#'
#' @importFrom R6 R6Class
#' @export
#' @format A metagene experiment manager

metagene <- R6Class("metagene",
  public = list(
    # Public members
    params = list(),
    regions = GRangesList(),
    coverages = list(),
    matrices = list(),
    design = data.frame(),
    
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
        self$params[["padding_size"]] <- padding_size
        self$params[["verbose"]] <- verbose
        self$params[["bam_files"]] <- bam_files
        self$params[["force_seqlevels"]] <- force_seqlevels
      
        # Prepare bam files
        private$print_verbose("Prepare bam files...")
        private$bam_handler <- Bam_Handler$new(bam_files, cores = cores)
      
        # Prepare regions
        private$print_verbose("Prepare regions...")
        self$regions <- private$prepare_regions(regions)
      
        # Parse bam files
        private$print_verbose("Parse bam files...\n")
        private$print_verbose("  coverages...\n")
        self$coverages <- private$produce_coverages()
    },
    get_bam_count = function(filename) {
        # Parameters validation are done by Bam_Handler object
        private$bam_handler$get_aligned_count(filename)
    },
    plot = function(design = NULL, regions_group = NULL, bin_count = 100,
                    bin_size = NULL, alpha = 0.05, sample_count = 1000,
		    range = c(-1,1), title = NULL) {
	private$check_plot_param(design = design, regions_group = regions_group,
			       bin_count = bin_count, bin_size = bin_size,
			       alpha = alpha, sample_count = sample_count,
			       range = range, title = title)
        if (!is.null(bin_size)) {
            width <- width(self$regions[[1]][1])
	    bin_count <- floor(width / bin_size)
	}
        self$design <- private$prepare_design(design)

        # 1. Get the correctly formatted matrices
        self$matrices <- private$produce_matrices(regions_group, bin_count)

        # 2. Calculate means and confidence intervals
        sample_size <- as.integer(min(unlist(sapply(self$matrices, sapply,
                                              sapply, nrow))))
        DF <- data.frame(group = character(), position = numeric(),
                       value = numeric(), qinf = numeric(), qsup = numeric())
        for (region in names(self$matrices)) {
            for (bam_file in names(self$matrices[[region]])) {
                group_name <- paste(bam_file, region, sep = "_")
                print(group_name)
                data <- self$matrices[[region]][[bam_file]][["input"]]
                ctrl <- self$matrices[[region]][[bam_file]][["ctrl"]]
                bootstrap_stat <- Bootstrap_Stat$new(data = data, ctrl = ctrl,
                                                sample_size = sample_size, 
                                                sample_count = sample_count,
                                                range = range)
                current_DF <- bootstrap_stat$get_statistics()
                current_DF <- cbind(rep(group_name, nrow(current_DF)), 
                                    current_DF)
                colnames(current_DF)[1] <- "group"
                DF <- rbind(DF, current_DF)
            }
        }
        
        #
        # 3. Test de Friedman
        #
        # Friedman only done when more than 1 curve is present
        # 
        friedman <- NULL
        if (length(unique(DF[["group"]])) > 1) {
            friedman <- mu.friedman.test(DF[["value"]], group=DF[["group"]], 
                             block=DF[["position"]])
        }
        
        # 4. Produce the graph
        #    DF <- metagene:::getDataFrame(bootstrap_result, 
        #                                     range=c(-1,1), binSize=1)
	if (is.null(title)) {
	    title <- paste(unique(DF[["group"]]), collapse=" vs ")
	}
        p <- private$plot_graphic(DF, title = title,
                                binSize = 1, friedman=friedman)
        print(p)
        return(list(DF = DF, friedman_test = friedman, graph = p))
    },
    export = function(bam_file, region, file) {
        region <- tools::file_path_sans_ext(basename(region))
        region <- self$regions[[region]]
        param <- Rsamtools::ScanBamParam(which = region)
        alignments <- GenomicAlignments::readGAlignments(bam_file, 
                                                            param = param)
        weight <- 1 - private$bam_handler$get_rpm_coefficient(bam_file)
        seqlevels(alignments) <- seqlevels(region)
        coverage <- GenomicAlignments::coverage(alignments, weight=weight)
        rtracklayer::export(coverage, file, "BED")
        invisible(coverage)
    },
    heatmap = function(region, bam_file, bin_size = 10) {
        # Check that bin_size is a positive integer
        if (!((is.numeric(bin_size) || is.integer(bin_size)) && 
                bin_size > 0 &&  as.integer(bin_size) == bin_size)) {
            stop("bin_size must be a positive integer")
        }
        
        region <- tools::file_path_sans_ext(basename(region))
        data <- private$get_matrix(region, bam_file, bin_size)
        heatmap.2(log2(data+1), dendrogram = "none", trace = "none",
                labCol = NA, labRow = NA, margins = c(2,2),
                xlab = "position", ylab = "log2(coverages)")
    }
  ),
  private = list(
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
            stop(paste0("regions must be either a vector of BED filenames, a ",
                "GRanges object or a GrangesList object"))
        }
        # Validation specific to regions as a vector
        if (is.character(regions) && !all(sapply(regions, file.exists))) {
            stop("regions must be a list of existing BED files")
        }
    },
    check_design = function(design) {
        if(!is.null(design) && !is.data.frame(design)) {
            stop("design must be a data.frame object")
        }
        if(!is.null(design) && dim(design)[2] < 2) {
            stop("design must have at least 2 columns")
        }
        if (!is.null(design) && !(is.character(design[,1]) ||
                                    is.factor(design[,1]))) {
            stop("The first column of design must be BAM filenames")
        }
        if (!is.null(design) && !all(apply(design[, -1, drop=FALSE],
                                        MARGIN=2, is.numeric))) {
            stop(paste0("All design column, except the first one, must be in ",
                            "numeric format"))
        }
    },
    check_plot_param = function(design, regions_group, bin_count, bin_size,
			       alpha, sample_count, range, title) {
        # Most parameters validation are done in Bootstrap_Stat constructor
        private$check_design(design)
        if (!is.null(regions_group) && !(is.vector(regions_group) ||
                        is.list(regions_group))) {
            stop("regions_group should be a list or a vector")
        }
        if (!is.null(regions_group) &&
                !all(unlist(regions_group) %in%  names(self$regions))) {
            stop(paste0("All elements in regions_group should be regions ",
                        "defined during the creation of metagene object"))
        }
        # At least one file must be used in the design
        if (!is.null(design) &&
                sum(rowSums(design[ , -1, drop=FALSE]) > 0) == 0) {
            stop("At least one BAM file must be used in the design")
        }
        # Test only BAM file used in the design
        if(!is.null(design) &&
            !all(apply(design[rowSums(design[ , -1, drop=FALSE]) > 0, 1,
                            drop=FALSE], MARGIN = 2, FUN=file.exists))) {
            stop("At least one BAM file does not exist")
        }
	# bin_count should be a numeric value without digits
        if (!is.null(bin_count)) {
	    if (!is.numeric(bin_count) || bin_count < 0 ||
		as.integer(bin_count) != bin_count) {
              stop("bin_count must be NULL or a positive integer")
          }
	}
	# bin_size should be a numeric value without digits
        if (!is.null(bin_size)) {
	    if (!is.numeric(bin_size) || bin_size < 0 ||
		as.integer(bin_size) != bin_size) {
              stop("bin_size must be NULL or a positive integer")
          }
	}
	# bin_size should be used only if all regions have the same width
	if (!is.null(bin_size)) {
          if (is.null(regions_group)) {
              regions_group <- names(self$regions)
          }
	  widths <- width(self$regions[names(self$regions) %in% regions_group])
	  widths <- unique(unlist(widths))
	  if (length(widths) != 1) {
	      msg <- "bin_size can only be used if all selected regions have"
	      msg <- paste(msg, "same width")
	      stop(msg)
	  } else {
              modulo <- widths %% bin_size
	      if (modulo != 0) {
                  msg <- paste0("width (", widths, ") is not a multiple of ")
	          msg <- paste0(msg, "bin_size (", bin_size, "), last bin ")
	          msg <- paste0(msg, "will be removed.")
		  warning(msg)
	      }
	  }
	}
    },
    print_verbose = function(to_print) {
        if (self$params[["verbose"]]) {
            cat(paste0(to_print, "\n"))
        }
    },
    prepare_design = function(design = NULL) {
        if (is.null(design)) {
            bam_files <- self$params[["bam_files"]]
            design <- data.frame(bam_files = bam_files)
            for (bam_file in names(self$coverages)) {
                colname <- file_path_sans_ext(basename(bam_file))
                design[[colname]] <- rep(0, length(bam_files))
                i <- bam_files == bam_file
                design[[colname]][i] <- 1
            }
        }
        design
    },
    produce_matrices = function(regions_group, bin_count) {
        matrices <- list()
        if (is.null(regions_group)) {
            regions_group <- names(self$regions)
        }
        for (region in regions_group) {
            matrices[[region]] <- list()
            for (design in colnames(self$design)[-1]) {
                matrices[[region]][[design]] <- list()

                i <- self$design[[design]] == 1
                bam_files <- as.character(self$design[,1][i])
                matrices[[region]][[design]][["input"]] <-
                private$get_matrix(region, bam_files, bin_count)

                i <- self$design[[design]] == 2
                bam_files <- as.character(self$design[,1][i])
                matrices[[region]][[design]][["ctrl"]] <-
                private$get_matrix(region, bam_files, bin_count)
         }
      }
      matrices
    },
    get_matrix = function(region, bam_files, bcount) {
        if (length(bam_files) == 0) {
            return(NULL)
        }
        # 1. Produce the matrices
        get_matrices <- function(bf, gr) {
	  grl <- split(gr, GenomeInfoDb::seqnames(gr))
	  i <- sapply(grl, length) > 0
	  do.call("rbind", lapply(grl[i], private$get_view_means, bam_file = bf,
				  bcount = bcount))
	}
        matrices <- lapply(bam_files, get_matrices, gr = self$regions[[region]])
        # 2. Calculate the means
	means <- Reduce("+", matrices) / length(matrices)
    },
    get_view_means = function(gr, bam_file, bcount) {
	chr <- unique(as.character(GenomeInfoDb::seqnames(gr)))
        gr <- intoNbins(gr, bcount)
        stopifnot(length(chr) == 1)
        views <- Views(self$coverages[[bam_file]][[chr]], start(gr), end(gr))
	matrix(viewMeans(views), ncol = bcount, byrow = TRUE)
    },
    prepare_regions = function(regions) {
        if (class(regions) == "character") {
            names <- sapply(regions,
                            function(x) file_path_sans_ext(basename(x)))
            regions <- private$parallel_job$launch_job(data = regions,
                                                        FUN = import)
            names(regions) <- names
        } else if (class(regions) == "GRanges") {
            regions <- GRangesList(regions = regions)
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
            start(x) <- start(x) - self$params$padding_size
            start(x)[start(x) < 0] <- 1
            end(x) <- end(x) + self$params$padding_size
            # Clean seqlevels
	    x <- sortSeqlevels(x)
#            seqlevels(x) <- unique(as.character(seqnames(x)))
            x
        }))
    },
    produce_coverages = function() {
	# TODO: add support for named bam files
	regions <- GenomicRanges::reduce(BiocGenerics::unlist(self$regions))
	res <- lapply(self$params[["bam_files"]],
		      FUN = private$bam_handler$get_normalized_coverage,
		      regions = regions,
		      force_seqlevels= self$params[["force_seqlevels"]])
        names(res) <- self$params[["bam_files"]]
        lapply(res, GenomeInfoDb::sortSeqlevels)
    },
    # Produce a plot with based on a data.frame
    #
    # Input:
    #    DF:       The data frame produced by the plot.getDataFrame function
    #    title:    The title of the graph
    #    binSize:  The number of nucleotides in each bin
    #    friedman: The Friedman test result or NULL if not test done
    # Ouput:
    #    The graph that is printed on the current device.
    plot_graphic = function(DF, title, binSize, friedman=NULL) {
        # Prepare y label
        if (binSize > 1) {
            yLabel <- paste("Mean RPM for each", binSize, "positions")
        } else {
            yLabel <- paste("Mean RPM for each position")
        }
        friedmanLabel <- ifelse(is.null(friedman), "", 
                            paste0("Friedman p-value \n", 
                            signif(friedman$p.value, digits = 6), " "))
        # TODO: add x label
        p <- ggplot(DF, aes(x=position, y=value, ymin=qinf, ymax=qsup)) +
            geom_ribbon(aes(fill=group), alpha=0.3) +
            geom_line(aes(color=group), size=1, litype=1, bg="transparent") +
            theme(panel.grid.major = element_line()) +
            theme(panel.grid.minor = element_line()) +
            theme(panel.background = element_blank()) +
            theme(panel.background = element_rect()) +
            theme_bw(base_size = 20) +
            theme(axis.title.x = element_blank()) +
            ylab(yLabel) + annotate("text", label = friedmanLabel, 
                            x = Inf, y = Inf, vjust=1, hjust=1, size=4) +
            ggtitle(title)
        p
    }
  )
) 

