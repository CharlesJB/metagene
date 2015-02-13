#' A class to manage metagene analysis.
#'
#' This class will allow to load, convert and normalize alignments and regions
#' files/data. Once the data is ready, the user can then chose to produce
#' metagene plots on the data (or a subset of the data).
#'
#' @section Constructor:
#' \describe{
#'   \item{}{\code{mg <- metagene$new(regions, bam_files, padding_size = 0,
#'                              cores = SerialParam(), verbose = FALSE)}}
#'   \item{regions}{Either a \code{list} of bed filenames, a \code{GRanges}
#'                  object or a \code{GRangesList} object.}
#'   \item{bam_files}{A list of bam filenames. The bam files must be indexed.
#'                    i.e.: if a file is named file.bam, there must be a file
#'                    named file.bam.bai in the same directory.}
#'   \item{padding_size}{The regions will be extended on each side by the value
#'                       of this parameter. The padding_size must be a 
#'                       non-negative integer. Default = 0.}
#'   \item{cores}{The number of cores available to parallelize the analysis.
#'                Either a positive integer or a \code{BiocParallelParam}.
#'                Default: SerialParam().}
#'   \item{verbose}{Print progression of the analysis. A logical constant. 
#'                  Default: \code{FALSE}.}
#' }
#'
#'  \code{metagene$new} returns a \code{metagene} object that contains the
#'  coverages for every bam files in the regions from the \code{regions} param.
#'
#' @return
#' \code{metagene$new} returns a \code{metagene} object which contains the
#' normalized coverage values for every regions and for every BAM files.
#'
#' @section Methods:
#' \describe{
#'   \item{}{\code{df <- mg$plot(design = NULL, regions_group = NULL,
#'               bin_size = 100, alpha = 0.05, sample_count = 1000)}}
#'   \item{design}{A \code{data.frame} that describe to experiment to plot. The
#'                 first column must be the filenames. The other columns (at
#'                 least one is required) represent how the files should be
#'                 grouped. All the files in the same group will be combined
#'                 before doing the statistical analysis. For each column,
#'                 there will be one line with its ribon in the metagene plot.
#'                 Default = NULL.
#'
#'                 0: Do not use file.
#'                 1: File is input.
#'                 2: File is control.}
#'   \item{region_group}{A list of region names to include in the analysis. If
#'                      \code{NULL}, all the regions used when creating the
#'                      \code{metagene} object will be used.
#'                      Default: \code{NULL}}
#'   \item{bin_size}{The size of bin to use before calculating the statistics.
#'                   larger \code{bin_size} will reduce the calculation time
#'                   and produce smoother curves. Default = 100.}
#'   \item{alpha}{The condifence interval (CI) to represent with a ribbon. Must
#'                be a value between 0 and 1. With a value of 0.05, the ribbon
#'                will represent 95% of the estimated values of the means.
#'                Default = 0.05.}
#'   \item{sample_count}{The number of draw to do during bootstrap analysis.
#'                       Default: 1000.}
#' }
#' \describe{
#'   \item{}{\code{mg$export(bam_file, region, file}}
#'   \item{bam_file}{The name of the bam file to export.}
#'   \item{region}{The name of the region to export.}
#'   \item{file}{The name of the ouput file.}
#' }
#' \describe{
#'   \item{}{\code{mg$heatmap(region, bam_file, bin_size}}
#'   \item{region}{The name of the region to export.}
#'   \item{bam_file}{The name of the bam file to export.}
#'   \item{bin_size}{The size of the bin to produce before creating heatmap.}
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
                          cores = SerialParam(), verbose = FALSE) {
      # Check params...
      private$check_param(regions = regions, bam_files = bam_files, 
                          padding_size = padding_size,
                          cores = cores, verbose = verbose)
      # Save params
      private$parallel_job <- Parallel_Job$new(cores)
      self$params[["padding_size"]] <- padding_size
      self$params[["verbose"]] <- verbose
      self$params[["bam_files"]] <- bam_files
      
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
      private$bam_handler$get_aligned_count(filename)
    },
    plot = function(design = NULL, regions_group = NULL, bin_size = 100,
                    alpha = 0.05, sample_count = 1000) {
      self$design <- private$prepare_design(design)

      # 1. Get the correctly formatted matrices
      self$matrices <- private$produce_matrices(regions_group, bin_size)

      # 2. Calculate means and confidence intervals
      sample_size = as.integer(min(sapply(self$matrices, sapply, sapply, nrow)))
      DF <- data.frame(group = character(), position = numeric(),
                       value = numeric(), qinf = numeric(), qsup = numeric())
      for (region in names(self$matrices)) {
        for (bam_file in names(self$matrices[[region]])) {
          group_name <- paste(bam_file, region, sep = "_")
          print(group_name)
          data <- self$matrices[[region]][[bam_file]][["input"]]
          ctrl <- self$matrices[[region]][[bam_file]][["ctrl"]]
          bootstrap_stat <- Bootstrap_Stat$new(data = data, ctrl = ctrl,
                                               sample_size = sample_size)
          current_DF <- bootstrap_stat$get_statistics()
          current_DF <- cbind(rep(group_name, nrow(current_DF)), current_DF)
          colnames(current_DF)[1] <- "group"
          DF <- rbind(DF, current_DF)
        }
      }
      # 3. Produce the graph
      #    DF <- metagene:::getDataFrame(bootstrap_result, 
      #                                     range=c(-1,1), binSize=1)
      private$plot_graphic(DF, paste(unique(DF[["group"]]), collapse=" vs "), 
                            binSize = 1)
      return(DF)
    },
    export = function(bam_file, region, file) {
      region <- tools::file_path_sans_ext(basename(region))
      region <- self$regions[[region]]
      param <- Rsamtools::ScanBamParam(which = region)
      alignments <- GenomicAlignments::readGAlignments(bam_file, param = param)
      weight <- 1 - private$bam_handler$get_rpm_coefficient(bam_file)
      seqlevels(alignments) <- seqlevels(region)
      coverage <- GenomicAlignments::coverage(alignments, weight=weight)
      rtracklayer::export(coverage, file, "BED")
      invisible(coverage)
    },
    heatmap = function(region, bam_file, bin_size = 10) {
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
                           cores, verbose) {
        # Check parameters validity
        if (!is.logical(verbose)) {
            stop("verbose must be a logicial value (TRUE or FALSE)")
        }
        if (!(is.numeric(padding_size) || is.integer(padding_size)) || 
            padding_size < 0 || as.integer(padding_size) != padding_size) {
            stop("padding_size must be a non-negative integer")
        }
        isBiocParallel = is(cores, "BiocParallelParam")
        isInteger = ((is.numeric(cores) || is.integer(cores)) && 
                         padding_size > 0 &&  as.integer(cores) == cores)
        if (!isBiocParallel && !isInteger) {
            stop(paste0("cores must be a positive numeric or ", 
                        "BiocParallelParam instance"))
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
    produce_matrices = function(regions_group, bin_size) {
      matrices <- list()
      if (is.null(regions_group)) {
        regions_group <- names(self$regions)
      }
      for (region in regions_group) {
        matrices[[region]] <- list()
         for (design in colnames(self$design)[-1]) {
           matrices[[region]][[design]] <- list()

           i <- self$design[[design]] == 1
           bam_files <- self$design[,1][i]
           matrices[[region]][[design]][["input"]] <-
             private$get_matrix(region, bam_files, bin_size)

           i <- self$design[[design]] == 2
           bam_files <- self$design[,1][i]
           matrices[[region]][[design]][["ctrl"]] <-
             private$get_matrix(region, bam_files, bin_size)
         }
      }
      matrices
    },
    get_matrix = function(region, bam_files, bin_size) {
      if (length(bam_files) == 0) {
        return(NULL)
      }
      res <- do.call(rbind, lapply(bam_files, function(bam_file)
        t(sapply(self$coverages[[bam_file]][[region]], as.numeric))
      ))
      private$bin_matrix(res, bin_size)
    },
    prepare_regions = function(regions) {
      if (class(regions) == "character") {
        names <- sapply(regions, function(x) file_path_sans_ext(basename(x)))
        regions <- private$parallel_job$launch_job(data = regions,
                                                        FUN = import)
        names(regions) <- names
      } else if (class(regions) == "GRanges") {
        regions <- GRangesList(regions = regions)
      } else if (class(regions) == "GRangesList") {
        regions <- regions
      }
      # TODO: Check if there is a id column in the mcols of every ranges.
      #       If not, add one by merging seqnames, start and end.
      GRangesList(lapply(regions, function(x) {
        # Add padding
        start(x) <- start(x) - self$params$padding_size
	start(x)[start(x) < 0] <- 1
        end(x) <- end(x) + self$params$padding_size
        # Clean seqlevels
        seqlevels(x) <- unique(as.character(seqnames(x)))
        x
      }))
    },
    produce_coverages = function() {
      res <- lapply(self$params[["bam_files"]], function(bam_file) {
        lapply(self$regions, function(regions) {
          private$bam_handler$get_normalized_coverage(bam_file, regions)
          })
      })
      names(res) <- self$params[["bam_files"]]
      res
    },
    # Bin matrix columns
    #
    # Input:
    #    data:        The matrix to bin.
    #    bin_size:    The number of nucleotides in each bin.
    #
    # OUPUT:
    #    A matrix with each column representing the mean of 
    #    bin_size nucleotides.
    bin_matrix = function(data, bin_size) {
      stopifnot(bin_size %% 1 == 0)
      stopifnot(bin_size > 0)
      
      # If the number of bins is too small, we warn the user.
      bin_count <- floor(ncol(data) / bin_size)
      if (bin_count < 5) {
        message <- paste0("Number of bins is very small: ", bin_count, "\n")
        message <- paste0(message, "  You should consider reducing the ", 
                        "bin_size value.")
        warning(message)
      }
      # If bin_size is not a multiple of ncol(data), we remove the bin that 
      # have a different size than the others.
      remainder <- ncol(data) %% bin_size
      if (remainder != 0) {
        message <- paste0("bin_size is not a multiple of the number of ", 
                        "columns in data.\n")
        message <- paste0(message, "  Columns ", ncol(data) - remainder, "-", 
                        ncol(data), " will be removed.")
        warning(message)
        data <- data[,1:(ncol(data)-remainder)]
      }
      
      splitMean <- function(x, bs) {
        if (bs < length(x)) {
          return(tapply(x, (seq_along(x)-1) %/% bs, mean))
        } else {
          return(x)
        }
      }
      return(unname(t(apply(data, 1, splitMean, bin_size))))
    },
    # Produce a plot with based on a data.frame
    #
    # Input:
    #    DF:       The data frame produced by the plot.getDataFrame function
    #    title:    The title of the graph
    #
    # Ouput:
    #    The graph that is printed on the current device.
    plot_graphic = function(DF, title, binSize) {
      # Prepare y label
      if (binSize > 1) {
        yLabel <- paste("Mean RPM for each", binSize, "positions")
      } else {
        yLabel <- paste("Mean RPM for each position")
      }
      # TODO: add x label
      p <- ggplot(DF, aes(x=position, y=value, ymin=qinf, ymax=qsup)) +
        geom_ribbon(aes(fill=group), alpha=0.3) +
        geom_line(aes(color=group),size=1,litype=1,bg="transparent")+
        theme(panel.grid.major = element_line())+
        theme(panel.grid.minor = element_line())+
        theme(panel.background = element_blank())+
        theme(panel.background = element_rect())+
        theme_bw(base_size = 20)+
        theme(axis.title.x = element_blank())+
        ylab(yLabel)+
        ggtitle(title)
      print(p)
    }
  )
) 
