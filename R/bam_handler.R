# Class used to manage bam files
Bam_Handler <- R6Class("Bam_Handler",
  public = list(
    bam_files = data.frame(),
    parameters = list(),
    initialize = function(bam_files, cores = SerialParam()) {
      # Check prerequisites 
      # All BAM file names must be of string type
      if (sum(sapply(bam_files, is.character)) != length(bam_files)) {
        stop("At least one BAM file name is not a character string.")
      }

      # All BAM files must exist
      if (sum(sapply(bam_files, file.exists)) != length(bam_files)) {
        stop("At least one BAM file does not exist.")
      }
      
      # Initialize the Bam_Handler object
      private$parallel_job <- metagene:::Parallel_Job$new(cores)
      self$parameters[["cores"]] <- cores
      self$bam_files <- data.frame(matrix(nrow = length(bam_files)))
      colnames(self$bam_files) <- "old"
      self$bam_files[["old"]] <- bam_files
      self$bam_files[["bam"]] <- sapply(bam_files, private$index_bam_file)
      self$bam_files[["aligned_count"]] <-
        sapply(self$bam_files[["bam"]], private$get_file_count)
    },
    get_aligned_count = function(bam_file) {
      # Check prerequisites
      # The bam file must be in the list of bam files used for initialization
      if (! bam_file %in% self$bam_files[["old"]]) {
        stop(paste0("Bam file ", bam_file, " not found."))
      }
      i <- self$bam_files[["old"]] == bam_file
      self$bam_files[["aligned_count"]][i]
    },
    get_rpm_coefficient = function(bam_file) {
      return(self$get_aligned_count(bam_file) / 1000000)
    }
  ),
  private = list(
      parallel_job = '',
      index_bam_file = function(bam_file) {
      if (file.exists(paste(bam_file, ".bai", sep=""))  == FALSE) {
        # If there is no index file, we sort and index the current bam file
        # TODO: we need to check if the sorted file was previously produced
        #       before doing the costly sort operation
        sorted_bam_file <- paste0(basename(bam_files), ".sorted")
        sortBam(bam_file, sorted_bam_file)
        sorted_bam_file <- paste0(sorted_bam_file, ".bam")
        indexBam(sorted_bam_file)
        bam_file <- sorted_bam_file
      }
      bam_file
    },
    get_file_count = function(bam_file) {
      param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery=FALSE))
      # To speed up analysis we split the file by chromosome
      cores <- private$parallel_job$get_core_count()
      chr <- scanBamHeader(bam_file)[[bam_file]]$targets
      chr <- GRanges(seqnames = names(chr), IRanges(1, chr))
      do.call(sum, private$parallel_job$launch_job(
                            data = suppressWarnings(split(chr, 1:cores)),
                            FUN = function(x) {
                              param = ScanBamParam(which = x);
                              countBam(bam_file, param = param)$records;
                            }))

    }
  )
)
