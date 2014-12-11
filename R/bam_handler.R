# Class used to manage bam files
bam_handler <- R6Class("bam_handler",
  public = list(
    initialize = function(bam_files) {
      # Check prerequisites 
      # All BAM file names must be of string type
      if (sum(sapply(bam_files, is.character)) != length(bam_files)) {
        stop("At least one BAM file name is not a valid name (a character string).")
      }

      # All BAM files must exist
      if (sum(sapply(bam_files, file.exists)) != length(bam_files)) {
        stop("At least one BAM file does not exist.")
      }
      
      private$bam_files <<- data.frame(bam <- bam_files, old <- bam_files, aligned_count <- 0)
      private$bam_files[["bam"]] <- sapply(bam_files, private$index_bam_file)
      private$bam_files[["aligned_count"]] <- sapply(private$bam_files[["bam"]], private$get_file_count)
    },
    get_aligned_counts = function(bam_file) {
      # Check prerequisites
      # The bam file must be in the list of bam files used for initialization
      if (! bam_file %in% private$bam_files[["old"]]) {
        stop(paste0("Bam file ", bam_file, " not found."))
      }
      i <- private$bam_files[["old"]] == bam_file
      private$bam_files[["aligned_count"]][i]
    },
    get_rpm_coefficient = function(bam_file) {
      self$get_aligned_counts(bam_file) / 1000000
    }
  ),
  private = list(
    index_bam_file = function(bam_file) {
      if (file.exists(paste(bam_file, ".bai", sep=""))  == FALSE) {
        # If there is no index file, we sort and index the current bam file
        # TODO: we need to check if the sorted file was previously produced before
        #       doing the costly sort operation
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
        countBam(bam_file, param=param)$records
    }
  )
)
