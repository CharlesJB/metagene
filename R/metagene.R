metagene <- R6Class("metagene",
  public = list(
    # Public members
    params = list(),
    bam_files = data.frame(),
    regions = GRangesList(),
    coverages = list(),
    views = NULL,
    matrices = list(),
    
    # Methods
    initialize = function(regions, bam_files, padding_size = 2000, BPPARAM = NULL,
                          verbose = FALSE) {
      # Check params...
      
      # Save params
      self$params$padding_size <- padding_size
      self$params$verbose <- verbose
      self$params$bpparam <- BPPARAM
      
      # Prepare bam files
      cat("Prepare bam files...\n")
      self$bam_files <- private$prepare_bam_files(bam_files)
      
      # Prepare regions
      cat("Prepare regions...\n")
      self$regions <- private$prepare_regions(regions)
      
      # Parse bam files
      cat("Parse bam files...\n")
      cat("  coverages...\n")
      self$coverages <- private$produce_coverages()
      
      # Get matrix
      #self$matrices <- private$produce_matrices() # We should do this in plot
    },
    print = function(...) {
      print_values <- function(max, values) {
        max_value <- min(max, length(values))
        for (i in 1:max_value) {
          cat(paste0("   ", values[i], "\n"))
        }
        if (length(values) - max_value > 1) {
          cat("     ...\n")
        }
        if (length(values) - max_value > 0 ) {
          cat(paste0("   ", values[length(values)], "\n"))
        }
      }
      cat("metagene object:\n")
      bam_files_count <- nrow(self$bam_files)
      regions_count <- length(self$regions)
      cat("\n")
      cat(paste0(bam_files_count*regions_count, " groups:"))
      cat("\n")
      cat(paste0(" ", bam_files_count, " bam files:\n"))
      print_values(8, self$bam_files$bam)
      cat("\n")
      cat(paste0(" ", regions_count, " regions:\n"))
      print_values(8, names(self$regions))
      cat("\n")
    },
    get_bam_count = function(filename) {
      self$bam_files$alignedCount[self$bam_files$bam == filename]
    },
    produce_matrices = function() {
      # Initialize the matrices list
      ncol <- median(as.numeric(sapply(self$coverages, function(x)
          unlist(sapply(x, function(y) sapply(y, length))))
        ))
      matrices <- list()
      for (bam_file in self$bam_files$bam) {
        matrices[[bam_file]] <- list()
        for (region in names(self$regions)) {
          nrow = length(self$regions[[region]])
          current_matrix <- private$parallel_job(
              data = self$coverages[[bam_file]][[region]],
              FUN = function(x) metagene:::scaleVector(as.numeric(x), ncol)
            )
          current_matrix <- do.call("rbind", current_matrix)
          colnames(current_matrix) <- as.character(1:ncol)
          rownames(current_matrix) <- mcols(self$regions[[region]])[["id"]]
          matrices[[bam_file]][[region]] <- current_matrix
          #matrices[[bam_file]][[region]] <- matrix(nrow = nrow, ncol = ncol)
        }
      }
      matrices
    }
  ),
  private = list(
    parallel_job = function(data, FUN, ...) {
      if (!is.null(self$params$bpparam)) {
        library(BiocParallel)
        BiocParallel:::bplapply(data, FUN, BPPARAM=self$params$bpparam, ...)
      } else {
        lapply(data, FUN, ...)
      }
    },
    prepare_bam_files = function(bamFiles) {
      # Check prerequisites
      
      # All BAM file names must be of string type
      if (sum(unlist(lapply(bamFiles, is.character))) != length(bamFiles)) {
        stop("At least one BAM file name is not a valid name (a character string).")
      }
      
      # All BAM files must exist
      if (sum(unlist(lapply(bamFiles, file.exists))) != length(bamFiles)) {
        stop("At least one BAM file does not exist.")
      }
      
      # This function will only index a file if there is no index file
      indexBamFiles <- function(bamFile) {
        if (file.exists(paste(bamFile, ".bai", sep=""))  == FALSE) {
          # If there is no index file, we sort and index the current bam file
          # TODO: we need to check if the sorted file was previously produced before
          #       doing the costly sort operation
          sortedBamFile <- basename(bamFile)
          sortedBamFile <- paste(sortedBamFile, ".sorted", sep="")
          sortBam(bamFile, sortedBamFile)
          sortedBamFile <- paste(sortedBamFile, ".bam", sep="")
          indexBam(sortedBamFile)
          bamFile <- sortedBamFile
          cat(bamFile)
        }
        return(bamFile)
      }
      results <- as.data.frame(unlist(lapply(bamFiles, indexBamFiles)))
      colnames(results) <- "bam"
      results$oldBam <- bamFiles
      results$bam <- as.character(results$bam)
      results$oldBam <- as.character(results$oldBam)
      
      # This function will calculate the number of aligned read in each bam file
      countAlignedReads <- function(bamFile) {
        return(countBam(bamFile, param=ScanBamParam(flag = scanBamFlag(isUnmappedQuery=FALSE)))$records)
      }
     
      #results$alignedCount <- unlist(lapply(bamFiles, countAlignedReads))
      results$alignedCount <- unlist(private$parallel_job(
        data = bamFiles, FUN = countAlignedReads))

      return(results)
    },
    prepare_regions = function(regions) {
      if (class(regions) == "character") {
        regions <-
          metagene:::prepareRegions(regions)
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
        end(x) <- end(x) + self$params$padding_size
        # Clean seqlevels
        seqlevels(x) <- unique(as.character(seqnames(x)))
        x
      }))
    },
    produce_coverages = function() {
      coverages <- list()
      get_coverage <- function(regions,bam_file) {
        stopifnot(class(regions) == "GRanges")
        # Get a list of SimpleRleList, one element per chromosome
        
        count <- self$get_bam_count(bam_file)
        coverage <- private$parallel_job(
          # TODO: We should split by worker count instead of chr
          data = GenomicRanges::split(regions, seqnames(regions)),
          FUN = metagene:::extract_coverage_by_regions,
          bam_file = bam_file, count = count
        )
        # Merge all the SimpleRleList into a single one
        do.call("c", unname(coverage))
      }
      coverages <- lapply(self$bam_files$bam, function(bam_file) {
        cat(paste0("    ", bam_file, "...\n"))
        sapply(names(self$regions), function(name) {
          current_regions <- self$regions[[which(names(self$regions) == name)]]
          cat(paste0("      ", name, "...\n"))
          get_coverage(current_regions, bam_file)
        })
      })
      names(coverages) <- self$bam_files$bam
      coverages
    }
  )
) 
