# Class used to manage parallell jobs
Parallel_Job <- R6Class("Parallel_Job",
  public = list(
    parameters = list(),
    BPPARAM = SerialParam(),
    initialize = function(cores = SerialParam()) {
      self$set_core_count(cores)
    },
    launch_job = function(data, FUN, ...) {
      BiocParallel:::bplapply(data, FUN, BPPARAM = self$BPPARAM, ...)
    },
    get_core_count = function() {
      cores <- self$parameters[["cores"]]
      ifelse(is.numeric(cores), cores, as.numeric(bpworkers(cores)))
    },
    set_core_count = function(cores) {
      # Note: cores can be numeric or BiocParallelParam instance
      #       BPPARAM is always a BiocParallelParam instance
      self$parameters[["cores"]] <- cores
      self$BPPARAM <- private$get_bpparam(cores)
    }
  ),
  private = list(
    get_bpparam = function(cores) {
      if (is.numeric(cores)) {
        # The number of cores has to be a positive integer
        if(as.integer(cores) != cores || cores <= 0) {
          stop("cores must be positive numeric or BiocParallelParam instance.")
        }
        if (cores == 1) {
          BPPARAM <- SerialParam()
        } else {
          BPPARAM <- MulticoreParam(workers = cores)
        }
      } else {
        # Must be one of the BiocParallelParam class
        if (class(cores) != "SerialParam"
            & class(cores) != "MulticoreParam"
            & class(cores) != "SnowParam"
            & class(cores) != "BatchJobParam"
            & class(cores) != "DoparParam") {
          stop("cores must be positive numeric or BiocParallelParam instance.")
        }
        BPPARAM <- cores
      }
      BPPARAM
    }
  )
)
# 
# # This function take care of checking the values of the cores parameter and
# # of launching the parallel code
# parallel_job = function(data, FUN, cores = SerialParam(), ...) {
#   if (is.numeric(cores)) {
#     # The number of cores has to be a positive integer
#     if(as.integer(cores) != cores || cores <= 0) {
#       stop("The number of cores has to be a positive integer.")
#     }
#     if (cores == 1) {
#       BPPARAM <- SerialParam()
#     } else {
#       BPPARAM <- MulticoreParam(workers = cores)
#     }
#   } else {
#     # Must be one of the BiocParallelParam class
#     if (class(cores) != "SerialParam"
#         & class(cores) != "MulticoreParam"
#         & class(cores) != "SnowParam"
#         & class(cores) != "BatchJobParam"
#         & class(cores) != "DoparParam") {
#       stop("Param cores must be numeric or BiocParallelParam instance.")
#     }
#     BPPARAM <- cores
#   }
#   BiocParallel:::bplapply(data, FUN, BPPARAM=BPPARAM, ...)
# }
