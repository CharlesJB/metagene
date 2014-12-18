# Class representing a abstract Stat which is used to calculate values and
# confidence intervals of a matrix in a column by column manner.

Stat <- R6Class("Stat",
                public = list(
    initialize = function(data, alpha = 0.05, average = "mean",
                          range = c(-1, 1), cores = SerialParam()) {
      
      # Check parameters validity
      if (!is.matrix(data) || sum(!is.na(data) == 0)) {
        stop("data must be a matrix with at least one value")
      }
      if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
        stop("alpha parameter must be a numeric between 0 and 1")
      }
      if (! average %in% c("mean", "median")) {
        stop("average parameter must be either 'mean' or 'median'")
      }
      if (!is.numeric(range) | length(range) != 2) {
        stop("range parameter must be a numeric of length 2")
      }
      
      # Save parameters
      if (average == "mean") {
        private$parameters[["average"]] <- base:::mean
      } else if (average == "median") {
        private$parameters[["average"]] <- median
      }
      private$parameters[["alpha"]] <- alpha
      private$parameters[["cores"]] <- cores
      private$parameters[["range"]] <- range
      
      private$parallel_job <- metagene:::Parallel_Job$new(cores)
      private$data <- data
      
      # Calculate the statistics
      private$statistics <- private$calculate_statistics()
      private$validate_statistics()
    },
    get_statistics = function() {
      private$statistics
    }
  ),
  private = list(
    statistics = data.frame(position = numeric(),
                            value = numeric(),
                            qinf = numeric(),
                            qsup = numeric()),
    parameters = list(),
    data = matrix(),
    parallel_job = '',
    validate_statistics = function() {
      error <- "Stat object is not in a valid state:"
      
      # Check class
      if (class(private$statistics) != "data.frame") {
        reason = "statistics must be a data.frame."
        stop(paste(error, reason))
      }
      
      # Check columns
      expected <- c("position", "value", "qinf", "qsup")
      if (!identical(colnames(private$statistics), expected)) {
        reason <- "invalid column names."
        stop(paste(error, reason))
      }
      if (nrow(private$statistics > 0)) {
        if (!all(sapply(private$statistics, class) == "numeric")) {
          reason <- "invalid column classes."
          stop(paste(error, reason))
        }
      }
    }
  )
)

# Class representing basic Stat. The value is either the mean or the median
# (selected by the user) and qsup and qinf represent the range of the data in
# the alpha/2:1-alpha/2 percentile.

Basic_Stat <- R6Class("Basic_Stat",
  inherit = Stat,
  public = list(
  ),
  private = list(
    calculate_statistics = function() {
      # Fetch relevant params
      alpha <- private$parameters[["alpha"]]
      data <- private$data
      average <- private$parameters[["average"]]
      range <- private$parameters[["range"]]
      
      # Prepare function
      calculate_statistic <- function(column_values) {
        value <- average(column_values)
        qinf <- unname(quantile(column_values, alpha/2))
        qsup <- unname(quantile(column_values, 1-(alpha/2)))
        c(value = value, qinf = qinf, qsup = qsup)
      }
      # Calculate results
      position <- seq(range[1], range[2], length.out = nrow(data))
      res <- data.frame(do.call(rbind, private$parallel_job$launch_job(
        data = split(data, rep(1:ncol(data), each = nrow(data))),
        FUN = calculate_statistic)))
      cbind(position, res, row.names = NULL)
    }
  )
)
