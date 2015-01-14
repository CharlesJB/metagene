# Class representing a abstract Stat which is used to calculate values and
# confidence intervals of a matrix in a column by column manner.

Stat <- R6Class("Stat",
  public = list(
    initialize = function(data, alpha = 0.05, average = "mean",
                          range = c(-1, 1), cores = SerialParam()) {

      private$check_param(data = data, alpha = alpha, average = average,
                          range = range)

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
    check_param = function(data, alpha, average, range, cores) {
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
    },
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
      position <- seq(range[1], range[2], length.out = ncol(data))
      res <- data.frame(do.call(rbind, private$parallel_job$launch_job(
        data = split(data, rep(1:ncol(data), each = nrow(data))),
        FUN = calculate_statistic)))
      cbind(position, res, row.names = NULL)
    }
  )
)

# Class representing the bootstrap stat.

Bootstrap_Stat <- R6Class("Bootstrap_Stat",
  inherit = Stat,
  public = list(
    initialize = function(data, alpha = 0.05, average = "mean",
                          range = c(-1, 1), cores = SerialParam(),
                          sample_count = 1000, sample_size = NA) {

      # Check parameters validity
      super$check_param(data = data, alpha = alpha, average = average,
                        range = range)
      if (!is.numeric(sample_count)
          || as.integer(sample_count) != sample_count
          || sample_count < 1) {
        stop("sample_count must be a positive integer.")
      }
      if (!is.na(sample_size)) {
        if (!is.numeric(sample_size)
            || as.integer(sample_size) != sample_count
            || sample_size < 1) {
          stop("sample_size must be a positive integer.")
        }
      }

      # Save parameters
      private$parameters[["sample_count"]] <- sample_count
      if (is.na(sample_size)) {
        sample_size <- nrow(data)
      }
      private$parameters[["sample_size"]] <- sample_size

      # Initialize and calculate statistic
      super$initialize(data = data, alpha = alpha, average = average,
                       range = range, cores = cores)
    }
  ),
  private = list(
    calculate_statistics = function() {
      # Fetch relevant params
      data <- private$data
      range <- private$parameters[["range"]]

      # Calculate results
      position <- seq(range[1], range[2], length.out = ncol(data))
      res <- data.frame(do.call(rbind, lapply(
        split(data, rep(1:ncol(data), each = nrow(data))),
        private$calculate_statistic)))
      cbind(position, res, row.names = NULL)
    },
    # Calculate the statistic for a single column
    calculate_statistic = function(column_values) {
      # Check param
      stopifnot(is.numeric(column_values))
      stopifnot(length(column_values) > 0)

      # Fetch relevant parameters
      alpha <- private$parameters[["alpha"]]
      average <- private$parameters[["average"]]

      # Calculate result
      replicates <- private$generate_draw_values(column_values)
      values <- private$calculate_replicate_values(replicates)

      res <- quantile(values, c(0.5, alpha/2, 1-(alpha/2)))
      names(res) <- c("value", "qinf", "qsup")
      res
    },
    generate_draw_values = function(column_values) {
      sample_count <- private$parameters[["sample_count"]]
      sample_size <- private$parameters[["sample_size"]]

      sample_data <- function(data) {
        data[sample(1:length(data), sample_size, replace=TRUE)]
      }
      replicate(sample_count, sample_data(column_values))
    },
    calculate_replicate_values = function(replicates) {
      workers <- private$parallel_job$get_core_count()
      average <- private$parameters[["average"]]
      sample_size <- private$parameters[["sample_size"]]

      groups <- suppressWarnings(split(replicates, 1:workers))
      groups <- lapply(groups, function(x) matrix(x, nrow = sample_size))
      get_average <- function(group, average) {
        apply(group, 2, average)
      }
      values <- private$parallel_job$launch_job(groups, get_average, average)
      do.call(c, values)
    }
  )
)