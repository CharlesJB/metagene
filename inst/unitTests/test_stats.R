## Test functions present in the stats.R file

### {{{ --- Test setup ---

if(FALSE) {
  library( "RUnit" )
  library( "metagene" )
}

### }}}

valid_data <- matrix(1:100, ncol = 5)

###################################################
## Test the Basic_Stats$new() function (initialize)
###################################################

base_msg <- "Basic_Stat initialize -"

## Invalid data class
test.basic_stat_initialize_invalid_data_class <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = c(1,2)), error=conditionMessage)
  exp <- "data must be a matrix with at least one value"
  msg <- paste(base_msg, "An invalid data class did not generate an exception with expected message." )
  checkIdentical(obs, exp, msg)
}

## Invalid data dimensions
test.basic_stat_initialize_invalid_data_dimension <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = matrix()), error=conditionMessage)
  exp <- "data must be a matrix with at least one value"
  msg <- paste(base_msg, "Invalid data dimensions did not generate an exception with expected message." )
  checkIdentical(obs, exp, msg)
}

## Invalid alpha class
test.basic_stat_initialize_invalid_data_dimension <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = matrix()), error=conditionMessage)
  exp <- "data must be a matrix with at least one value"
  msg <- paste(base_msg, "Invalid data dimensions did not generate an exception with expected message." )
  checkIdentical(obs, exp, msg)
}

## Invalid average value
test.basic_stat_initialize_invalid_average_value <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data, average = "abc"), error=conditionMessage)
  exp <- "average parameter must be either 'mean' or 'median'"
  msg <- paste(base_msg, "Invalid average value did not generate an exception with expected message." )
  checkIdentical(obs, exp, msg)
}

## Invalid range class
test.basic_stat_initialize_invalid_range_class <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data, range = "abc"), error=conditionMessage)
  exp <- "range parameter must be a numeric of length 2"
  msg <- paste(base_msg, "Invalid range class did not generate an exception with expected message." )
  checkIdentical(obs, exp, msg)
}

## Invalid range length
test.basic_stat_initialize_invalid_range_length <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data, range = 1), error=conditionMessage)
  exp <- "range parameter must be a numeric of length 2"
  msg <- paste(base_msg, "Invalid range class did not generate an exception with expected message." )
  checkIdentical(obs, exp, msg)
}

## Invalid cores numeric values - zero
test.basic_stat_initialize_zero_core_number <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data, cores = 0), error=conditionMessage)
  exp <- "cores must be positive numeric or BiocParallelParam instance."
  msg <- paste(base_msg, "A zero core number argument did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Invalid cores numeric values - negative value
test.basic_stat_initialize_negative_core_number <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data, cores = -1), error=conditionMessage)
  exp <- "cores must be positive numeric or BiocParallelParam instance."
  msg <- paste(base_msg, "A negative core number argument did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Invalid cores non-numeric values - string
test.basic_stat_initialize_string_core_number <- function() {
  obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data, cores = "1"), error=conditionMessage)
  exp <- "cores must be positive numeric or BiocParallelParam instance."
  msg <- paste(base_msg, "A string core number argument did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

###################################################
## Test the Basic_Stats$get_statistics() function
###################################################

base_msg <- "Basic_Stat get_statistics -"

## Valid case
test.basic_stat_get_statistics_valid_case <- function() {
  basic_stat <- metagene:::Basic_Stat$new(data = valid_data)
  obs <- basic_stat$get_statistics()
  msg <- paste(base_msg, "Statistics do not have correct class.")
  checkEquals(class(obs), "data.frame", msg)
  msg <- paste(base_msg, "Statistics do not have the correct dimensions.")
  checkIdentical(as.numeric(dim(obs)), c(5, 4), msg)
  msg <- paste(base_msg, "Statistics do not have the right colnames.")
  checkIdentical(colnames(obs), c("position", "value", "qinf", "qsup"), msg)
}

###################################################
## Test the Bootstrap_Stats$new() function (initialize)
###################################################

base_msg <- "Bootstrap_Stat  (initialize) -"

## Negative sample size
test.bootstrap_stat_get_statistics_negative_sample_size <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, sample_size = -1),
                  error=conditionMessage)
  exp <- "sample_size must be a positive integer."
  msg <- paste(base_msg, "A negative sample size did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Zero sample size
test.bootstrap_stat_get_statistics_zero_sample_size <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, sample_size = 0),
                  error=conditionMessage)
  exp <- "sample_size must be a positive integer."
  msg <- paste(base_msg, "A negative sample size did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Not numeric sample size
test.bootstrap_stat_get_statistics_non_numeric_sample_size <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, sample_size = "1"),
                  error=conditionMessage)
  exp <- "sample_size must be a positive integer."
  msg <- paste(base_msg, "A negative sample size did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Decimal sample size
test.bootstrap_stat_get_statistics_decimal_sample_size <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, sample_size = 1.1),
                  error=conditionMessage)
  exp <- "sample_size must be a positive integer."
  msg <- paste(base_msg, "A negative sample size did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Negative sample count
test.bootstrap_stat_get_statistics_negative_sample_count <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, sample_count = -1),
                  error=conditionMessage)
  exp <- "sample_count must be a positive integer."
  msg <- paste(base_msg, "A negative sample count did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Zero sample count
test.bootstrap_stat_get_statistics_zero_sample_count <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, sample_count = 0),
                  error=conditionMessage)
  exp <- "sample_count must be a positive integer."
  msg <- paste(base_msg, "A negative sample count did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Not numeric sample count
test.bootstrap_stat_get_statistics_non_numeric_sample_count <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, sample_count = "1"),
                  error=conditionMessage)
  exp <- "sample_count must be a positive integer."
  msg <- paste(base_msg, "A negative sample count did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Decimal sample count
test.bootstrap_stat_get_statistics_decimal_sample_count <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, sample_count = 1.1),
                  error=conditionMessage)
  exp <- "sample_count must be a positive integer."
  msg <- paste(base_msg, "A negative sample count did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

## Not boolean debug
test.bootstrap_stat_initialize_not_boolean_debug <- function() {
  obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, debug = 1),
                  error=conditionMessage)
  exp <- "debug must be TRUE or FALSE."
  msg <- paste(base_msg, "A not boolean debug did not generate an exception with expected message.")
  checkIdentical(obs, exp, msg)
}

###################################################
## Test the Bootstrap_Stats$get_statistics() function
###################################################

base_msg <- "Bootstrap_Stat get_statistics -"

## Valid case
test.bootstrap_stat_get_statistics_valid_case <- function() {
  bootstrap_stat <- metagene:::Bootstrap_Stat$new(data = valid_data)
  obs <- bootstrap_stat$get_statistics()
  msg <- paste(base_msg, "Statistics do not have correct class.")
  checkEquals(class(obs), "data.frame", msg)
  msg <- paste(base_msg, "Statistics do not have the correct dimensions.")
  checkIdentical(as.numeric(dim(obs)), c(5, 4), msg)
  msg <- paste(base_msg, "Statistics do not have the right colnames.")
  checkIdentical(colnames(obs), c("position", "value", "qinf", "qsup"), msg)
}

## Valid case - debug

## Valid case
test.bootstrap_stat_get_statistics_valid_case_debug <- function() {
  bootstrap_stat <- metagene:::Bootstrap_Stat$new(data = valid_data, debug = TRUE)
  obs <- bootstrap_stat$get_statistics()
  msg <- paste(base_msg, "Debug result do not have correct class.")
  checkEquals(class(obs), "list", msg)
  msg <- paste(base_msg, "Debug result do not have the correct length.")
  checkEquals(length(obs), 3, msg)
  msg <- paste(base_msg, "Statistics do not have the correct class")
  checkEquals(class(obs$statistics), "data.frame", msg)
  msg <- paste(base_msg, "Statistics do not have the correct dimensions.")
  checkIdentical(as.numeric(dim(obs$statistics)), c(5, 4), msg)
  msg <- paste(base_msg, "Statistics do not have the right colnames.")
  checkIdentical(colnames(obs$statistics), c("position", "value", "qinf", "qsup"), msg)
  msg <- paste(base_msg, "Values list do not have the right class")
  checkIdentical(class(obs$values), "list", msg)
  msg <- paste(base_msg, "Values list do not have the right length")
  checkEquals(length(obs$values), 5, msg)
  msg <- paste(base_msg, "Values do not have the right class")
  checkIdentical(unique(sapply(obs$values, class)), "numeric", msg)
  msg <- paste(base_msg, "Values do not have the right length")
  checkEquals(unique(sapply(obs$values, length)), 1000, msg)
  msg <- paste(base_msg, "Replicates list do not have the right class")
  checkIdentical(class(obs$replicates), "list", msg)
  msg <- paste(base_msg, "Replicates list do not have the right length")
  checkEquals(length(obs$replicates), 5, msg)
  msg <- paste(base_msg, "Replicates do not have the right class")
  checkIdentical(unique(sapply(obs$replicates, class)), "matrix", msg)
  msg <- paste(base_msg, "Replicates do not have the right number of rows")
  checkEquals(unique(sapply(obs$replicates, nrow)), 20, msg)
  msg <- paste(base_msg, "Replicates do not have the right number of columns")
  checkEquals(unique(sapply(obs$replicates, ncol)), 1000, msg)
}
