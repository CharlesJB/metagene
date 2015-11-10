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

## Invalid data class
test.basic_stat_initialize_invalid_data_class <- function() {
    obs <- tryCatch(metagene:::Basic_Stat$new(data = c(1,2)),
                    error = conditionMessage)
    exp <- "data must be a matrix with at least one value"
    checkIdentical(obs, exp)
}

## Invalid data dimensions
test.basic_stat_initialize_invalid_data_dimension <- function() {
    obs <- tryCatch(metagene:::Basic_Stat$new(data = matrix()),
                    error = conditionMessage)
    exp <- "data must be a matrix with at least one value"
    checkIdentical(obs, exp)
}

## Invalid alpha class
test.basic_stat_initialize_invalid_data_dimension <- function() {
    obs <- tryCatch(metagene:::Basic_Stat$new(data = matrix()),
                    error = conditionMessage)
    exp <- "data must be a matrix with at least one value"
    checkIdentical(obs, exp)
}

## Invalid average value
test.basic_stat_initialize_invalid_average_value <- function() {
    obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data,
                                              average = "abc"),
                    error = conditionMessage)
    exp <- "average parameter must be either 'mean' or 'median'"
    checkIdentical(obs, exp)
}

## Invalid range class
test.basic_stat_initialize_invalid_range_class <- function() {
    obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data, range = "abc"),
                    error = conditionMessage)
    exp <- "range parameter must be a numeric of length 2"
    checkIdentical(obs, exp)
}

## Invalid range length
test.basic_stat_initialize_invalid_range_length <- function() {
    obs <- tryCatch(metagene:::Basic_Stat$new(data = valid_data, range = 1),
                    error = conditionMessage)
    exp <- "range parameter must be a numeric of length 2"
    checkIdentical(obs, exp)
}

###################################################
## Test the Basic_Stats$get_statistics() function
###################################################

## Valid case
test.basic_stat_get_statistics_valid_case <- function() {
    basic_stat <- metagene:::Basic_Stat$new(data = valid_data)
    obs <- basic_stat$get_statistics()
    checkEquals(class(obs), "data.frame")
    checkIdentical(as.numeric(dim(obs)), c(5, 4))
    checkIdentical(colnames(obs), c("position", "value", "qinf", "qsup"))
}

###################################################
## Test the Bootstrap_Stats$new() function (initialize)
###################################################

## Negative sample size
test.bootstrap_stat_get_statistics_negative_sample_size <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data,
                                                  sample_size = -1),
                  error = conditionMessage)
    exp <- "sample_size must be a positive integer."
    checkIdentical(obs, exp)
}

## Zero sample size
test.bootstrap_stat_get_statistics_zero_sample_size <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data,
                                                  sample_size = 0),
                    error = conditionMessage)
    exp <- "sample_size must be a positive integer."
    checkIdentical(obs, exp)
}

## Not numeric sample size
test.bootstrap_stat_get_statistics_non_numeric_sample_size <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data,
                                                  sample_size = "1"),
                    error = conditionMessage)
    exp <- "sample_size must be a positive integer."
    checkIdentical(obs, exp)
}

## Decimal sample size
test.bootstrap_stat_get_statistics_decimal_sample_size <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data,
                                                  sample_size = 1.1),
                    error = conditionMessage)
    exp <- "sample_size must be a positive integer."
    checkIdentical(obs, exp)
}

## Negative sample count
test.bootstrap_stat_get_statistics_negative_sample_count <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data,
                                                  sample_count = -1),
                    error = conditionMessage)
    exp <- "sample_count must be a positive integer."
    checkIdentical(obs, exp)
}

## Zero sample count
test.bootstrap_stat_get_statistics_zero_sample_count <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data,
                                                  sample_count = 0),
                    error = conditionMessage)
    exp <- "sample_count must be a positive integer."
    checkIdentical(obs, exp)
}

## Not numeric sample count
test.bootstrap_stat_get_statistics_non_numeric_sample_count <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data,
                                                  sample_count = "1"),
                    error = conditionMessage)
    exp <- "sample_count must be a positive integer."
    checkIdentical(obs, exp)
}

## Decimal sample count
test.bootstrap_stat_get_statistics_decimal_sample_count <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data,
                                                  sample_count = 1.1),
                    error = conditionMessage)
    exp <- "sample_count must be a positive integer."
    checkIdentical(obs, exp)
}

## Not boolean debug
test.bootstrap_stat_initialize_not_boolean_debug <- function() {
    obs <- tryCatch(metagene:::Bootstrap_Stat$new(data = valid_data, debug = 1),
                    error = conditionMessage)
    exp <- "debug must be TRUE or FALSE."
    checkIdentical(obs, exp)
}

###################################################
## Test the Bootstrap_Stats$get_statistics() function
###################################################

## Valid case
test.bootstrap_stat_get_statistics_valid_case <- function() {
    bootstrap_stat <- metagene:::Bootstrap_Stat$new(data = valid_data)
    obs <- bootstrap_stat$get_statistics()
    checkEquals(class(obs), "data.frame")
    checkIdentical(as.numeric(dim(obs)), c(5, 4))
    checkIdentical(colnames(obs), c("position", "value", "qinf", "qsup"))
}

## Valid case - debug
test.bootstrap_stat_get_statistics_valid_case_debug <- function() {
    bootstrap_stat <- metagene:::Bootstrap_Stat$new(data = valid_data,
                                                    debug = TRUE)
    obs <- bootstrap_stat$get_statistics()
    checkEquals(class(obs), "list")
    checkEquals(length(obs), 3)
    checkEquals(class(obs$statistics), "data.frame")
    checkIdentical(as.numeric(dim(obs$statistics)), c(5, 4))
    checkIdentical(colnames(obs$statistics),
                   c("position", "value", "qinf", "qsup"))
    checkIdentical(class(obs$values), "list")
    checkEquals(length(obs$values), 5)
    checkIdentical(unique(sapply(obs$values, class)), "numeric")
    checkEquals(unique(sapply(obs$values, length)), 1000)
    checkIdentical(class(obs$replicates), "list")
    checkEquals(length(obs$replicates), 5)
    checkIdentical(unique(sapply(obs$replicates, class)), "matrix")
    checkEquals(unique(sapply(obs$replicates, nrow)), 20)
    checkEquals(unique(sapply(obs$replicates, ncol)), 1000)
}
