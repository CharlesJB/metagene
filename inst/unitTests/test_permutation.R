## Test functions present in the permutation.R file

### {{{ --- Test setup ---

if(FALSE) {
  library( "RUnit" )
  library( "metagene" )
}

### }}}


###################################################
## Test the permutation_test() function
###################################################

mg <- get_demo_metagene()
mg$produce_matrices()
m1 <- mg$get_matrices()$list1$align1_rep1$input
m2 <- mg$get_matrices()$list1$align2_rep1$input

## Valid use
test.permutation_test_valid_use <- function() {
    sample_size <- min(nrow(m1), nrow(m2))
    sample_count <- 100
    FUN = function(a, b) { mean(a) - mean(b) }
    res <- permutation_test(m1, m2, sample_size = sample_size,
                            sample_count = sample_count, FUN = FUN)
    checkTrue(is.numeric(res))
    checkTrue(length(res) ==  sample_count)
}

## Valid sample_size max value
test.permutation_test_valid_sample_size_max_value <- function() {
    sample_size <- (nrow(m1) + nrow(m2)) / 2
    sample_count <- 100
    FUN = function(a, b) { mean(a) - mean(b) }
    res <- permutation_test(m1, m2, sample_size = sample_size,
                            sample_count = sample_count, FUN = FUN)
    checkTrue(is.numeric(res))
    checkTrue(length(res) ==  sample_count)
}

## Invalid m1 class
test.permutation_test_invalid_m1_class <- function() {
    obs <- tryCatch(permutation_test(1), error=conditionMessage)
    exp <- "is.matrix(matrix1) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid m2 class
test.permutation_test_invalid_m2_class <- function() {
    obs <- tryCatch(permutation_test(m1, 1), error=conditionMessage)
    exp <- "is.matrix(matrix2) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid dim m1 not equal dim m2
test.permutation_test_invalid_dim_m1_not_equal_dim_m2 <- function() {
    obs <- tryCatch(permutation_test(m1, matrix()), error=conditionMessage)
    exp <- "ncol(matrix1) == ncol(matrix2) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid sample_size class
test.permutation_test_invalid_sample_size_class <- function() {
    obs <- tryCatch(permutation_test(m1, m2, sample_size = "a"),
                    error=conditionMessage)
    exp <- "is.numeric(sample_size) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid sample_size value too large
test.permutation_test_invalid_sample_size_value_too_large <- function() {
    sample_size <- (nrow(m1) + nrow(m2)) / 2 + 1
    obs <- tryCatch(permutation_test(m1, m2, sample_size = sample_size,
                                     sample_count = 10),
                    error=conditionMessage)
    exp <- "sample_size <= ((nrow(matrix1) + nrow(matrix2))/2) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid sample_count class
test.permutation_test_invalid_sample_count_class <- function() {
    obs <- tryCatch(permutation_test(m1, m2, sample_size = 100,
                                     sample_count = "a"),
                    error=conditionMessage)
    exp <- "is.numeric(sample_count) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid fun class
test.permutation_test_invalid_fun_class <- function() {
    obs <- tryCatch(permutation_test(m1, m2, sample_size = 10,
                                     sample_count = 100, FUN = "a"),
                    error=conditionMessage)
    exp <- "is.function(FUN) is not TRUE"
    checkIdentical(obs, exp)

}
