#' Perform a permutation test on 2 matrices
#'
#' The goal of this function is to calculate the values of a test performed
#' by \code{FUN} after each of \code{sample_count} permutations.
#'
#' Each round of the permutation test, two new matrices will be randomly
#' sampled from using the combination of the two original matrices. The means
#' of each columns will be calculated to produce the vectors that will be sent
#' \code{FUN}.
#' 
#' @param matrix1 The first matrix.
#' @param matrix2 The second matrix.
#' @param sample_size The number of element to draw for each matrix.
#' @param sample_count The number of permutations.
#' @param FUN The function to use to compare the 2 matrices. First two params
#'            must be \code{numeric} vector and the return must be a single
#'            \code{numeric} value.
#' @param ... Extra param for \code{FUN}.
#'
#' @return A \code{vector} of numeric corresponding to the result of each
#'         permutation.
#'
#' @examples
#' # Get some matrices
#' mg <- get_demo_metagene()
#' mg$produce_matrices()
#' m1 <- mg$get_matrices()$list1$align1_rep1$input
#' m2 <- mg$get_matrices()$list1$align2_rep1$input
#'
#' # Perform permutation test
#' sample_size <- min(nrow(m1), nrow(m2))
#' FUN = function(a, b) { mean(a) - mean(b) } # Dummy function for demo purpose
#' permutation_results <- permutation_test(m1, m2, sample_size = sample_size,
#'                                         sample_count = 1000, FUN = FUN)
#'
#' @export
permutation_test <- function(matrix1, matrix2, sample_size, sample_count, FUN,
                             ...) {
    stopifnot(is.matrix(matrix1))
    stopifnot(is.matrix(matrix2))
    stopifnot(ncol(matrix1) == ncol(matrix2))
    stopifnot(is.numeric(sample_size))
    stopifnot(is.numeric(sample_count))
    stopifnot(sample_size <= ((nrow(matrix1) + nrow(matrix2)) / 2))
    stopifnot(is.function(FUN))

    # We combine to 2 original matrices to create the pool for the permutations
    new_matrix <- rbind(matrix1, matrix2)

    # We create 2 matrices of index with sample_count number of columns and
    # sample_size number of rows. The same column of the 2 matrices cannot
    # contain the same index (i.e.: a sampling without replacement).
    i1 <- matrix(nrow = sample_size, ncol = sample_count)
    i2 <- matrix(nrow = sample_size, ncol = sample_count)
    for (i in seq(1:sample_count)) {
        #j <- sample(1:nrow(new_matrix), sample_size * 2, prob = prob)
        j <- sample(1:nrow(new_matrix), sample_size * 2)
        i1[,i] <- split(j, 1:2)[[1]]
        i2[,i] <- split(j, 1:2)[[2]]
    }

    # We generate a matrix of profiles corresponding to the mean of the columns
    # of the matrices generated during each permutation round.
    get_means <- function(i) {
        replicates <- apply(i, 2, function(j) new_matrix[j,])
        colmeans <- function(x) {
            colMeans(matrix(x, ncol = ncol(new_matrix)))
        }
        apply(replicates, 2, colmeans)
    }
    m1 <- get_means(i1)
    m2 <- get_means(i2)

    # We calculate the scores for each combination of profiles.
    stopifnot(identical(dim(m1), dim(m2)))
    vapply(1:ncol(m1), function(x) FUN(m1[,x], m2[,x], ...), numeric(1))
}
