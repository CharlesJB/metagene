#' Perform a permutation test on 2 tables
#'
#' The goal of this function is to calculate the values of a test performed
#' by \code{FUN} after each of \code{sample_count} permutations.
#'
#' Each round of the permutation test, two new matrices will be randomly
#' sampled from using the combination of the two original tables. The means
#' of each columns will be calculated to produce the vectors that will be sent
#' \code{FUN}.
#' 
#' @param table1 The first table.
#' @param table2 The second table.
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
#' \dontrun{
#' # Get some tables
#' mg <- get_demo_metagene()
#' mg$produce_table()
#' tab <- mg$get_table()
#' tab <- tab[which(tab$region == "list1"),]
#' tab1 <- tab[which(tab$design == "align1_rep1"),]
#' tab2 <- tab[which(tab$design == "align2_rep2"),]
#'
#' # Perform permutation test
#' sample_size <- min(nrow(tab1), nrow(tab2))
#' FUN = function(a, b) { mean(a) - mean(b) } # Dummy function for demo purpose
#' # A sample_count >= 1000 should be used in a real analysis
#' permutation_results <- permutation_test(m1, m2, sample_size = sample_size,
#'                                         sample_count = 10, FUN = FUN)
#' }
#'
#' @export
permutation_test <- function(table1, table2, sample_size, sample_count, FUN,
                             ...) {
    stopifnot(is.data.table(table1))
    stopifnot(is.data.table(table2))
	stopifnot(!identical(table1, table2))
    stopifnot(ncol(table1) == ncol(table2))
    stopifnot(is.numeric(sample_size))
    stopifnot(is.numeric(sample_count))
	stopifnot(sample_size > 0)
	stopifnot(sample_count > 0)
	bincount <- length(unique(table1$bin))
    stopifnot(sample_size <= ((nrow(table1)/bincount + nrow(table2)/bincount)/2))
    stopifnot(is.function(FUN))

    # We combine to 2 original matrices to create the pool for the permutations
    new_table <- rbind(table1, table2)
	print(new_table)

	#to new_table to new_matrix in order to get value in row and bin in columns
	#to facilitate following sample and mean computations
	new_matrix <- matrix(new_table$value, ncol=bincount, byrow=TRUE)
	print(new_matrix)
	
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
	print(i2)
	print(i2)

	
	# We generate a matrix of profiles corresponding to the mean by bin
	# from the matrices i1 and i2 generated during each permutation round.
    get_means <- function(i) {
		# "2" indicates column
		# apply gets row of new_table selected during the sampling
        replicates <- apply(i, 2, function(j) new_matrix[j,])
		print(replicates)
        #compute the mean of each replicates
		colmeans <- function(x) {
            colMeans(matrix(x, ncol = sample_count))
        }
        apply(replicates, 2, colmeans)
    }
    m1 <- get_means(i1)
	print(m1)
    m2 <- get_means(i2)
	print(m2)

    # We calculate the scores for each combination of profiles.
    stopifnot(identical(dim(m1), dim(m2)))
    x <- vapply(1:ncol(m1), function(x) FUN(m1[,x], m2[,x], ...), numeric(1))
	print(x)
	x
}
