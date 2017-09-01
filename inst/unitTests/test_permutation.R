## Test functions present in the permutation.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
}

### }}}
####################################################
### Test the metagene$permutation_test() function
####################################################

test.metagene_permutation_test_valid <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	ratio_normalized_intersect <- 
		perm_fun(tab1[, .(moy=mean(value)), by=bin]$moy, tab2[, .(moy=mean(value)), by=bin]$moy)
	#print(paste("ratio_normalized_intersect :", ratio_normalized_intersect))
	permutation_results <- permutation_test(tab1, tab2, sample_size = 50,
										sample_count = 1000, FUN = perm_fun)

	#print(paste("p-value :", sum(ratio_normalized_intersect >= permutation_results) / length(permutation_results)))
	print(TRUE)
}

test.metagene_permutation_test_unvalid_table1_table2_are_the_same <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	# ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab1, sample_size = 50,
										sample_count = 1000, FUN = perm_fun),
					error = conditionMessage)
	exp <- "!identical(table1, table2) is not TRUE"
	checkIdentical(obs, exp)
}

test.metagene_permutation_test_table1_or_2_unvalid_class <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- c(1:50)
	tab2 <- c(1:50)
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	# ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = 50,
										sample_count = 1000, FUN = perm_fun),
					error = conditionMessage)
	exp <- "is.data.table(table1) is not TRUE"
	checkIdentical(obs, exp)
	
	tab1 <- tab[which(tab$design == "align1"),]
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = 50,
										sample_count = 1000, FUN = perm_fun),
					error = conditionMessage)
	exp <- "is.data.table(table2) is not TRUE"
	checkIdentical(obs, exp)
}

test.metagene_permutation_test_table1_wilder_than_table2 <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),-1]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	#ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = 50,
										sample_count = 1000, FUN = perm_fun),
					error = conditionMessage)
	exp <- "ncol(table1) == ncol(table2) is not TRUE"
	checkIdentical(obs, exp)
}

test.metagene_permutation_test_sample_size_unvalid_class <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	#ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = "test",
										sample_count = 1000, FUN = perm_fun),
					error = conditionMessage)
	exp <- "is.numeric(sample_size) is not TRUE"
	checkIdentical(obs, exp)
}

test.metagene_permutation_test_sample_size_unvalid_value <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	#ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = -10,
										sample_count = 1000, FUN = perm_fun),
					error = conditionMessage)
	exp <- "sample_size > 0 is not TRUE"
	checkIdentical(obs, exp)
}

test.metagene_permutation_test_sample_size_too_big <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	#ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = dim(tab1)[1],
										sample_count = 1000, FUN = perm_fun),
					error = conditionMessage)
	exp <- "sample_size <= ((nrow(table1)/bincount + nrow(table2)/bincount)/2) is not TRUE"
	checkIdentical(obs, exp)
}

test.metagene_permutation_test_sample_count_unvalid_class <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	#ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = 50,
										sample_count = "test", FUN = perm_fun),
					error = conditionMessage)
	exp <- "is.numeric(sample_count) is not TRUE"
	checkIdentical(obs, exp)
}

test.metagene_permutation_test_sample_count_unvalid_value <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	#ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = 50,
										sample_count = -10, FUN = perm_fun),
					error = conditionMessage)
	exp <- "sample_count > 0 is not TRUE"
	checkIdentical(obs, exp)
}

test.metagene_permutation_test_FUN_invalid_class <- function() {
	mg <- get_demo_metagene()
	design <- get_demo_design()
	mg$produce_table(design = design)
	tab <- mg$get_table()
	tab <- tab[which(tab$region == "list1"),]
	tab1 <- tab[which(tab$design == "align1"),]
	tab2 <- tab[which(tab$design == "align2"),]
	
	library(similaRpeak)
	perm_fun <- function(profile1, profile2) {
		similarity(profile1, profile2)[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
	}
	#ratio_normalized_intersect <- perm_fun(tab1$value, tab2$value)
	obs <- tryCatch(permutation_test(tab1, tab2, sample_size = 50,
										sample_count = 1000, FUN = "test"),
					error = conditionMessage)
	exp <- "is.function(FUN) is not TRUE"
	checkIdentical(obs, exp)
}




###################################################
## Test the permutation_test() function
###################################################

#mg <- get_demo_metagene()
#mg$produce_table()
#m1 <- mg$get_table()$list1$align1_rep1$input
#m2 <- mg$get_table()$list1$align2_rep1$input
#
### Valid use
#test.permutation_test_valid_use <- function() {
#    sample_size <- min(nrow(m1), nrow(m2))
#    sample_count <- 100
#    FUN = function(a, b) { mean(a) - mean(b) }
#    res <- permutation_test(m1, m2, sample_size = sample_size,
#                            sample_count = sample_count, FUN = FUN)
#    checkTrue(is.numeric(res))
#    checkTrue(length(res) ==  sample_count)
#}
#
### Valid sample_size max value
#test.permutation_test_valid_sample_size_max_value <- function() {
#    sample_size <- (nrow(m1) + nrow(m2)) / 2
#    sample_count <- 100
#    FUN = function(a, b) { mean(a) - mean(b) }
#    res <- permutation_test(m1, m2, sample_size = sample_size,
#                            sample_count = sample_count, FUN = FUN)
#    checkTrue(is.numeric(res))
#    checkTrue(length(res) ==  sample_count)
#}
#
### Invalid m1 class
#test.permutation_test_invalid_m1_class <- function() {
#    obs <- tryCatch(permutation_test(1), error=conditionMessage)
#    exp <- "is.matrix(matrix1) is not TRUE"
#    checkIdentical(obs, exp)
#}
#
### Invalid m2 class
#test.permutation_test_invalid_m2_class <- function() {
#    obs <- tryCatch(permutation_test(m1, 1), error=conditionMessage)
#    exp <- "is.matrix(matrix2) is not TRUE"
#    checkIdentical(obs, exp)
#}
#
### Invalid dim m1 not equal dim m2
#test.permutation_test_invalid_dim_m1_not_equal_dim_m2 <- function() {
#    obs <- tryCatch(permutation_test(m1, matrix()), error=conditionMessage)
#    exp <- "ncol(matrix1) == ncol(matrix2) is not TRUE"
#    checkIdentical(obs, exp)
#}
#
### Invalid sample_size class
#test.permutation_test_invalid_sample_size_class <- function() {
#    obs <- tryCatch(permutation_test(m1, m2, sample_size = "a"),
#                    error=conditionMessage)
#    exp <- "is.numeric(sample_size) is not TRUE"
#    checkIdentical(obs, exp)
#}
#
### Invalid sample_size value too large
#test.permutation_test_invalid_sample_size_value_too_large <- function() {
#    sample_size <- (nrow(m1) + nrow(m2)) / 2 + 1
#    obs <- tryCatch(permutation_test(m1, m2, sample_size = sample_size,
#                                     sample_count = 10),
#                    error=conditionMessage)
#    exp <- "sample_size <= ((nrow(matrix1) + nrow(matrix2))/2) is not TRUE"
#    checkIdentical(obs, exp)
#}
#
### Invalid sample_count class
#test.permutation_test_invalid_sample_count_class <- function() {
#    obs <- tryCatch(permutation_test(m1, m2, sample_size = 100,
#                                     sample_count = "a"),
#                    error=conditionMessage)
#    exp <- "is.numeric(sample_count) is not TRUE"
#    checkIdentical(obs, exp)
#}
#
### Invalid fun class
#test.permutation_test_invalid_fun_class <- function() {
#    obs <- tryCatch(permutation_test(m1, m2, sample_size = 10,
#                                     sample_count = 100, FUN = "a"),
#                    error=conditionMessage)
#    exp <- "is.function(FUN) is not TRUE"
#    checkIdentical(obs, exp)
#}
