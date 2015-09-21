## Test functions present in the  Utilities.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
}

### }}}

###################################################
## Test the metagene:::intoNbins() function
###################################################

## Valid gr default n
test.intonbins_valid_gr_default_n <- function() {
    gr <- rtracklayer::import(get_demo_regions()[1])
    gr_bins <- metagene:::intoNbins(gr)
    checkTrue(class(gr_bins) == "GRanges")
    checkTrue(length(gr_bins) == length(gr) * 10)
    checkTrue(length(BiocGenerics::setdiff(gr, gr_bins)) == 0)
}

## Valid gr valid n
test.intonbins_valid_gr_valid_n <- function() {
    gr <- rtracklayer::import(get_demo_regions()[1])
    n <- 100
    gr_bins <- metagene:::intoNbins(gr, n)
    checkTrue(class(gr_bins) == "GRanges")
    checkTrue(length(gr_bins) == length(gr) * n)
    checkTrue(length(BiocGenerics::setdiff(gr, gr_bins)) == 0)
}

## Valid gr valid n length equal width
test.intonbins_valid_gr_valid_n_length_equal_gr_width <- function() {
    gr <- rtracklayer::import(get_demo_regions()[1])
    n <- unique(width(gr))
    gr_bins <- metagene:::intoNbins(gr, n)
    checkTrue(class(gr_bins) == "GRanges")
    checkTrue(length(gr_bins) == length(gr) * n)
    checkTrue(length(BiocGenerics::setdiff(gr, gr_bins)) == 0)
}

## Invalid gr class
test.intonbins_invalid_gr_class <- function() {
    obs <- tryCatch(metagene:::intoNbins(1), error = conditionMessage)
    exp <- "class(gr) == \"GRanges\" is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid gr empty
test.intonbins_invalid_gr_empty <- function() {
    obs <- tryCatch(metagene:::intoNbins(GRanges()), error = conditionMessage)
    exp <- "length(gr) > 0 is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid n class
test.intonbins_invalid_n_class <- function() {
    gr <- rtracklayer::import(get_demo_regions()[1])
    n <- "a"
    obs <- tryCatch(metagene:::intoNbins(gr, n), error = conditionMessage)
    exp <- "is.numeric(n) is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid n zero
test.intonbins_invalid_n_zero <- function() {
    gr <- rtracklayer::import(get_demo_regions()[1])
    n <- 0
    obs <- tryCatch(metagene:::intoNbins(gr, n), error = conditionMessage)
    exp <- "n > 0 is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid n negative
test.intonbins_invalid_n_negative <- function() {
    gr <- rtracklayer::import(get_demo_regions()[1])
    n <- -1
    obs <- tryCatch(metagene:::intoNbins(gr, n), error = conditionMessage)
    exp <- "n > 0 is not TRUE"
    checkIdentical(obs, exp)
}

## Invalid n greater than width gr
test.intonbins_invalid_n_greater_than_width_gr <- function() {
    gr <- rtracklayer::import(get_demo_regions()[1])
    n <- unique(width(gr)) + 1
    obs <- tryCatch(metagene:::intoNbins(gr, n), error = conditionMessage)
    exp <- "all 'width(gr)' must be >= 'n'"
    checkIdentical(obs, exp)
}
