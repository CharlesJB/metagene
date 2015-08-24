#' Get BAM filenames for demo
#' 
#' @return A vector of BAM filenames
#' 
#' @examples
#' bam_files <- get_demo_bam_files()
get_demo_bam_files <- function() {
    c(system.file("extdata/align1_rep1.bam", package="metagene"),
        system.file("extdata/align1_rep2.bam", package="metagene"),
        system.file("extdata/align2_rep1.bam", package="metagene"),
        system.file("extdata/align2_rep2.bam", package="metagene"),
        system.file("extdata/ctrl.bam", package="metagene"))
}

#' Get regions filenames for demo
#' 
#' @return A vector of regions filenames
#' 
#' @examples
#' regions <- get_demo_regions()
get_demo_regions <- function() {
    c(system.file("extdata/list1.bed", package="metagene"),
        system.file("extdata/list2.bed", package="metagene"))
}

get_not_indexed_bam_file <- function() {
    system.file("extdata/not_indexed.bam", package = "metagene")
}

get_different_seqnames_bam_file <- function() {
    system.file("extdata/different_header.bam", package = "metagene")
}

get_coverage_bam_file <- function() {
    system.file("extdata/coverage.bam", package = "metagene")
}

get_coverage_region <- function() {
    system.file("extdata/list_coverage.bed", package = "metagene")
}
