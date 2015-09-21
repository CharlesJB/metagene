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

#' Get a demo metagene object
#'
#' @return A metagene object
#'
#' @examples
#' mg <- get_demo_metagene()
get_demo_metagene <- function() {
    regions <- get_demo_regions()
    bam_files <- get_demo_bam_files()
    metagene$new(regions = regions, bam_files = bam_files)
}

#' Get a demo design object
#'
#' @return A \code{data.frame} corresponding to a valid design.
#'
#' @examples
#' mg <- get_demo_design()
get_demo_design <- function() {
    file_design <- system.file("extdata/design.txt", package = "metagene")
    read.table(file_design, header = TRUE, stringsAsFactors = FALSE)
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

get_narrowpeak_region <- function() {
    system.file("extdata/list1.narrowPeak", package = "metagene")
}

get_broadpeak_region <- function() {
    system.file("extdata/list1.broadPeak", package = "metagene")
}
