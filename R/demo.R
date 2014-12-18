get_demo_bam_files <- function() {
  c(system.file("extdata/align1_rep1.bam", package="metagene"),
    system.file("extdata/align1_rep2.bam", package="metagene"),
    system.file("extdata/align2_rep1.bam", package="metagene"),
    system.file("extdata/align2_rep2.bam", package="metagene"))
}

get_not_indexed_bam_file <- function() {
  system.file("extdata/not_indexed.bam", package = "metagene")
}