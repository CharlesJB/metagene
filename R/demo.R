get_demo_bam_files <- function() {
  c(system.file("extdata/align1_rep1.bam", package="metagene"),
    system.file("extdata/align1_rep2.bam", package="metagene"),
    system.file("extdata/align2_rep1.bam", package="metagene"),
    system.file("extdata/align2_rep2.bam", package="metagene"))
}