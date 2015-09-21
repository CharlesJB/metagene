## Test functions present in the demo.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
}

### }}}


###################################################
## Test the get_demo_bam_files() function
###################################################

test.get_demo_bam_files <- function() {
    bam_files <- get_demo_bam_files()
    checkTrue(length(bam_files) == 5)
    checkTrue(all(file.exists(bam_files)))
}

###################################################
## Test the get_demo_regions() function
###################################################

test.get_demo_regions <- function() {
    regions <- get_demo_regions()
    checkTrue(length(regions) == 2)
    checkTrue(all(file.exists(regions)))
}

###################################################
## Test the get_demo_metagene() function
###################################################

test.get_demo_metagene <- function() {
    mg <- get_demo_metagene()
    checkTrue(all(class(mg) == c("metagene", "R6")))
    bam_files <- mg$get_params()$bam_files
    regions <- names(mg$get_regions())
    checkTrue(length(bam_files) == 5)
    checkTrue(all(file.exists(bam_files)))
}

###################################################
## Test the get_demo_design() function
###################################################

## get_demo_design
test.get_demo_design <- function() {
    obs <- get_demo_design()
    samples <- c("align1_rep1.bam", "align1_rep2.bam", "align2_rep1.bam",
                 "align2_rep2.bam", "ctrl.bam")
    exp <- data.frame(Samples = samples, align1 = c(1L, 1L, 0L, 0L, 2L),
                      align2 = c(0L, 0L, 1L, 1L, 2L))
    exp$Samples <- as.character(exp$Samples)
    checkIdentical(obs, exp)
}
