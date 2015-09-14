## Test functions present in the demo.R file

### {{{ --- Test setup ---

if(FALSE) {
  library( "RUnit" )
  library( "metagene" )
}

### }}}

###################################################
## Test the get_demo_design() function
###################################################

base_msg <- "get_demo_design -"

## get_demo_design
test.get_demo_design <- function() {
  obs <- get_demo_design()
  samples <- c("align1_rep1.bam", "align1_rep2.bam", "align2_rep1.bam",
	       "align2_rep2.bam", "ctrl.bam")
  exp <- data.frame(Samples = samples, align1 = c(1L, 1L, 0L, 0L, 2L),
		    align2 = c(0L, 0L, 1L, 1L, 2L))
  exp$Samples <- as.character(exp$Samples)
  msg <- paste(base_msg, "A valid function call did not give the expected results." )
  checkIdentical(obs, exp, msg)
}
