## Test functions present in the  calcultate_intensity.R file

### {{{ --- Test setup ---

if(FALSE) {
	library( "RUnit" )
	library( "MetaFeatures" )
}

### }}}

###################################################
## Test the prepareBamFiles() function
###################################################

## Zero core should not be accepted as an argument
test.prepareBamFiles_zero_core_number<- function() {
	obs <- tryCatch(MetaFeatures:::prepareBamFiles(c("fileA", "fileB"), 0), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="A zero core number argument did not generate an exception with expected message.")
}

## Negative core number should not be accepted as an argument
test.prepareBamFiles_negative_core_number<- function() {
	obs <- tryCatch(MetaFeatures:::prepareBamFiles(c("fileA", "fileB"), -1), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="A negative core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.prepareBamFiles_not_integer_core_number<- function() {
	obs <- tryCatch(MetaFeatures:::prepareBamFiles(c("fileA", "fileB"),  2.22), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="A decimal core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.prepareBamFiles_not_integer_core_number<- function() {
	obs <- tryCatch(MetaFeatures:::prepareBamFiles(c("fileA", "fileB"),  "NotAInteger"), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="A generic text used as a core number did not generate an exception with expected message.")
}

## All bam files must be in string format
test.prepareBamFiles_file_name_not_in_string_format<- function() {
	obs <- tryCatch(MetaFeatures:::prepareBamFiles(c(1, 2), 2), error=conditionMessage)
	exp <- "At least one BAM file name is not a valid name (not a character string)."
	checkEquals(obs, exp, msg="Integers used as files argument did not generate an exception with expected message.")	
}

## All bam files must exist
test.prepareBamFiles_not_existing_files<- function() {
	obs <- tryCatch(MetaFeatures:::prepareBamFiles(c("NotExistingFile", "NotExistingFile2"), 2), error=conditionMessage)
	exp <- "At least one BAM file does not exist."
	checkEquals(obs, exp, msg="Not existing BAM file used as argument did not generate an exception with expected message.")	
}

###################################################
## Test the getGenes() function
###################################################

## An invalid specie should not be accepted as an argument
test.getGenes_not_valid_specie<- function() {
	obs <- tryCatch(MetaFeatures:::getGenes("tomato"), error=conditionMessage)
	exp <- "Incorrect parameter for specie name.\nCurrently supported species are \"human\" and \"mouse\"."
	checkIdentical(obs, exp, msg="An invalid specie argument did not generate an exception with expected message.")
}
