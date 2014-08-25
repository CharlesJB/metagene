# Created by Astrid Deschenes
# 2014-07-28

## Test functions present in the parseFeatures.R file

### {{{ --- Test setup ---

if(FALSE) {
	library( "RUnit" )
	library( "metagene" )
}

### }}}

###################################################
## Test the parseFeatures() function
###################################################

## Zero core should not be accepted as an argument
test.parseFeatures_zero_core_number<- function() {
	obs <- tryCatch(metagene:::parseFeatures(c("fileA", "fileB"), cores = 0), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="parseFeatures() - A zero core number argument did not generate an exception with expected message.")
}

## Negative core number should not be accepted as an argument
test.parseFeatures_negative_core_number<- function() {
	obs <- tryCatch(metagene:::parseFeatures(c("fileA", "fileB"), cores = -1), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="parseFeatures() - A negative core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.parseFeatures_not_integer_core_number<- function() {
	obs <- tryCatch(metagene:::parseFeatures(c("fileA", "fileB"),  cores = 2.22), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="parseFeatures() - A decimal core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.parseFeatures_string_core_number<- function() {
	obs <- tryCatch(metagene:::parseFeatures(c("fileA", "fileB"),  cores = "NotAInteger"), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="parseFeatures() - A generic text used as a core number did not generate an exception with expected message.")
}

## All bam files must be in string format
test.parseFeatures_file_name_not_in_string_format<- function() {
	obs <- tryCatch(metagene:::parseFeatures(c(1, 2), cores = 1), error=conditionMessage)
	exp <- "At least one BAM file name is not a valid name (a character string)."
	checkEquals(obs, exp, msg="parseFeatures() - Integers used as files argument did not generate an exception with expected message.")	
}

## All bam files must exist
test.parseFeatures_not_existing_files<- function() {
	obs <- tryCatch(metagene:::parseFeatures(c("NotExistingFile", "NotExistingFile2"), 2), error=conditionMessage)
	exp <- "At least one BAM file does not exist."
	checkEquals(obs, exp, msg="parseFeatures() - Not existing BAM file used as argument did not generate an exception with expected message.")	
}

## Zero maximum distance should not be accepted as an argument
test.prepareFeatures_zero_max_distance<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::parseFeatures(c(temp_file), maxDistance=0), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The maximum distance has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="parseFeatures() - A zero maximum distance argument did not generate an exception with expected message.")
}

## Negative maximum distance should not be accepted as an argument
test.prepareFeatures_negative_max_distance<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::parseFeatures(c(temp_file), maxDistance=-2), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The maximum distance has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="parseFeatures() - A negative maximum distance argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a maximum distance argument 
test.prepareFeatures_not_integer_max_distance<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::parseFeatures(c(temp_file), maxDistance=2.33), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The maximum distance has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="parseFeatures() - A decimal maximum distance argument did not generate an exception with expected message.")
}
