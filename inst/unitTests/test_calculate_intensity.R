## Test functions present in the  calcultate_intensity.R file

### {{{ --- Test setup ---

if(FALSE) {
	library( "RUnit" )
	library( "MetaMineR" )
}

### }}}

###################################################
## Test the prepareBamFiles() function
###################################################

## Zero core should not be accepted as an argument
test.prepareBamFiles_zero_core_number<- function() {
	obs <- tryCatch(MetaMineR:::prepareBamFiles(c("fileA", "fileB"), 0), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="A zero core number argument did not generate an exception with expected message.")
}

## Negative core number should not be accepted as an argument
test.prepareBamFiles_negative_core_number<- function() {
	obs <- tryCatch(MetaMineR:::prepareBamFiles(c("fileA", "fileB"), -1), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="A negative core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.prepareBamFiles_not_integer_core_number<- function() {
	obs <- tryCatch(MetaMineR:::prepareBamFiles(c("fileA", "fileB"),  2.22), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="A decimal core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.prepareBamFiles_not_integer_core_number<- function() {
	obs <- tryCatch(MetaMineR:::prepareBamFiles(c("fileA", "fileB"),  "NotAInteger"), error=conditionMessage)
	exp <- "The number of cores has to be a positive integer."
	checkIdentical(obs, exp, msg="A generic text used as a core number did not generate an exception with expected message.")
}

## All bam files must be in string format
test.prepareBamFiles_file_name_not_in_string_format<- function() {
	obs <- tryCatch(MetaMineR:::prepareBamFiles(c(1, 2), 2), error=conditionMessage)
	exp <- "At least one BAM file name is not a valid name (a character string)."
	checkEquals(obs, exp, msg="Integers used as files argument did not generate an exception with expected message.")	
}

## All bam files must exist
test.prepareBamFiles_not_existing_files<- function() {
	obs <- tryCatch(MetaMineR:::prepareBamFiles(c("NotExistingFile", "NotExistingFile2"), 2), error=conditionMessage)
	exp <- "At least one BAM file does not exist."
	checkEquals(obs, exp, msg="Not existing BAM file used as argument did not generate an exception with expected message.")	
}

###################################################
## Test the getGenes() function
###################################################

## An invalid specie should not be accepted as an argument
test.getGenes_not_valid_specie<- function() {
	obs <- tryCatch(MetaMineR:::getGenes(specie="tomato"), error=conditionMessage)
	exp <- "Incorrect parameter for specie name.\nCurrently supported species are \"human\" and \"mouse\"."
	checkIdentical(obs, exp, msg="An invalid specie argument did not generate an exception with expected message.")
}

## Human specie should return a data set with sepcific header
test.getGenes_human<- function() {
	data<-MetaMineR:::getGenes("human")
	checkEquals(names(data), c("feature", "strand", "space", "start_position", "end_position"),  "The returned dataset does not have the expected column names")
	checkTrue(length(data[, 1]) > 0, "The returned dataset has not observation.")
}

## Mouse specie should return a data set with sepcific header
test.getGenes_mouse<- function() {
	data<-MetaMineR:::getGenes("mouse")
	checkEquals(names(data), c("feature", "strand", "space", "start_position", "end_position"),  "The returned dataset does not have the expected column names")
	checkTrue(length(data[, 1]) > 0, "The returned dataset has not observation.")
}

###################################################
## Test the prepareFeatures() function
###################################################

## An invalid specie should not be accepted as an argument
test.prepareFeatures_not_valid_specie<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::prepareFeatures(temp_file, "tomato"), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "Incorrect parameter for specie name.\nCurrently supported species are \"human\" and \"mouse\"."
	checkIdentical(obs, exp, msg="An invalid specie argument did not generate an exception with expected message.")
}

## Zero core should not be accepted as an argument
test.prepareFeatures_zero_core_number<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::prepareFeatures(temp_file, cores=0), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The number of cores has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="A zero core number argument did not generate an exception with expected message.")
}

## Negative core number should not be accepted as an argument
test.prepareFeatures_negative_core_number<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::prepareFeatures(temp_file, cores=-2), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The number of cores has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="A negative core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.prepareFeatures_not_numeric_core_number<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::prepareFeatures(temp_file, cores=2.22), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The number of cores has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="A decimal core number argument did not generate an exception with expected message.")
}

## Zero maximum distance should not be accepted as an argument
test.prepareFeatures_zero_max_distance<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::prepareFeatures(temp_file, maxDistance=0), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The maximum distance has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="A zero maximum distance argument did not generate an exception with expected message.")
}

## Negative maximum distance should not be accepted as an argument
test.prepareFeatures_negative_max_distance<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::prepareFeatures(temp_file, maxDistance=-2), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The maximum distance has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="A negative maximum distance argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a maximum distance argument 
test.prepareFeatures_not_integer_max_distance<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::prepareFeatures(temp_file, maxDistance=2.33), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The maximum distance has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="A decimal maximum distance argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a core number argument 
test.prepareFeatures_not_integer_max_distance<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::prepareFeatures(temp_file, maxDistance="NotAInteger"), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The maximum distance has to be a positive numeric with no decimals."
	checkIdentical(obs, exp, msg="A generic text used as a maximum distance did not generate an exception with expected message.")
}

## All features files must be in string format
test.prepareFeatures_file_name_not_in_string_format<- function() {
	obs <- tryCatch(MetaMineR:::prepareFeatures(c(1, 2)), error=conditionMessage)
	exp <- "At least one features file name is not a valid name (a character string)."
	checkEquals(obs, exp, msg="Integers used as features files argument did not generate an exception with expected message.")	
}

## All features files must exist
test.prepareFeatures_not_existing_files<- function() {
	obs <- tryCatch(MetaMineR:::prepareFeatures(c("NotExistingFile", "NotExistingFile2")), error=conditionMessage)
	exp <- "At least one features file does not exist."
	checkEquals(obs, exp, msg="Not existing features files used as argument did not generate an exception with expected message.")	
}

###################################################
## Test the extractReadsInRegion() function
###################################################

## Zero starting position should not be accepted as an argument
test.extractReadsInRegion_zero_start<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::extractReadsInRegion(temp_file, "X", 0, 100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The starting position has to be a positive integer."
	checkIdentical(obs, exp, msg="A zero maximum distance argument did not generate an exception with expected message.")
}

## Negative starting position should not be accepted as an argument
test.extractReadsInRegion_negative_start<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::extractReadsInRegion(temp_file, "X", -4, 100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The starting position has to be a positive integer."
	checkIdentical(obs, exp, msg="A negative maximum distance argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a starting position argument 
test.extractReadsInRegion_not_integer_start<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::extractReadsInRegion(temp_file, "X", 2.33, 100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The starting position has to be a positive integer."
	checkIdentical(obs, exp, msg="A decimal maximum distance argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a starting position argument 
test.extractReadsInRegion_not_integer_start<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::extractReadsInRegion(temp_file, "X", "NotAInteger", 100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The starting position has to be a positive integer."
	checkIdentical(obs, exp, msg="A generic text used as a starting position did not generate an exception with expected message.")
}

## Zero ending position should not be accepted as an argument
test.extractReadsInRegion_zero_end<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::extractReadsInRegion(temp_file, "X", 10, 0), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The ending position has to be a positive integer."
	checkIdentical(obs, exp, msg="A zero ending position argument did not generate an exception with expected message.")
}

## Negative ending position should not be accepted as an argument
test.extractReadsInRegion_negative_end<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::extractReadsInRegion(temp_file, "X", 3, -100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
exp <- "The ending position has to be a positive integer."
checkIdentical(obs, exp, msg="A negative ending position argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an ending position argument 
test.extractReadsInRegion_not_integer_end<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::extractReadsInRegion(temp_file, "X", 3, 33.2), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The ending position has to be a positive integer."
	checkIdentical(obs, exp, msg="A decimal ending position argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a starting position argument 
test.extractReadsInRegion_not_integer_end<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(MetaMineR:::extractReadsInRegion(temp_file, "X", 20, "NotAInteger"), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The ending position has to be a positive integer."
	checkIdentical(obs, exp, msg="A generic text used as an ending position did not generate an exception with expected message.")
}

