# Created by Astrid Deschenes
# 2014-08-05

## Test functions present in the parseRegions.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
}

### }}}

###################################################
## Test the parseRegions() function
###################################################

## An invalid specie should not be accepted as an argument
test.parseRegions_not_valid_specie<- function() {
    obs <- tryCatch(metagene:::parseRegions(specie="tomato"), error=conditionMessage)
    exp <- "Incorrect parameter for specie name.\nCurrently supported species are: \"mouse\", \"human\"."
    checkIdentical(obs, exp, msg="parseRegions() - An invalid specie argument did not generate an exception with expected message.")
}

## Zero core should not be accepted as an argument
test.parseRegions_zero_core_number<- function() {
    obs <- tryCatch(metagene:::parseRegions(bamFiles=c("fileA", "fileB"), cores = 0), error=conditionMessage)
    exp <- "The number of cores has to be a positive integer."
    checkIdentical(obs, exp, msg="parseRegions() - A zero core number argument did not generate an exception with expected message.")
}

## Negative core number should not be accepted as an argument
test.parseRegions_negative_core_number<- function() {
    obs <- tryCatch(metagene:::parseRegions(bamFiles=c("fileA", "fileB"), cores = -1), error=conditionMessage)
    exp <- "The number of cores has to be a positive integer."
    checkIdentical(obs, exp, msg="parseRegions() - A negative core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.parseRegions_not_integer_core_number<- function() {
    obs <- tryCatch(metagene:::parseRegions(bamFiles=c("fileA", "fileB"),  cores = 2.22), error=conditionMessage)
    exp <- "The number of cores has to be a positive integer."
    checkIdentical(obs, exp, msg="parseRegions() - A decimal core number argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an core number argument 
test.parseRegions_string_core_number<- function() {
    obs <- tryCatch(metagene:::parseRegions(bamFiles=c("fileA", "fileB"),  cores = "NotAInteger"), error=conditionMessage)
    exp <- "The number of cores has to be a positive integer."
    checkIdentical(obs, exp, msg="parseRegions() - A generic text used as a core number did not generate an exception with expected message.")
}

## All bam files must be in string format
test.parseRegions_file_name_not_in_string_format<- function() {
    obs <- tryCatch(metagene:::parseRegions(bamFiles=c(1, 2), cores = 1), error=conditionMessage)
    exp <- "At least one BAM file name is not a valid name (a character string)."
    checkEquals(obs, exp, msg="parseRegions() - Integers used as files argument did not generate an exception with expected message.")	
}

## All bam files must exist
test.parseRegions_not_existing_files<- function() {
    obs <- tryCatch(metagene:::parseRegions(bamFiles=c("NotExistingFile", "NotExistingFile2"), 2), error=conditionMessage)
    exp <- "At least one BAM file does not exist."
    checkEquals(obs, exp, msg="parseRegions() - Not existing BAM file used as argument did not generate an exception with expected message.")	
}

## Negative padding size should not be accepted as an argument
test.parseRegions_negative_padding_size<- function() {
    temp_file<-tempfile()
    file.create(temp_file)
    obs <- tryCatch(metagene:::parseRegions(bamFiles = c(temp_file), paddingSize=-2), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
    exp <- "The padding size has to be a positive numeric with no decimals."
    checkIdentical(obs, exp, msg="parseRegions() - A negative padding size argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a padding size argument 
test.parseRegions_not_integer_padding_size<- function() {
    temp_file<-tempfile()
    file.create(temp_file)
    obs <- tryCatch(metagene:::parseRegions(bamFiles = c(temp_file), paddingSize=2.33), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
    exp <- "The padding size has to be a positive numeric with no decimals."
    checkIdentical(obs, exp, msg="parseRegions() - A decimal padding size argument did not generate an exception with expected message.")
}
