# Created by Astrid Deschenes
# 2014-07-28

## Test functions present in the parseHelpers.R file

### {{{ --- Test setup ---

if(FALSE) {
	library( "RUnit" )
	library( "metagene" )
}

### }}}

###################################################
## Test the extractReadsInRegion() function
###################################################

## Zero starting position should not be accepted as an argument
test.extractReadsInRegion_zero_start<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::extractReadsInRegion(temp_file, "X", 0, 100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The starting position has to be a positive integer."
	checkIdentical(obs, exp, msg="extractReadsInRegion() - A zero maximum distance argument did not generate an exception with expected message.")
}

## Negative starting position should not be accepted as an argument
test.extractReadsInRegion_negative_start<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::extractReadsInRegion(temp_file, "X", -4, 100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The starting position has to be a positive integer."
	checkIdentical(obs, exp, msg="extractReadsInRegion() - A negative maximum distance argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a starting position argument 
test.extractReadsInRegion_not_integer_start<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::extractReadsInRegion(temp_file, "X", 2.33, 100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The starting position has to be a positive integer."
	checkIdentical(obs, exp, msg="extractReadsInRegion() - A decimal maximum distance argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a starting position argument 
test.extractReadsInRegion_not_integer_start<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::extractReadsInRegion(temp_file, "X", "NotAInteger", 100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The starting position has to be a positive integer."
	checkIdentical(obs, exp, msg="extractReadsInRegion() - A generic text used as a starting position did not generate an exception with expected message.")
}

## Zero ending position should not be accepted as an argument
test.extractReadsInRegion_zero_end<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::extractReadsInRegion(temp_file, "X", 10, 0), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The ending position has to be a positive integer."
	checkIdentical(obs, exp, msg="extractReadsInRegion() - A zero ending position argument did not generate an exception with expected message.")
}

## Negative ending position should not be accepted as an argument
test.extractReadsInRegion_negative_end<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::extractReadsInRegion(temp_file, "X", 3, -100), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The ending position has to be a positive integer."
	checkIdentical(obs, exp, msg="extractReadsInRegion() - A negative ending position argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as an ending position argument 
test.extractReadsInRegion_not_integer_end<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::extractReadsInRegion(temp_file, "X", 3, 33.2), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The ending position has to be a positive integer."
	checkIdentical(obs, exp, msg="A decimal ending position argument did not generate an exception with expected message.")
}

## Something other than an integer number should not be accepted as a starting position argument 
test.extractReadsInRegion_not_integer_end<- function() {
	temp_file<-tempfile()
	file.create(temp_file)
	obs <- tryCatch(metagene:::extractReadsInRegion(temp_file, "X", 20, "NotAInteger"), error=conditionMessage, finally={if (file.exists(temp_file)){file.remove(temp_file)}})
	exp <- "The ending position has to be a positive integer."
	checkIdentical(obs, exp, msg="extractReadsInRegion() - A generic text used as an ending position did not generate an exception with expected message.")
}
