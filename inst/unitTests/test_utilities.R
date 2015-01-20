## Test functions present in the  Utilities.R file

### {{{ --- Test setup ---

if(FALSE) {
	library( "RUnit" )
	library( "metagene" )
}

### }}}


###################################################
## Test the getGenes() function
###################################################

## An invalid specie should not be accepted as an argument
test.getGenes_not_valid_specie<- function() {
	obs <- tryCatch(metagene:::getGenes(specie="tomato"), error=conditionMessage)
	exp <- "Incorrect parameter for specie name.\nCurrently supported species are: \"mouse\", \"human\"."
	checkIdentical(obs, exp, msg="getGenes() - An invalid specie argument did not generate an exception with expected message.")
}

## Human specie should return a data set with sepcific header
test.getGenes_human<- function() {
	data<-metagene:::getGenes("human")
	checkTrue(class(data) == "GRanges")
	checkTrue(length(data$feature) > 0, "getGenes() - The returned dataset has not observation.")
}

## Mouse specie should return a data set with sepcific header
test.getGenes_mouse<- function() {
	data<-metagene:::getGenes("mouse")
	checkTrue(class(data) == "GRanges")
	checkTrue(length(data$feature) > 0, "getGenes() - The returned dataset has not observation.")
}
