## Test functions present in the  calcultate_intensity.R file

###################################################
## Test the prepareBamFiles() function
###################################################

## Zero core should not be accepted as an argument
test.prepareBamFiles_zero_core_number<- function() {
	checkException(prepareBamFiles(c("fileA", "fileB"), 0))
}

## Negative core number should not be accepted as an argument
test.prepareBamFiles_negative_core_number<- function() {
	checkException(prepareBamFiles(c("fileA", "fileB"), -1))
}

## Something other than an integer number should not be accepted as an core number argument 
test.prepareBamFiles_not_integer_core_number<- function() {
	checkException(prepareBamFiles(c("fileA", "fileB"),  2.22))
	checkException(prepareBamFiles(c("fileA", "fileB"), "NotAInteger"))
}