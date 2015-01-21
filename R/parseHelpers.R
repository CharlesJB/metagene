# Created by Charles Joly Beauparlant
# 2013-11-26

# Scale the values of a vector to fit with predetermined size
#
# Input:
#    values:    The values to scale.
#    domain:    The range to fit the value to.
#
# Output:
#    A vector with the scaled data
scaleVector <- function (values, domain)
{
    to_return <- numeric(domain)
    last_end <- 0
    if (length(values) < domain) {
        ratio <- domain/length(values)
        for (i in seq(1, length(values))) {
            current_end <- round(i * ratio)
            to_return[(last_end+1):current_end] <- values[i]
            last_end <- current_end
        }
    }
    else if (length(values) > domain) {
        ratio <- length(values)/domain
        for (i in seq(1, domain)) {
            current_end <- round(i * ratio)
            to_return[i] <- mean(values[(last_end+1):current_end])
            last_end <- current_end
        }
    }
    else {
        return(values)
    }
    return(to_return)
}
