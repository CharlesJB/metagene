# Created by Charles Joly Beauparlant
# 2013-11-26

# Create a graph
#
# Input:
#    matricesGroups:    A list with every groupsto include in the graph. A
#                       group is one or more combination of featureGroup
#                       designGroup. The names must correspond to the names
#                       of the matrix object returned from the parse function.
#    data:              The object returned by the parse functions.
#    binSize:           The number of nucleotides in each bin for the
#                       bootstrap step.
#    alpha:             Confidence interval.
#    sampleSize:       Number of time each bin will be resampled (should be
#                       at least 1000).
#    cores:             Number of cores for parallel processing (require
#                       parallel package).
# Ouput:
#     The data.frame used to produce the graph
plotMatrices <- function(matricesGroups, data, binSize=100, alpha=0.05, sampleSize=1000, cores=1) {
    # 1. Extract and/or combine relevant matrices
    for (matrixName in names(matricesGroups)) {
        conditions <- unlist(names(data$matrix) %in% matricesGroups[[matrixName]])
        matricesGroups[[matrixName]] <- do.call(rbind, data$matrix[conditions])
    }
    # 2. Bootstrap
    bootstrap <- function(x, bootstrapCores=1) {
        return(bootstrapAnalysis(x, binSize=binSize, alpha=alpha, sampleSize=sampleSize, cores=bootstrapCores))
    }
    bootstrapResults <- applyOnGroups(matricesGroups, cores=1, FUN=bootstrap, bootstrapCores=cores)

    # 3. Prepare data.frame
    DF <- getDataFrame(bootstrapResults, range=data$range)

    # 4. Create graph
    plotGraphic(DF, paste(names(matricesGroups), collapse=" vs "))
    return(DF)
}

# Perform the bootstrap analysis
#
# Input:
#    currentMatrix:    The matrix to use for the bootstrap analysis
#    binSize:          The number of nucleotides in each bin for the
#                      bootstrap step.
#    alpha:            Confidence interval.
#    sampleSize:       Number of time each bin will be resampled (should be
#                      at least 1000).
#    cores:            Number of cores for parallel processing (require
#                      parallel package).
#
# Output:
#    A list with the mean of the bootstraped vector and the quartile of order
#    alpha/2 and (1-alpha/2)
bootstrapAnalysis <- function(currentMatrix, binSize, alpha, sampleSize, cores=1) {
    binnedMatrix <- binMatrix(currentMatrix, binSize)
    if (cores > 1) {
        library(parallel)
        bootResults <- mclapply(1:ncol(binnedMatrix), function(x) binBootstrap(binnedMatrix[,x], alpha=alpha, sampleSize=sampleSize), mc.cores=cores)
    } else {
        bootResults <- lapply(1:ncol(binnedMatrix), function(x) binBootstrap(binnedMatrix[,x], alpha=alpha, sampleSize=sampleSize))
    }
    bootResults <- do.call(rbind, bootResults)
    toReturn <- list()
    toReturn$mean <- unlist(bootResults[,1])
    toReturn$qinf <- as.numeric(unlist(bootResults[,2]))
    toReturn$qsup <- as.numeric(unlist(bootResults[,3]))
    return(toReturn)
}

# Bin matrix columns
#
# Input:
#    data:       The matrix to bin.
#    binSize:    The number of nucleotides in each bin.
#
# OUPUT:
#    A matrix with each column representing the mean of binSize nucleotides.
binMatrix <- function(data,binSize)
{
    stopifnot(binSize %% 1 == 0)
    stopifnot(binSize > 0)
    # Not sure if best solution. If the binSize is not a multiple of number of
    # row in data, I silently trim data to make it even
    if (ncol(data) %% binSize != 0) {
        remainder <- ncol(data) %% binSize
        if (remainder == 1) {
        # If remainder is one, we remove it at the end
            data <- data[,-ncol(data)]
        } else if (remainder %% 2 != 0) {
        # If remainder is odd, we remove 1 more at the end
            toRemove <- floor(remainder / 2)
            data <- data[,(toRemove+1):(ncol(data)-toRemove-1)]
	} else {
        # If remainder is even, we remove the same quantity on both end
            toRemove <- remainder / 2
            data <- data[,(toRemove+1):(ncol(data)-toRemove)]
        }
    }

    splitSum <- function(x, bs) {
        if (bs < length(x)) {
            return(tapply(x, (seq_along(x)-1) %/% bs, sum))
        } else {
            return(x)
        }
    }
    return(unname(t(apply(data, 1, splitSum, binSize))))
}

# Estimate mean and confidence interval of a column using bootstrap.
#
# Input:
#    data:          A vector representing a column from the binned matrix.
#    alpha:         Confidence interval.
#    sampleSize:    Number of time each bin will be resampled (should be at
#                   least 1000).
#    cores:         Number of cores for parallel processing (require parallel
#                   package).
#
# OUPUT:
#    A list with the mean of the bootstraped vector and the quartile of order
#    alpha/2 and (1-alpha/2).
binBootstrap <- function(data, alpha, sampleSize, cores=1)
{
    # TODO: it would probably be more efficient to parallelize here
    size <- length(data)
    # TODO: try to sample data directly
    S <- matrix(replicate(sampleSize, data[sample(1:length(data),size,replace=TRUE)]), nrow=sampleSize)
    # TODO: We should also be able to plot non-bootstrapped mean
    mean <- mean(sapply(1:sampleSize, function(i){mean(S[i,])}))
    qinf <- quantile(sapply(1:sampleSize, function(i){mean(S[i,])}), prob=alpha/2)
    qsup <- quantile(sapply(1:sampleSize, function(i){mean(S[i,])}), prob=(1-alpha/2))
    liste <- list(mean=mean, qinf=qinf, qsup=qsup)
    return(liste)
}

# Convert the bootstrapped data into a data.frame
#
# Input:
#    bootstapData:    Data produced during the bootstrap analysis
#
# Output:
#    A data.frame with the condensed results from the main data structure
#    Columns:
#         * Groups: name of current group
#         * distances: the number of bin for each entry
#         * means: the means to plot
#         * qinf: the lower end of the confidence interval
#         * qsup: the higher end of the confidence interval
getDataFrame <- function(bootstrapData, range) {
    grid = seq(range[1], range[2],length=length(bootstrapData[[1]]$mean))
    DF = data.frame (
        Groups <- factor(rep(names(bootstrapData), each=length(grid))),
        distances <- rep(grid, length(bootstrapData)),
        means <- c(sapply(1:length(bootstrapData), function(x) bootstrapData[[x]]$mean)),
        qinf <-  c(sapply(1:length(bootstrapData), function(x) bootstrapData[[x]]$qinf)),
        qsup <-  c(sapply(1:length(bootstrapData), function(x) bootstrapData[[x]]$qsup))
    )
    colnames(DF) <- c("Groups", "distances", "means", "qinf", "qsup")
    return(DF)
}

# Produce a plot with based on a data.frame
#
# Input:
#    DF:       The data frame produced by the plot.getDataFrame function
#    title:    The title of the graph
#
# Ouput:
#    The graph that is printed on the current device.
plotGraphic <- function(DF, title) {
    # TODO: add x label
    p <- ggplot(DF, aes(x=distances, y=means, ymin=qinf, ymax=qsup)) +
    geom_ribbon(aes(fill=Groups), alpha=0.3) +
    geom_line(aes(color=Groups),size=1,litype=1,bg="transparent")+
    theme(panel.grid.major = element_line())+
    theme(panel.grid.minor = element_line())+
    theme(panel.background = element_blank())+
    theme(panel.background = element_rect())+
    theme_bw()+
    theme(axis.title.x = element_blank())+
    ylab("Mean of RPM for each 100 bps")+
    ggtitle(title)
    print(p)
}
