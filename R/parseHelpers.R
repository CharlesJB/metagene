# Created by Charles Joly Beauparlant
# 2013-11-26

# Apply a function on every groups of the main data structure
#
# Input:
#    groups:    The main data structure
#    cores:     Number of cores for parallel processing (require parallel
#               package).
#    FUN:       The function to apply on every groups
#    ...:       Extract arguments for FUN
#
# Output:
#    The result of the function applied
applyOnGroups <- function(groups, cores=1, FUN, ...) {
    if (cores > 1) {
        library(parallel)
        return(mclapply(groups, function(x) FUN(x, ...), mc.cores=cores))
    } else {
        return(lapply(groups, function(x) FUN(x, ...)))
    }
}

# Parse bed files and convert them in a list of data.frames
# Input:
#    regions:    A vector of bed file names corresponding to the regions to
#                include in the analysis. The file name (minus the extension)
#                will be used as the name of the region.
#    cores:      Number of cores for parallel processing (require parallel
#                package).
#
# Output:
#    A list of data.frame. One data.frame by group of features.
#    The names of each element of the list correspond to the name of the group.
prepareRegions <- function(regions, cores=1) {
    # 1. Parse the bed files
    readBedFile <- function(bedFileName) {
        currentBed <- read.table(bedFileName, header=FALSE, stringsAsFactors=FALSE)
        # We only need the infos of the first three columns
        colnames(currentBed)[1:3] <- c("space", "start_position", "end_position")
        # TODO: should we check is start_position is always smaller than end_position?
        return(currentBed)
    }
    toReturn <- applyOnGroups(regions, cores=cores, FUN=readBedFile)

    # 2. Add the names
    getBaseName <- function(bedFileName) {
        return(sub("^([^.]*).*", "\\1", basename(bedFileName)))
    }
    names(toReturn) <- unlist(applyOnGroups(regions, cores=cores, FUN=getBaseName))

    return(toReturn)
}

# Get the genomic position of the paddings
#
# Input:
#    regionsGroups:    The output of prepareRegions or prepareFeatures
#                      functions.
#    side:             The side of the padding to prepare. Either 'left' or
#                      right'.
#    paddingSize:      The length padding we want to add on each side of each
#                      regions.
#    cores:            Number of cores for parallel processing (require
#                      parallel package). The number of cores has to be a 
#                      positive integer.
#
# Output:
#    A list of data.frames, each data.frame corresponding to the paddings for
#    a group of features.
prepareRegionsPaddings <- function(regionsGroups, side, paddingSize=2000, cores=1) {
    # Check prerequisites

    # Side must be either 'left' or 'right'
    if (! side %in% c("left", "right")) {
        stop("Side argument must be either 'left' or 'right'.")
    }

    # The number of cores has to be a positive integer
    if(!is.numeric(cores) || cores <= 0) {
        stop("The number of cores has to be a positive integer.")
    }

    getPaddings <- function(currentRegions) {
        currentPaddings <- currentRegions
        if (side == "left") {
            currentPaddings$start_position <- currentRegions$start_position - paddingSize
            currentPaddings$end_position <- currentRegions$start_position
            currentPaddings$start_position[currentPaddings$start_position < 0] <- 0
        } else {
            currentPaddings$start_position <- currentRegions$end_position
            currentPaddings$end_position <- currentRegions$end_position + paddingSize
        }
        return(currentPaddings)
    }

    toReturn <- applyOnGroups(regionsGroups, cores=cores, FUN=getPaddings)
    names(toReturn) <- names(regionsGroups)
    return(toReturn)
}

# Parse multiple bamFiles
#
# Input:
#    bamFiles:          A vector of bam files
#    featuresGroups:    A list of data.frame. One data.frame by group of
#                       features.
#                       The names of each element of the list correspond to
#                       the name of the group.
#    cores:             Number of cores for parallel processing (require
#                       parallel package).
#
# Output:
#    A list of list of list that contains the raw counts for every features
#    groups:
#        * First level:     One entry by features group
#        * Second level:    One entry per bamFile in current features group
#        * Third level:     One entry per feature in current features group
parseBamFiles <- function(bamFiles, featuresGroups, cores=1) {
    parseFeatures <- function(features, currentName, bamFiles, cores) {
        print(currentName)
        raw.counts <- lapply(bamFiles, parseBamFile, features=features, cores=cores)
        names(raw.counts) <- bamFiles
        return(raw.counts)
    }
    raw.data <- lapply(1:length(featuresGroups), function(x) parseFeatures(featuresGroups[[x]], names(featuresGroups)[x], bamFiles=bamFiles, cores=cores))
    names(raw.data) <- names(featuresGroups)
    return(raw.data)
}

# Parse a single bam file
#
# Input:
#    bamFile:     The name of the bam file to parse. Must be sorted and indexed.
#    features:    A data.frame with the infos for every features to parse
#                 Must have the folowing columns: feature, strand, space,
#                 start_position and end_position
#    cores:       Number of cores for parallel processing (require parallel
#                 package).
#
# Output:
#    A list with an element for every feature to parse. Each element contains
#    a vector of reads expression.
parseBamFile <- function(bamFile, features, cores=1) {
    print(paste("Current bam:", bamFile))
    extractReadsDensity <- function(feature, bamFile) {
        # Extract raw counts
        currentReads <- extractReadsInRegion(bamFile, feature$space, feature$start_position, feature$end_position)
        # TODO: sometimes there are NA, need to find why and replace this current hack
                currentReads <- currentReads[!is.na(currentReads$pos),]
        vectorResult <- convertReadsToDensity(currentReads, feature)

        # If on negative strand, invert the current vector
        if (feature$start_position > feature$end_position) {
            vectorResult <- rev(vectorResult)
        }
        return(vectorResult)
    }
    return(applyOnGroups(1:nrow(features), cores=cores, FUN=function(x) extractReadsDensity(features[x,], bamFile=bamFile)))
}

# Extract reads from BAM file that overlap with a specified genomic region.
#
# Input:
#    bamFile:    Path to the bam file.
#    chr:        Current chromosome.
#    start:      Starting position of the current region.
#    end:        Ending position of the current region.
#
# Prerequisites:
# The BAM file must exist.
# The chromosome name must be in character format.
# The starting position has to be a positive integer.
# The ending position has to be a positive integer.
#
# Output:
#    A data.frame containing every reads overlapping the current genomic
#    region:
#        * rname
#        * pos
#        * qwidth
extractReadsInRegion <- function(bamFile, chr, start, end) {

    # Check prerequisites

    # The BAM file name must be of string type
    if (!is.character(bamFile)) {
        stop("The BAM file name is not a valid name (a character string).")
    }

    # The BAM file name must exist
    if (! file.exists(bamFile)) {
        stop("The BAM file does not exist.")
    }

    #The chromosome name must be of character type
    if (!is.character(chr)) {
        stop("The chromosome name is not a valid name (a character string).")
    }

    # The starting position has to be a positive integer
    if(!is.numeric(start) || start <= 0) {
        stop("The starting position has to be a positive integer.")
    }

    # The ending position has to be a positive integer
    if(!is.numeric(end) || end <= 0) {
        stop("The ending position has to be a positive integer.")
    }

    if (start > end) {
        tmp <- start
        start <- end
        end <- tmp
    }
    df <- data.frame()
    which <- GRanges(seqnames=Rle(chr), ranges=IRanges(start, end))
    what <- c("rname", "pos", "qwidth")
    param <- ScanBamParam(which=which, what=what)
    bam <- scanBam(as.character(bamFile), param=param)
    bam <- unname(bam)
    return(do.call("DataFrame", bam))
}

# Convert a list of read in a vector of positions.
#
# Input:
#    currentReads:        The list of read to parse.
#    currentFeature:      The feature to parse.
#
# Output:
#    A vector with the coverage of every positions calculated from the reads
#    around the max distance from TSS.
convertReadsToDensity <- function(currentReads, currentFeature) {
    maxSize <- abs(currentFeature$end_position - currentFeature$start_position)
    start <- min(currentFeature$start_position, currentFeature$end_position)
    vectorResult <- numeric(maxSize)
    if (nrow(currentReads) > 0) {
        positions <- unlist(mapply(function(x,y) seq(x, x+y), currentReads$pos - currentFeature$start_position, currentReads$qwidth-1))
        # TODO: add unit test -> first position < 0 used to return FALSE and skip current feature (when it was &&)
        positions <- positions[positions > 0 & positions <= maxSize]
        vectorResult <- tabulate(positions, nbins=maxSize)
    }
    return(vectorResult)
}

# Convert the raw counts into reads per million aligned (rpm)
#
# Input:
#    rawCounts:              The data structure returned by the parseBamFiles
#                            function.
#    bamFilesDescription:    The data.frame obtained with the prepareBamFiles
#                            function.
#    cores:                  Number of cores for parallel processing (require
#                            parallel package).
#
# Output:
#    A list of list of list that contains the rpm for every features groups:
#    First level:     One entry by features group
#    Second level:    One entry per bamFile in current features group
#    Third level:     One entry per feature in current features group
rawCountsToRPM <- function(rawCounts, bamFilesDescription, cores=1) {
    convertToRPM <- function(bamFile, bamFilesDescription, cores) {
        alignedCount <- bamFilesDescription[bamFilesDescription$bam == names(bamFile),]$alignedCount
        bamFile <- bamFile[[1]]
        return(applyOnGroups(bamFile, cores, FUN=function(x) { x / (alignedCount / 1000000) }))
    }
    rpmCounts <- list()
    for (i in 1:length(rawCounts)) {
        currentFeaturesGroupName <- names(rawCounts)[i]
        featuresGroup <- rawCounts[[i]]
        rpmCounts[[currentFeaturesGroupName]] <-
            lapply(1:length(featuresGroup), function(x) convertToRPM(featuresGroup[x], bamFilesDescription=bamFilesDescription, cores=cores))
        names(rpmCounts[[currentFeaturesGroupName]]) <- names(featuresGroup)
    }
    return(rpmCounts)
}

# Distribute the bam filenames in their respective groups.
#
# Input:
#    featuresGroupsNames:    The name of the group of features in the current
#                            analysis.
#    bamFiles:               A vector of bam filename(s).
#    design:                 A matrix explaining the relationship between
#                            multiple samples.
#                                * One line per samples.
#                                * One column per group of samples. For
#                                  example, biological replicates and
#                                  corresponding controls are in the same group.
#                                      1: treatment file(s)
#                                      2: control file(s)
#
# Output:
#    A list of list of list.
#    * First level is the name of the group.
#    * Second level is the metadata of the group:
#        featureName:    The name of the group of features.
#        designName:     The name of the column in the design file.
#        bamFiles:       The list of bam files associated with this
#                        combination of featureName/designName
#    * Third level is the actual list of bam files.
#    Note1: If there are groups of feature and groups of bam file, every
#           possible combination of feature groups / bam groups will be created.
#           Otherwise, the is a group for each group of features.
#    Note2: The names of the elements the bam list are the same as the ones
#           in the design file, which are also the original bam file names.
#           The content of the each element of the list of bam file is the
#           sorted bam file name, which should be use when parsing bam.
prepareGroups <- function(featuresGroupsNames, bamFiles, design=NULL) {
    allGroups <- list()
    for (featureName in featuresGroupsNames) {
        # If there is no design file, we put every bamFile in each groups
        if (is.null(design)) {
            allGroups[[featureName]] <- list()
            allGroups[[featureName]]$featureName <- featureName
            allGroups[[featureName]]$bamFiles <- as.list(as.character(bamFiles$bam))
            names(allGroups[[featureName]]$bamFiles) <- as.character(bamFiles$oldBam)
        } else {
            for (i in 2:ncol(design)) {
                name <- paste(featureName, colnames(design)[i], sep="_")
                bam <- design[design[,i] != 0,][,1]
                designName <- colnames(design)[i]
                # We want the sorted name of the bam files as the content of the list
                j <- numeric()
                for (bamFile in bam) {
                    j <- append(j, which(bamFiles$oldBam == bamFile))
                }
                sortedBamNames <- as.character(bamFiles[j,]$bam)
                # Create the actual element of the list
                allGroups[[name]] <- list()
                allGroups[[name]]$featureName <- featureName
                allGroups[[name]]$designName <- designName
                allGroups[[name]]$bamFiles <- as.list(sortedBamNames)
                names(allGroups[[name]]$bamFiles) <- bam
            }
        }
    }
    return(allGroups)
}

# Substract controls from a single group
#
# Input:
#    group:           A group extracted from the main data structure
#    data.rpm:        The normalized data for every bam files.
#    design:          The line from matrix explaining the relationship
#                     between current samples.
#    controlCores:    Number of cores for parallel processing (require
#                     parallel package).
#
# Output:
#    The group extracted from the main data structure from which the controls
#    were substracted then deleted.
removeControls <- function(group, data.rpm, design=NULL, controlCores=1) {
    if (!is.null(design)) {
        designName <- group$designName

        # 1. Extract relevant design columns
        currentDesign <- design[,c(1,which(colnames(design) == designName))]
        treatmentNames <- as.character(currentDesign[currentDesign[,2] == 1, 1])
        controlNames <- as.character(currentDesign[currentDesign[,2] == 2, 1])

        # 2 Get the relevant group of rpm
        current.rpm <- data.rpm[[group$featureName]]

        # 3. Merge controls
        # TODO
        if (length(controlNames) > 1) {
            print("TODO")
            # uuu <- data.frame(zzz <- c(1,2,3), xxx <- c(4,5,6), yyy <- c(7,8,9))
        } else {
            mergedControls <- current.rpm[[controlNames[1]]]
        }

        # 4. Substract merged control from every treatment samples
	if (length(controlNames) == 0) {
            group$noCTRL <- current.rpm[treatmentNames]
        } else {
		substractFeature <- function(treatment, i) {
		    currentFeatureTreatment <- unlist(treatment[i])
		    currentFeatureControl <- unlist(mergedControls[i])
		    newValues <- currentFeatureTreatment - currentFeatureControl
		    newValues[newValues < 0] <- 0
		    return(newValues)
		}
		substractControl <- function(currentTreatment) {
		    newTreatment <- applyOnGroups(1:length(currentTreatment), cores=controlCores, FUN=function(x) substractFeature(currentTreatment, x))
		    return(newTreatment)
		}
		treatment.rpm <- current.rpm[names(current.rpm) %in% treatmentNames]
		group$noCTRL <- lapply(treatment.rpm, substractControl)
		names(group$noCTRL) <- treatmentNames
		return(group)
	}
    } else {
        current.rpm <- data.rpm[[group$featureName]]
        group$noCTRL <- current.rpm
    }
    return(group)
}

# Convert list of vectors in matrix for a single group
#
# Input:
#    group:    Correspond to an element in the data subsection of the main
#              data structure.
#    level:    The names of the element to merge.
#
# Output:
#    The group in input with an extra element named matrix (group$matrix).
mergeMatrix <- function(group, level) {
    mergeMatrixFixedLengthFeatures <- function(currentGroup) {
        return(lapply(currentGroup[[level]], function(x) do.call(rbind, x)))
    }
    group$matrix <- do.call(rbind, mergeMatrixFixedLengthFeatures(group))
    return(group)
}

# Resize the vectors of every group so they have the same length
#
# Input:
#    group:         A list that contains a list of list of vectors
#    level:         The names of the element of the group to merge.
#    domain:        The target length for the vectors
#    scaleCores:    Number of cores for parallel processing (require parallel
#                   package).
#
# Output:
#    The same group that was used in input with an extra element named scaled.
scaleVectors <- function(group, level, domain, scaleCores=1) {
    scaleBamFile <- function(bamFile, domain, scaleCores=1) {
        if (scaleCores > 1) {
            return(mclapply(bamFile, scaleVector, domain, mc.cores=scaleCores))
        } else {
            return(lapply(bamFile, scaleVector, domain))
        }
    }
    group$scaled <- lapply(group[[level]], scaleBamFile, domain=domain, scaleCores=scaleCores)
    names(group$scaled) <- names(group[[level]])
    return(group)
}

# Scale the values of a vector to fit with predetermined size
#
# Input:
#    values:    The values to scale.
#    domain:    The range to fit the value to.
#
# Output:
#    A vector with the scaled data
scaleVector <- function(values, domain) {
        to_return <- numeric(domain)
        if (length(values) < domain) {
                ratio <- domain / length(values)
                for (i in seq(1, length(values))) {
                        current_start <- round((i - 1) * ratio + 1)
                        current_end <- round(i * ratio)
                        to_return[current_start:current_end] <- values[i]
                }
        }
        else if (length(values) > domain) {
                ratio <- length(values) / domain
                for (i in seq(1, domain)) {
                        current_start <- round((i - 1) * ratio + 1)
                        current_end <- round(i * ratio)
                        to_return[i] <- mean(values[current_start:current_end])
                }

        }
        else {
                return(values)
        }
        return(to_return)
}
