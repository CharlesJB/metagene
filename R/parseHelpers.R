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

# Parse bed files and convert them in a GRangesList.
# Input:
#    regions:    A vector of bed file names corresponding to the regions to
#                include in the analysis. The file name (minus the extension)
#                will be used as the name of the region.
#    cores:      Number of cores for parallel processing (require parallel
#                package).
#
# Output:
#    A GRangesList. One GRanges by group of features.
#    The names of each GRanges of the list correspond to the name of the group.
prepareRegions <- function(regions, cores=1) {
    # 1. Parse the bed files
    toReturn <- applyOnGroups(regions, cores=cores, FUN=import)

    # 2. Add the names
    getBaseName <- function(bedFileName) {
        return(sub("^([^.]*).*", "\\1", basename(bedFileName)))
    }
    names(toReturn) <- unlist(applyOnGroups(regions, cores=cores, FUN=getBaseName))

    # 3. Convert to GRangesList
    return(GRangesList(toReturn))
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
    parseFeatures <- function(features, bamFiles, cores) {
        raw.counts <- lapply(bamFiles, parseBamFile, features=features, cores=cores)
        names(raw.counts) <- bamFiles
        return(raw.counts)
    }
    raw.data <- lapply(featuresGroups, parseFeatures, bamFiles = bamFiles, cores = cores)
    names(raw.data) <- names(featuresGroups)
    return(raw.data)
}

# Parse a single bam file
#
# Input:
#    bamFile:     The name of the bam file to parse. Must be sorted and indexed.
#    features:    A GRanges corresponding to the regions to parse.
#    cores:       Number of cores for parallel processing (require parallel
#                 package).
#
# Output:
#    A list with an element for every feature to parse. Each element contains
#    a vector of reads expression.
parseBamFile <- function(bamFile, features, cores=1) {
    print(paste("Current bam:", bamFile))
    param <- ScanBamParam(which = features)
    currentReads <- readGAlignments(bamFile, param = param)
    coverages <- coverage(currentReads)
    toReturn <- lapply(coverages[features], as.numeric)
    # Reverse the coverage values of features on minus strand
    if (length(features) > 0) {
        idx <- as.logical(strand(features) == "-")
        toReturn[idx] <- lapply(toReturn[idx], rev)
    }
    return(toReturn)
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
