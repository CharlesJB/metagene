NCIS.internal <- function (chip.pos, input.pos, shift.size = 100, min.binsize = 100,
    max.binsize = 20000, binsize.shift = 100, min.stop.binsize = 100,
    chr.vec = NULL, chr.len.vec = NULL) {
    if (is.null(chr.vec)) {
        chip.name <- names(chip.pos)
        input.name <- names(input.pos)
        chr.vec <- intersect(chip.name, input.name)
        if (length(chr.vec) < length(chip.name)) {
            cat("Control sample doesn't have chromosome", setdiff(chip.name,
                chr.vec), "which are in ChIP sample. These chromosomes are ignored.\n")
        }
        if (length(chr.vec) < length(input.name)) {
            cat("ChIP sample doesn't have chromosome", setdiff(chip.name,
                chr.vec), "which are in control sample. These chromosomes are ignored.\n")
        }
    }
    nchip <- sum(sapply(chr.vec, function(x) sapply(chip.pos[[x]],
        length)))
    ninput <- sum(sapply(chr.vec, function(x) sapply(input.pos[[x]],
        length)))
    r.seq.depth <- nchip/ninput
    sizevec <- rep(c(1, 2, 5), times = 3) * 10^rep(2:4, each = 3)
    sizevec <- sizevec[sizevec <= max.binsize]
    sizevec <- sizevec[sizevec >= min.binsize]
    norm.est <- rep(1000, length(sizevec))
    if (!is.null(chr.len.vec) & length(setdiff(chr.vec, names(chr.len.vec))) ==
        0) {
        chr.end.max <- chr.len.vec[chr.vec]
    }
    else {
        chr.end.max <- sapply(chr.vec, function(chr) max(c(max(chip.pos[[chr]][["+"]]),
            max(chip.pos[[chr]][["-"]]), max(input.pos[[chr]][["+"]]),
            max(input.pos[[chr]][["-"]]))))
    }
    binsize.est <- -1
    for (si in 1:length(sizevec)) {
        binsize <- sizevec[si]
        bindata <- bin.data(chip.pos, input.pos, binsize, shift.size = shift.size,
            chr.vec = chr.vec, chr.end.max = chr.end.max)
        res <- est.norm.med.search(bindata$chip, bindata$input)
        if (binsize < binsize.shift) {
            norm.est[si] <- res
        }
        else {
            bindata <- bin.data(chip.pos, input.pos, binsize,
                shift.size = shift.size, shift.half.size = TRUE,
                chr.vec = chr.vec, chr.end.max = chr.end.max)
            res2 <- est.norm.med.search(bindata$chip, bindata$input)
            norm.est[si] <- (res + res2)/2
        }
        if (si > 1 & binsize.est < 0) {
            if (norm.est[si] >= norm.est[si - 1]) {
                est <- norm.est[si - 1]
                binsize.est <- sizevec[si - 1]
            }
            if (si == length(sizevec) & binsize.est < 0) {
                est <- norm.est[si]
                binsize.est <- sizevec[si]
            }
        }
        if (binsize.est > 0 & binsize >= min.stop.binsize) {
            break
        }
    }
    return(list(est = est, binsize.est = binsize.est, r.seq.depth = r.seq.depth,
        pi0 = est/r.seq.depth))
}

bin.data <- function (chip.pos, input.pos, binsize, shift.size = 100, shift.half.size = FALSE,
    zero.filter = TRUE, by.strand = FALSE, chr.vec = NULL, chr.end.max = NULL,
    by.chr = FALSE)
{
    if (is.null(chr.vec)) {
        chip.name <- names(chip.pos)
        input.name <- names(input.pos)
        chr.vec <- intersect(chip.name, input.name)
        if (length(chr.vec) < length(chip.name)) {
            cat("Control sample doesn't have chromosome", setdiff(chip.name,
                chr.vec), "which are in ChIP sample. These chromosomes are ignored.\n")
        }
        if (length(chr.vec) < length(input.name)) {
            cat("ChIP sample doesn't have chromosome", setdiff(input.name,
                chr.vec), "which are in control sample. These chromosomes are ignored.\n")
        }
    }
    if (is.null(chr.end.max)) {
        chr.end.max <- sapply(chr.vec, function(chr) max(c(max(chip.pos[[chr]][["+"]]),
            max(chip.pos[[chr]][["-"]]), max(input.pos[[chr]][["+"]]),
            max(input.pos[[chr]][["-"]]))))
    }
       chr.len <- ceiling((chr.end.max + shift.size)/binsize)                                                            [26/120]
    chip.f <- list()
    chip.r <- list()
    input.f <- list()
    input.r <- list()
    if (!by.strand) {
        chip <- list()
        input <- list()
    }
    for (chr in chr.vec) {
        if (shift.half.size) {
            bk <- c(0, seq(from = round(binsize/2), by = binsize,
                length.out = chr.len[[chr]] + 1))
        }
        else {
            bk <- seq(from = 0, by = binsize, length.out = chr.len[[chr]] +
                1)
        }
        bk.f <- bk - shift.size
        bk.r <- c(0, bk + shift.size)
        chip.f[[chr]] <- hist(chip.pos[[chr]][["+"]], breaks = bk.f,
            plot = FALSE)$counts
        chip.r[[chr]] <- hist(chip.pos[[chr]][["-"]], breaks = bk.r,
            plot = FALSE)$counts[-1]
        input.f[[chr]] <- hist(input.pos[[chr]][["+"]], breaks = bk.f,
            plot = FALSE)$counts
        input.r[[chr]] <- hist(input.pos[[chr]][["-"]], breaks = bk.r,
            plot = FALSE)$counts[-1]
        if (zero.filter) {
            ind <- chip.f[[chr]] + chip.r[[chr]] + input.f[[chr]] +
                input.r[[chr]] > 0
            chip.f[[chr]] <- chip.f[[chr]][ind]
            chip.r[[chr]] <- chip.r[[chr]][ind]
            input.f[[chr]] <- input.f[[chr]][ind]
            input.r[[chr]] <- input.r[[chr]][ind]
        }
        if (!by.strand) {
            chip[[chr]] <- chip.f[[chr]] + chip.r[[chr]]
            input[[chr]] <- input.f[[chr]] + input.r[[chr]]
        }
            }
    if (by.strand) {
        if (by.chr) {
            return(list(chip.f = chip.f, chip.r = chip.r, input.f = input.f,
                input.r = input.r))
        }
        else {
            return(list(chip.f = unlist(chip.f, use.names = FALSE),
                chip.r = unlist(chip.r, use.names = FALSE), input.f = unlist(input.f,
                  use.names = FALSE), input.r = unlist(input.r,
                  use.names = FALSE)))
        }
    }
    else {
        if (by.chr) {
            return(list(chip = chip, input = input))
        }
        else {
            return(list(chip = unlist(chip, use.names = FALSE),
                input = unlist(input, use.names = FALSE)))
        }
    }
}
