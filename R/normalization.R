read.BAM <- function(bam_file) {
    stopifnot(file.exists(bam_file))
    bf <- Rsamtools::BamFile(bam_file)
    lengths <- GenomeInfoDb::seqlengths(bf)

    read.chr <- function(chr) {
        gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(1, lengths[[chr]]))
        flag_p <- Rsamtools::scanBamFlag(isNotPassingQualityControls = FALSE,
                                        isMinusStrand = FALSE)
        flag_m <- Rsamtools::scanBamFlag(isNotPassingQualityControls = FALSE,
                                         isMinusStrand = TRUE)
        param_p <- Rsamtools::ScanBamParam(flag = flag_p, which = gr)
        param_m <- Rsamtools::ScanBamParam(flag = flag_m, which = gr)
        plus <- GenomicAlignments::readGAlignments(bam_file, param = param_p)
        minus <- GenomicAlignments::readGAlignments(bam_file, param = param_m)
        list("-" = BiocGenerics::end(minus),
             "+" = BiocGenerics::start(plus) - 1L)
    }

    chr <- GenomeInfoDb::seqlevels(Rsamtools::BamFile(bam_file))
    res <- lapply(chr, read.chr)
    names(res) <- chr

    # Remove empty chr
    is_zero <- function(x) {
        all(vapply(x, function(y) length(y) == 0, logical(1)))
    }
    i <- vapply(res, is_zero, logical(1))
    res[i] <- NULL

    res
}
