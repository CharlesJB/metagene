# Created by Charles Joly Beauparlant
# 2013-11-26

# Split a GRanges into N bins
#
# Originally posted by Martin Morgan:
# https://stat.ethz.ch/pipermail/bioconductor/2012-September/047923.html
#
# param gr A GRanges with only one seqnames value.
# param n Number of bins to produce.
#
# return
#   A GRanges object splitted into N bins
#
# examples
#   gr <- GRanges("chr1", IRanges(c(100, 300), c(200, 500))
#   gr <- intoNbins(gr)
intoNbins <- function(gr, n = 10) {
    stopifnot(class(gr) == "GRanges")
    stopifnot(length(gr) > 0)
    stopifnot(is.numeric(n))
    stopifnot(n > 0)
    if (any(GenomicRanges::width(gr) < n)) {
        stop("all 'width(gr)' must be >= 'n'")
    }
    d <- width(gr) / n
    dd <- cumsum(rep(d, each=n))
    mask <- logical(n); mask[1] <- TRUE
    dd <- dd - rep(dd[mask], each=n)

    starts <- round(rep(GenomicRanges::start(gr), each=n) + dd)
    ends <- c(starts[-1], 0) - 1L
    ends[rev(mask)] <- end(gr)

    gr <- gr[rep(seq_along(gr), each=n)]
    GenomicRanges::ranges(gr) <- IRanges(starts, ends)
    gr
}

#' Extract Entrez genes promoters from TxDb object.
#'
#' @param txdb A valid \code{TxDb} object.
#' @param upstream The number of nucleotides upstream of TSS.
#' @param downstream The number of nucleotides downstream of TSS.
#'
#' @return A \code{GRanges} object that contains the promoters infos.
#'
#' @examples
#' \dontrun{
#'     # require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'     txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'     promoters_hg19 <- get_promoters_txdb(txdb)
#' }
#'
#' @import GenomeInfoDb
get_promoters_txdb <- function(txdb, upstream = 1000, downstream = 1000) {
    stopifnot(is(txdb, "TxDb"))
    GenomicFeatures::promoters(GenomicFeatures::genes(txdb),
                               upstream = upstream, downstream = downstream)
}
