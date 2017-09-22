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
# A GRanges object splitted into N bins
#
# examples
# gr <- GRanges("chr1", IRanges(c(100, 300), c(200, 500))
# gr <- intoNbins(gr)
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

exon_by_gene_with_observed_transcripts <- function (adb, quantification_files){
	stopifnot(is(adb, "TxDb") | is(adb, "EnsDb"))
	stopifnot(all(tolower(tools::file_ext(quantification_files)) == 'tsv'))
	
	message('Please wait, the process will end in about one minute')
	#retrieve of all transcript_id found in quantification_files (no duplicate)
	transcript_id <- unique(unlist(map(quantification_files, 
					~ fread(.x)[TPM > 0]$transcript_id)))
	transcript_id <- str_replace(transcript_id, '\\.[0-9]+', '')
	
	if(substr(transcript_id[1],1,3) == 'ENS') {
		edb <- EnsDb.Hsapiens.v86
		ebt <- exonsBy(edb, by='tx')
		all_exons_id <- unlist(ebt[which(names(ebt) %in% transcript_id)])$exon_id
		
		#all_exons_id contains ENSE id
		#but unlist(ebg)$exon_id contains only numbers
		
		select(edb, key='ENST00000000233', keytype='TXID', columns = c('GENEID', 'TXID', 'EXONID', 'EXONSEQSTART', 'EXONSEQEND', 'EXONIDX', 'SEQNAME'))
		
		
		gr3 <-
		GRanges(seqnames = c("chr1", "chr2"),
          ranges = IRanges(start=c(1, 3), end=c(4, 9)),
          strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
		
		
		ebg <- exonsBy(edb, by='gene')
		selected <- which(unlist(ebg)$exon_id %in% all_exons_id)
		col_gene_names <- rep(names(ebg), lengths(ebg))
		col_gene_names2 <- col_gene_names[selected]
		ebg2 <- unlist(ebg)[selected]
		split(ebg2, col_gene_names2)
	}
	
	if(substr(transcript_id[1],1,3) == '###') {
		txdb <- org.Hs.eg.db
		ebt <- exonsBy(txdb, by='tx')
		all_exons_id <- unlist(ebt[which(names(ebt) %in% transcript_id)])$exon_id
		
		ebg <- exonsBy(txdb, by='gene')
		selected <- which(unlist(ebg)$exon_id %in% all_exons_id)
		col_gene_names <- rep(names(ebg), lengths(ebg))
		col_gene_names2 <- col_gene_names[selected]
		ebg2 <- unlist(ebg)[selected]
		split(ebg2, col_gene_names2)
	}
}







