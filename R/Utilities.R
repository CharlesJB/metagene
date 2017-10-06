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


#' Extract exons by genes for which data are available in quantification files
#'
#' @param adb A valid \code{EnsDb} object.
#' @param quantification_files the quantification files. A vector of pathes.
#' @param downstream The number of nucleotides downstream of TSS.
#'
#' @return A \code{GRangesList} object containing exons by genes for which 
#'                                data are available in quantification files.
#'
#' @examples
#' \dontrun{
#'        require(EnsDb.Hsapiens.v86)
#'        edb <- EnsDb.Hsapiens.v86
#'        quantification_files <- 'file_path'
#'         ebgwot <- exon_by_gene_with_observed_transcripts(edb, 
#'                                                        quantification_files)
#'         bed_file_content_gr <- GRanges("chr16",ranges = IRanges(start=23581002, 
#'                                                                end=23596356))
#'         bed_file_filter(ebgwot, bed_file_content_gr)
#' }
#'
#' @import EnsDb.Hsapiens.v86
exon_by_gene_with_observed_transcripts <- function (adb, quantification_files){
    stopifnot((is(adb, "TxDb") | is(adb, "EnsDb")))
    stopifnot(all(tolower(tools::file_ext(quantification_files)) == 'tsv'))
    
    message(paste0('Please wait, the process will end in about one minute ',
                                'depending on quantification_files size'))
    #retrieve of all transcript_id found in quantification_files (no duplicate)
    transcript_id <- unique(unlist(map(quantification_files, 
                    ~ fread(.x)[TPM > 0]$transcript_id)))
    transcript_id <- str_replace(transcript_id, '\\.[0-9]+', '')
    
    if(is(adb, "EnsDb") & substr(transcript_id[1],1,3) == 'ENS') {
        edb <- EnsDb.Hsapiens.v86    
        slct <- unique(ensembldb::select(edb, key=all_exons_id, 
                                    keytype='EXONID', 
                                    columns = c('GENEID', 'EXONID', 
                                            'EXONSEQSTART', 'EXONSEQEND',
                                            'SEQSTRAND','EXONIDX', 'TXID', 
                                            'SEQNAME')))
        
        slct <- slct[which(slct$TXID %in% transcript_id),]
        slct$SEQSTRAND <- str_replace(slct$SEQSTRAND, '-1', '-')
        slct$SEQSTRAND <- str_replace(slct$SEQSTRAND, '1', '+')
        slct$SEQSTRAND <- str_replace(slct$SEQSTRAND, 'NA', '*')
        
        gr <-
        GRanges(seqnames = paste0('chr',slct$SEQNAME),
          ranges = IRanges(start=slct$EXONSEQSTART, end=slct$EXONSEQEND),
          strand = slct$SEQNAME, exon_id = slct$EXONID, 
                                                    gene_id = slct$GENEID)
        return(split(gr, gr$gene_id))
    } else if(is(adb, "TxDb") & is.numeric(transcript_id[1])) {
        #TODO
    } else {
        message(paste0('AnnotationDataBase (adb) or quantification_files not',
                    ' supported. Only EnsDB AnnotationDataBase',
                    ' and corresponding transcripts quantification files', 
                    ' with ENST ID are supported.'))
    }
}

#' Extract a list of ranges defined by the bed_file_content_gr argument from the 
#'    exon_by_gene_with_observed_transcripts GRangesList. Equivalent to the 
#'    exonsByOverlaps of GenomicFeatures.
#'
#' @param exon_by_gene_with_observed_transcripts A \code{GRangesList} object 
#'        provided by the exon_by_gene_with_observed_transcripts function.
#' @param bed_file_content_gr A \code{GRanges} object containing ranges of 
#'        interest.
#' @param reduce If the returned \code{GRanges} object will be reduce or not.
#'
#' @return A \code{GRanges} object that contains exons by genes selected.
#'
#' @examples
#' \dontrun{
#'        require(EnsDb.Hsapiens.v86)
#'        edb <- EnsDb.Hsapiens.v86
#'        quantification_files <- 'file_path'
#'         ebgwot <- exon_by_gene_with_observed_transcripts(edb, 
#'                                                        quantification_files)
#'         bed_file_content_gr <- GRanges("chr16",ranges = IRanges(start=23581002, 
#'                                                                end=23596356))
#'         bed_file_filter(ebgwot, bed_file_content_gr)
#' }
#'
bed_file_filter <- function (exon_by_gene_with_observed_transcripts, 
                                bed_file_content_gr, reduce = TRUE) {
    ebgwot <- exon_by_gene_with_observed_transcripts
    if(reduce){
        IRanges::reduce(unlist(ebgwot)[(
            unlist(start(ebgwot)) >= start(bed_file_content_gr) & 
            unlist(end(ebgwot)) <= end(bed_file_content_gr) & 
            as.character(unlist(seqnames(ebgwot))) == as.character(
                                                seqnames(bed_file_content_gr))
        ),])
    } else {
        unlist(ebgwot)[(
            unlist(start(ebgwot)) >= start(bed_file_content_gr) & 
            unlist(end(ebgwot)) <= end(bed_file_content_gr) & 
            as.character(unlist(seqnames(ebgwot))) == as.character(
                                                seqnames(bed_file_content_gr))
        ),]
    }
}


#' Transforms the bed_file_filter function output into a file.BED readable by
#'                                                                metagene.
#'
#' @param bed_file_filter_result A \code{GRanges} object : the output of 
#'        bed_file_filter function.
#' @param file the name of the output file without the extension 
#' @param path The path where the function will write the file
#'
#' @return output of write function
#'
#' @examples
#' \dontrun{
#'        require(EnsDb.Hsapiens.v86)
#'        edb <- EnsDb.Hsapiens.v86
#'        quantification_files <- 'file_path'
#'         ebgwot <- exon_by_gene_with_observed_transcripts(edb, 
#'                                                        quantification_files)
#'         bed_file_content_gr <- GRanges("chr16",ranges = IRanges(start=23581002, 
#'                                                                end=23596356))
#'         bffr <- bed_file_filter(ebgwot, bed_file_content_gr)
#'         write_bed_file_filter_result(bffr, file='test','./')
#' }
#'
write_bed_file_filter_result <- function(bed_file_filter_result, 
                                            file = 'file_name', path = './'){
    string <- ''
    bffr <- bed_file_filter_result
    for (i in 1:length(bffr)){
        string <- paste0(string, as.vector(seqnames(bffr))[i],'\t',
                            start(bffr)[i],'\t',
                            end(bffr)[i],'\t',
                            '.\t.\t',
                            as.vector(strand(bffr))[i],
                            '\n')
    }
    write(string, paste0(path, file, '.BED'))
}
