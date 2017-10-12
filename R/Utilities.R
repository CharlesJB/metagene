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
#'    require(EnsDb.Hsapiens.v86)
#'    edb <- EnsDb.Hsapiens.v86
#'    quantification_files <- 'file_path'
#'    ebgwot <- exon_by_gene_with_observed_transcripts(edb, 
#'                                                     quantification_files)
#'    bed_file_content_gr <- GRanges("chr16",ranges = IRanges(start=23581002, 
#'                                                            end=23596356))
#'    bed_file_filter(ebgwot, bed_file_content_gr)
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
#'   require(EnsDb.Hsapiens.v86)
#'   edb <- EnsDb.Hsapiens.v86
#'   quantification_files <- 'file_path'
#'   ebgwot <- exon_by_gene_with_observed_transcripts(edb, 
#'                                                    quantification_files)
#'   bed_file_content_gr <- GRanges("chr16",ranges = IRanges(start=23581002, 
#'                                                            end=23596356))
#'   bffr <- bed_file_filter(ebgwot, bed_file_content_gr)
#'   write_bed_file_filter_result(bffr, file='test','./')
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

#' Is is a function designed to remove values <= to 'gaps_threshold'. 
#' Nucleotides local and global positions, bins, size of regions/genes and exons
#' will be recalculated. To use on metagene's table during RNA-seq analysis. 
#' Not made for ChIP-Seq analysis or to apply on matagene's data_frame. A
#' similar function is implemented in produce_data_frame() with same arguments.
#' The unique goal of this function is to allow permutation_test which match
#' the plot created using avoid_gaps, bam_name and gaps_threshold arguments
#' in the produce_data_frame function.
#'
#' @param table A data.table from produce_table(...) function of metagene.
#' @param bam_name A reference bam_name to allow the same removal 
#'                 (position in bam) of values for other bam file.
#' @param gaps_threshold A threshold under which values will be removed.    
#'
#' @return A data.table with values <= to 'gaps_threshold' removed
#'
#' @examples
#' \dontrun{
#'  bam_files <- c(
#'  system.file("extdata/c_al4_945MLM_demo_sorted.bam", package="metagene"),
#'  system.file("extdata/c_al3_362PYX_demo_sorted.bam", package="metagene"),
#'  system.file("extdata/n_al4_310HII_demo_sorted.bam", package="metagene"),
#'  system.file("extdata/n_al3_588WMR_demo_sorted.bam", package="metagene"))
#'  region <- c(
#'  system.file("extdata/ENCFF355RXX_DPM1less.bed", package="metagene"),
#'  system.file("extdata/ENCFF355RXX_NDUFAB1less.bed", package="metagene"),
#'  system.file("extdata/ENCFF355RXX_SLC25A5less.bed", package="metagene"))
#'  mydesign <- matrix(c(1,1,0,0,0,0,1,1),ncol=2, byrow=FALSE)
#'  mydesign <- cbind(c("c_al4_945MLM_demo_sorted.bam",
#'                      "c_al3_362PYX_demo_sorted.bam",
#'                      "n_al4_310HII_demo_sorted.bam",
#'                      "n_al3_588WMR_demo_sorted.bam"), mydesign)
#'  colnames(mydesign) <- c('Samples', 'cyto', 'nucleo')
#'  mydesign <- data.frame(mydesign)
#'  mydesign[,2] <- as.numeric(mydesign[,2])-1
#'  mydesign[,3] <- as.numeric(mydesign[,3])-1
#'
#'  mg <- metagene$new(regions = region, bam_files = bam_files, 
#'                                                            assay = 'rnaseq')
#'  mg$produce_table(flip_regions = FALSE, bin_count = 100, 
#'                                design = mydesign, normalization = 'RPM')
#'  mg$produce_data_frame(avoid_gaps = TRUE, 
#'                        bam_name = "c_al4_945MLM_demo_sorted", 
#'                        gaps_threshold = 10)
#'  mg$plot()
#'  tab <- mg$get_table()
#'  tab <- avoid_gaps_update(tab, 
#'         bam_name = 'c_al4_945MLM_demo_sorted', gaps_threshold = 10)
#'  tab0 <- mg$get_table()
#'  tab1 <- tab0[which(tab0$design == "cyto"),]
#'  tab2 <- tab0[which(tab0$design == "nucleo"),]
#'  
#'  library(similaRpeak)
#'  perm_fun <- function(profile1, profile2) {
#'    sim <- similarity(profile1, profile2)
#'    sim[["metrics"]][["RATIO_NORMALIZED_INTERSECT"]]
#'  }
#'  
#'  ratio_normalized_intersect <- 
#'    perm_fun(tab1[, .(moy=mean(value)), by=bin]$moy, 
#'             tab2[, .(moy=mean(value)), by=bin]$moy)
#'  ratio_normalized_intersect
#'  
#'  permutation_results <- permutation_test(tab1, tab2, sample_size = 2,
#'                                  sample_count = 1000, FUN = perm_fun)
#'  hist(permutation_results, 
#'           main="ratio_normalized_intersect (1=total overlapping area)")
#'  abline(v=ratio_normalized_intersect, col = 'red')
#'  sum(ratio_normalized_intersect >= permutation_results) / 
#'         length(permutation_results)
#' }
#'

avoid_gaps_update <- function(table, bam_name, gaps_threshold = 0){
    new_table <- data.table::copy(table)
    
    bin_count <- max(unique(new_table$bin))
    message(paste('Gaps deletion is calibrated on data from',
                    'the bam file name provided as argument "bam_name" in',
                    '"produce_data_frame()" method. Otherwise, the first bam',
                    'will be used as default'))
            
    #how_namy_by_exon_by_design
    nb_nuc_removed <- new_table[value <= gaps_threshold 
                            & bam == bam_name, length(value),
                        by=c('exon', 'region')]
    
    #assignment of new exonsize
    for (i in 1:length(nb_nuc_removed$V1)){
        #selected = lines of the ith region and exon of nb_nuc_removed
        selected <- which(
            new_table$region == nb_nuc_removed$region[i] &
            new_table$exon == nb_nuc_removed$exon[i])
        #retrieve the exonsize value of the ith region and exon
        original_exonsize <- unique(new_table$exonsize[selected])
        #replace former exonsixe
        new_exonsize <- original_exonsize-nb_nuc_removed$V1[i]
        new_table$exonsize[selected] <- new_exonsize
    }
    
    nb_nuc_removed_by_gene <- new_table[value <= gaps_threshold 
                            & bam == bam_name, length(value),
                        by=c('region')]
    #assignment of new regionsize/genesize
    for (i in 1:length(unique(nb_nuc_removed_by_gene$region))){
        #selected = lines of the ith region of nb_nuc_removed
        selected <- which(
            new_table$region == nb_nuc_removed_by_gene$region[i])
        #retrieve the regionsize value of the ith region and exon
        original_regionsize <- unique(new_table$regionsize[selected])
        #replace former regionsize
        new_regionsize <- (original_regionsize
                            - nb_nuc_removed_by_gene$V1[i])
        new_table$regionsize[selected] <- new_regionsize
    }
    
	#assignment of new regionstartnuc (allow flip nuctot after gaps removal)
	for (i in 1:length(unique(nb_nuc_removed_by_gene$region))){
        #selected = lines of the ith region of nb_nuc_removed
        selected <- which(
            new_table$region == nb_nuc_removed_by_gene$region[i])
        #retrieve the regionsize value of the ith region and exon
        original_regionstartnuc <- unique(new_table$regionstartnuc[selected])
        #replace former regionsize
        new_regionstartnuc <- (original_regionstartnuc
                            - nb_nuc_removed_by_gene$V1[i])
        new_table$regionstartnuc[selected] <- new_regionstartnuc
    }
	
    ### removal of zero values
    ## stop if all bam haven't the same amount of lines in table
    stopifnot(length(unique(new_table[, .N, by=bam]$N)) == 1)
    bam_line_count <- tab[bam == bam_name, .N]
    #lines_to_remove for bam_name
    lines_to_remove <- which(new_table$bam == bam_name &
                                    new_table$value <= gaps_threshold)
    # %% provide the idx for the first bam
    lines_to_remove <- (lines_to_remove %% bam_line_count)
    #to avoid 0 if there is a x %% x = 0
    lines_to_remove <- replace(lines_to_remove, 
                            which(lines_to_remove == 0), 
                            bam_line_count)
    bam_count <- length(unique(new_table$bam))
    #lines_to_remove for all bam
    lines_to_remove <- unlist(map((0:(bam_count-1)), 
                        ~ lines_to_remove + bam_line_count * .x))
    new_table <- new_table[-lines_to_remove,]

    
    #reinitialization of nuctot before flip in next section to 
    # clear gaps in nuctot number seauence
    new_table$nuctot <- rep(1:length(which(
                    new_table$bam == bam_name)),
                    times = length(unique(new_table$bam)))
    
    #reorder the nuc and nuctot variables
    ascending = function(nuc) {
            nuc[1] < nuc[2]
    }
    are_genes_unflipped <- unlist(lapply(map(unique(new_table$region),
                        ~ new_table[which(new_table$region == .x & 
                                    new_table$bam == new_table$bam[1]),]$nuc)
                                        , ascending))
    if(!all(are_genes_unflipped)){
        
        flip_by_bam_n_region <- map2(rep(unique(new_table$bam), 
                    each=length(unique(new_table$region))), 
            rep(unique(new_table$region),
                    times=length(unique(new_table$bam))), 
            ~which(new_table$bam == .x & new_table$region == .y 
                        & new_table$strand == '-'))
        
        not_empty_idx <- which(map(flip_by_bam_n_region, 
                                                ~length(.x)) > 0) 
        if (length(not_empty_idx) > 0){
            map(flip_by_bam_n_region[not_empty_idx],
                            ~ (new_table$nuc[.x] <- length(.x):1))
            map(flip_by_bam_n_region[not_empty_idx],
                            ~ (new_table$nuctot[.x] <- 
                                    max(new_table$nuctot[.x]):
                                        min(new_table$nuctot[.x])))
        }
        
        unflip_by_bam_n_region <- map2(rep(unique(new_table$bam), 
                    each=length(unique(new_table$region))), 
            rep(unique(new_table$region), 
                    times=length(unique(new_table$bam))), 
            ~which(new_table$bam == .x & new_table$region == .y 
                        & (new_table$strand == '+' | 
                            new_table$strand == '*')))
        not_empty_idx <- which(map(unflip_by_bam_n_region, 
                                                ~length(.x)) > 0) 
        if (length(not_empty_idx) > 0){
            map(unflip_by_bam_n_region[not_empty_idx],
                            ~ (new_table$nuc[.x] <- 1:length(.x)))
            map(unflip_by_bam_n_region[not_empty_idx],
                            ~ (new_table$nuctot[.x] <- 
                                    min(new_table$nuctot[.x]):
                                        max(new_table$nuctot[.x])))
        }
    } else { # if all(are_genes_unflipped)
        by_bam_n_region <- map2(rep(unique(new_table$bam), 
                    each=length(unique(new_table$region))), 
            rep(unique(new_table$region), 
                    times=length(unique(new_table$bam))), 
            ~which(new_table$bam == .x & new_table$region == .y))
        not_empty_idx <- which(map(by_bam_n_region, 
                                                ~length(.x)) > 0) 
        if (length(not_empty_idx) > 0){
            map(by_bam_n_region[not_empty_idx], 
                            ~ (new_table$nuc[.x] <- 1:length(.x)))
            map(by_bam_n_region[not_empty_idx], 
                            ~ (new_table$nuctot[.x] <- 
                                    min(new_table$nuctot[.x]):
                                        max(new_table$nuctot[.x])))
        }
    }
    if(!is.null(bin_count)){
        #reinitialization of region/gene_size to be able to rebuild 
        #bin column
        length_by_region_n_bam <- new_table[,length(nuc),
                                            by=c('region','bam')]$V1
        new_table$regionsize <- rep(length_by_region_n_bam, 
                                    times=length_by_region_n_bam)
        #rebuild the correct bin column
        col_bins <- trunc((new_table$nuc/(new_table$regionsize+1))
                                *bin_count)+1
        new_table$bin <- as.integer(col_bins)
    }
    return(new_table)
}