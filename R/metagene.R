#' A class to manage metagene analysis.
#'
#' This class will allow to load, convert and normalize alignments and regions
#' files/data. Once the data is ready, the user can then chose to produce
#' metagene plots on the data (or a subset of the data).
#'
#' @section Constructor:
#' \describe{
#'    \item{}{\code{mg <- metagene$new(regions, bam_files, padding_size = 0,
#'                            cores = SerialParam(), verbose = FALSE,
#'                            force_seqlevels = FALSE, paired_end = FALSE,
#'                            assay = 'chipseq'))}}
#'    \item{regions}{Either a \code{vector} of BED, narrowPeak or broadPeak
#'                    filenames, a \code{GRanges} object or a 
#'                    \code{GRangesList} object.}
#'    \item{bam_files}{A \code{vector} of BAM filenames. The BAM files must be
#'                    indexed. i.e.: if a file is named file.bam, there must
#'                    be a file named file.bam.bai or file.bai in the same 
#'                    directory.}
#'    \item{padding_size}{The regions will be extended on each side by the
#'                        value of this parameter. The padding_size must be a
#'                        non-negative integer. Default = 0.}
#'    \item{cores}{The number of cores available to parallelize the analysis.
#'                Either a positive integer or a \code{BiocParallelParam}.
#'                Default: \code{SerialParam()}.}
#'    \item{verbose}{Print progression of the analysis. A logical constant.
#'                    Default: \code{FALSE}.}
#'    \item{force_seqlevels}{If \code{TRUE}, Remove regions that are not found
#'                in bam file header. Default: \code{FALSE}. TRUE and FALSE
#'                respectively correspond to pruning.mode = "coarse"
#'                and "error" in ?seqinfo.}
#'    \item{paired_end}{If \code{TRUE}, metagene will deal with paired-ended 
#'                data. If \code{FALSE}, single-ended data are expected. 
#'                Default: \code{FALSE}}
#'    \item{assay}{\code{'chipseq'} or \code{'rnaseq'}, the two available 
#'                options. Default: \code{'chipseq'}}
#' }
#'
#'    \code{metagene$new} returns a \code{metagene} object that contains the
#'        coverages for every BAM files in the regions from the \code{regions}
#'        param.
#'
#' @return
#' \code{metagene$new} returns a \code{metagene} object which contains the
#' normalized coverage values for every regions and for every BAM files.
#'
#' @section Methods:
#' \describe{
#'    \item{}{\code{mg$plot(region_names = NULL, design_names = NULL,
#'                title = NULL, x_label = NULL)}}
#'    \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                        all the regions are returned. Default: \code{NULL}.}
#'    \item{design_names}{The names of the experiments to extract. If a design
#'            was added to the \code{metagene} object, \code{design_names}
#'            correspond to the column names in the design, otherwise
#'            \code{design_names} corresponds to the BAM name or the BAM
#'            filename. If \code{NULL}, all the experiments are
#'            returned. Default: \code{NULL}.}
#'    \item{title}{A title to add to the graph. If \code{NULL}, will be
#'                        automatically created. Default: NULL}
#'    \item{x_label}{X-axis label to add to the metagene plot. If \code{NULL},
#'                    metagene will use generic label. Default: \code{NULL}.}
#' }
#' \describe{
#'    \item{}{\code{mg$produce_table(design, bin_count, noise_removal,
#'                normalization, flip_regions, bin_size = NULL}}
#'    \item{design}{A \code{data.frame} that describe to experiment to plot.
#'            see \code{plot} function for more details. \code{NA} can 
#'            be used keep previous design value. Default: \code{NA}.}
#'    \item{bin_count}{The number of bin to create. \code{NA} can be used to
#'                    keep previous bin_count value. A bin_count value of 100
#'                    will be used if no value is specified. Default:
#'                    \code{NA}.}
#'    \item{noise_removal}{The algorithm to use to remove control(s). Possible
#'                        values are \code{NA}, \code{NULL} or "NCIS". By
#'                        default, value is \code{NULL}. Use \code{NA} keep
#'                        previous \code{noise_removal} value (i.e. if
#'                        \code{produce_table} was called before). See
#'                        Liand and Keles 2012 for the NCIS algorithm.}
#'    \item{normalization}{The algorithm to use to normalize samples. Possible
#;                        values are \code{NA}, \code{NULL} or "RPM". By
#'                        default, value is \code{NULL} and no normalization
#'                        will be performed. Use \code{NA} keep
#'                        previous \code{normalization} value (i.e. if
#'                        \code{produce_table} was called before).}
#'    \item{flip_regions}{Should regions on negative strand be flip_regions?
#'                        Default: \code{FALSE}.}
#'    \item{bin_size}{Deprecated.}
#' }
#' \describe{
#'    \item{}{\code{mg$produce_data_frame(alpha = 0.05, sample_count = 1000, 
#'                                avoid_gaps = FALSE, gaps_threshold = 0)}}
#'    \item{alpha}{The range of the estimation to be shown with the ribbon.
#'                \code{1 - alpha / 2} and \code{alpha / 2} will be used.
#'                Default: 0.05.}
#'    \item{sample_count}{The number of draw to do in the bootstrap
#'                        calculation. Default: 1000.}
#'    \item{avoid_gaps}{Provide the possibility to remove values = 0 and refit
#'                    the data_frame for this suppression.
#'                    Default : \code{FALSE}.}
#'    \item{gaps_threshold}{It works with avoid_gaps argument. It lets to remove
#'                    values <= at gaps_threshold. Default : 0.}
#' }
#' \describe{
#'    \item{}{mg$get_params()}
#' }
#' \describe{
#'    \item{}{mg$get_design()}
#' }
#' \describe{
#'    \item{}{mg$get_regions(region_names = NULL)}
#'    \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                        all the regions are returned. Default: \code{NULL}.}
#' }
#' \describe{
#'    \item{}{mg$get_table = function()}
#' }
#' \describe{
#'    \item{}{mg$get_matrices = function()}
#' }
#' \describe{
#'    \item{}{mg$get_data_frame(region_names = NULL, design_names = NULL)}
#'    \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                        all the regions are returned. Default: \code{NULL}.}
#'    \item{design_names}{The names of the experiments to extract. If a design
#'            was added to the \code{metagene} object, \code{design_names}
#'            correspond to the column names in the design, otherwise
#'            \code{design_names} corresponds to the BAM name or the BAM
#'            filename. If \code{NULL}, all the experiments are
#'            returned. Default: \code{NULL}.}
#' }
#' \describe{
#'    \item{}{get_plot = function()}
#' }
#' \describe{
#'    \item{}{get_raw_coverages = function(filenames)}
#'    \item{filenames}{The name of the file to extract raw coverages. Can be
#'                    the filename with the extension of the name of the bam
#'                    file (if a named bam files was used during the creation
#'                    of the metagene object). If \code{NULL}, returns the
#'                    coverage of every bam files. Default: \code{NULL}.}
#' }
#' \describe{
#'    \item{}{get_normalized_coverages = function(filenames)}
#'    \item{filenames}{The name of the file to extract normalized coverages 
#'            (in RPM). Can be the filename with the extension of 
#'            the name of the bam file (if a named bam files was used during
#'            the creation of the metagene object). If \code{NULL},
#'            returns the coverage every bam files. Default:
#'            \code{NULL}.}
#' }
#' \describe{
#'    \item{}{\code{mg$export(bam_file, region, file)}}
#'    \item{bam_file}{The name of the bam file to export.}
#'    \item{region}{The name of the region to export.}
#'    \item{file}{The name of the ouput file.}
#' }
#' \describe{
#'    \item{}{\code{mg$add_design(design = NULL, check_bam_files = FALSE)}}
#'    \item{design}{A \code{data.frame} that describe to experiment to plot.
#'            See \code{plot} function for more details. \code{NA} can be
#'            used keep previous design value. Default: \code{NA}.}
#'    \item{check_bam_files}{Force check that all the bam files from the first
#'                            columns of the design are present in current
#'                            metagene object. Default: \code{FALSE}}
#' }
#' \describe{
#'    \item{}{\code{mg$unflip_regions()}}
#' }
#'
#' \describe{
#'    \item{}{\code{mg$flip_regions()}}
#' }
#' \describe{
#'    \item{}{\code{mg$unflip_regions()}}
#' }
#'
#' @examples
#' region <- get_demo_regions()[1]
#' bam_file <- get_demo_bam_files()[1]
#' mg <- metagene$new(regions = region, bam_files = bam_file)
#' \dontrun{
#'    df <- metagene$plot()
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom data.table data.table
#' @export
#' @format A metagene experiment manager

metagene <- R6Class("metagene",
    public = list(
    # Methods
        initialize = function(regions, bam_files, padding_size = 0,
                                cores = SerialParam(), verbose = FALSE,
                                force_seqlevels = FALSE, paired_end = FALSE,
                                assay = 'chipseq') {
            # Check params...
            private$check_param(regions = regions, bam_files = bam_files,
                                padding_size = padding_size,
                                cores = cores, verbose = verbose,
                                force_seqlevels = force_seqlevels, 
                                assay = assay)

            # Save params
            private$parallel_job <- Parallel_Job$new(cores)
            private$params[["padding_size"]] <- padding_size
            private$params[["verbose"]] <- verbose
            if (is.null(names(bam_files))) {
                new_names <- tools::file_path_sans_ext(basename(bam_files))
                names(bam_files) <- new_names
            }
            private$params[["bin_count"]] <- NULL
            private$params[["bam_files"]] <- bam_files
            private$params[["force_seqlevels"]] <- force_seqlevels
            private$params[["flip_regions"]] <- FALSE
            private$params[["assay"]] <- tolower(assay)
            private$params[["df_needs_update"]] <- TRUE
            private$params[["df_arguments"]] <- ""
            private$params[["table_needs_update"]] <- TRUE
            
            # Prepare bam files
            private$print_verbose("Prepare bam files...")
            private$bam_handler <- Bam_Handler$new(bam_files, 
                                        cores = cores,
                                        paired_end = paired_end)

            # Prepare regions
            private$print_verbose("Prepare regions...")
            private$regions <- private$prepare_regions(regions)

            # Parse bam files
            private$print_verbose("Parse bam files...\n")
            private$print_verbose("coverages...\n")
            private$coverages <- private$produce_coverages()        
        },
        get_bam_count = function(filename) {
            # Parameters validation are done by Bam_Handler object
            private$bam_handler$get_aligned_count(filename)
        },
        get_params = function() {
            private$params
        },
        get_design = function() {
            private$design
        },
        get_regions = function(region_names = NULL) {
            if (is.null(region_names)) {
                private$regions
            } else {
                new_names <- tools::file_path_sans_ext(
                                            basename(region_names))
                region_names <- new_names
                stopifnot(all(region_names %in% names(private$regions)))
                private$regions[region_names]
            }
        },
        get_table = function() {
            if (length(private$table) == 0) { 
                return(NULL)
            }
            if (private$params[['assay']] == 'chipseq'){
                return(copy(private$table))
            } else if (private$params[['assay']] == 'rnaseq'){
                if('bin' %in% colnames(private$table)){
                    returnedCol <- c("region","exon","bam","design","nuc","bin",
                                    "nuctot","exonsize","regionstartnuc",
                                    "regionsize","value","strand")
                } else {
                    returnedCol <- c("region","exon","bam","design","nuc",
                                    "nuctot","exonsize","regionstartnuc",
                                    "regionsize","value","strand")
                }
                return(copy(private$table[,returnedCol,with=FALSE]))
            }
        },
        get_matrices = function() {
            if (is.null(self$get_table())){
                return(NULL)
            }
            if (private$params[['assay']] == 'chipseq') {
                matrices <- list()
                nbcol <- private$params[["bin_count"]]
                nbrow <- vapply(self$get_regions(), length, numeric(1))
                for (regions in names(self$get_regions())) {
                    matrices[[regions]] <- list()
                    for (design_name in colnames(self$get_design())[-1]) {
                        matrices[[regions]][[design_name]] <- list()
                        matrices[[regions]][[design_name]][["input"]] <- 
                                matrix(private$table[region == regions & 
                                design == design_name,]$value, 
                                nrow=nbrow[regions], ncol=nbcol, byrow=TRUE)
                    }
                }
                return (matrices)
            } else {
                stop(paste('unsupported function for assay of type',
                        private$params[['assay']],
                        '. Only available for "chipseq" assay.')) 
            }
        },
        get_data_frame = function(region_names = NULL, design_names = NULL) {
            if (nrow(private$df) == 0) {
                NULL
            } else if (is.null(region_names) & is.null(design_names)) {
                return(copy(private$df))
            } else {
                if (!is.null(region_names)) {
                    stopifnot(is.character(region_names))
                    stopifnot(all(region_names %in% 
                                    unique(private$table$region)))
                } else {
                    region_names <- names(private$regions)
                }
                if (!is.null(design_names)) {
                    stopifnot(is.character(design_names))
                    stopifnot(all(design_names %in% 
                                    unique(private$table$design)))
                } else {
                    design_names <- colnames(private$design)[-1]
                }
                i <- (private$df$region %in% region_names &
                                    private$df$design %in% design_names)
                return(copy(private$df[i,]))
            }
        },
        get_plot = function() {
            if (is.character(private$graph)) {
                NULL
            } else {
                private$graph
            }
        },
        get_raw_coverages = function(filenames = NULL) {
            if (is.null(filenames)) {
                private$coverages
            } else {
                stopifnot(is.character(filenames))
                stopifnot(length(filenames) > 0)
                bam_names <- private$get_bam_names(filenames)
                stopifnot(length(bam_names) == length(filenames))
                private$coverages[bam_names]
            }
        },
        get_normalized_coverages = function(filenames = NULL) {
            normalize_coverage <- function(filename) {
            count <- private$bam_handler$get_aligned_count(filename)
                weight <- 1 / (count / 1000000)
                coverages[[filename]] <- coverages[[filename]] * weight
        }
        coverages <- self$get_raw_coverages(filenames)
            coverage_names <- names(coverages)
            coverages <-
                private$parallel_job$launch_job(data = coverage_names,
                                                FUN = normalize_coverage)
            names(coverages) <- coverage_names
            coverages
        },
        add_design = function(design, check_bam_files = FALSE) {
            private$design = private$fetch_design(design, check_bam_files)
            private$params[["table_needs_update"]] <- TRUE
        },
        produce_table = function(design = NA, bin_count = NA, bin_size = NULL,
                                noise_removal = NA, normalization = NA,
                                flip_regions = FALSE) {
            if (!is.null(bin_size)) {
                warning("bin_size is now deprecated. Please use bin_count.")
                bin_size <- NULL
            }

            design <- private$fetch_design(design)
            private$check_produce_table_params(bin_count = bin_count,
                                                bin_size = bin_size,
                                                design = design,
                                                noise_removal = noise_removal,
                                                normalization = normalization,
                                                flip_regions = flip_regions)
            bin_count <- private$get_param_value(bin_count, "bin_count")
            noise_removal <- private$get_param_value(noise_removal,
                                                    "noise_removal")
            normalization <- private$get_param_value(normalization,
                                                    "normalization")
            coverages <- private$coverages

            #addition of private$params[["table_needs_update"]] comes from 
            #troubles in table update when adding a design with the add_design
            #function and changing nothing else in produce_table parameters
            if (private$table_need_update(design = design,
                                    bin_count = bin_count,
                                    bin_size = bin_size,
                                    noise_removal = noise_removal,
                                    normalization = normalization) |
                                    private$params[["table_needs_update"]]) {
                
                if (!is.null(normalization)) {
                    coverages <- private$normalize_coverages(coverages)
                    message('Normalization done')
                }
                
                if (private$params[['assay']] == 'rnaseq'){
                
                    # here the word 'gene' = 'region'
                    bam_files_names <- names(private$params[["bam_files"]])

                    # useful variables for caculations of 
                    # standard data table structure (std_dt_struct)
                    gene_count <- length(private$regions)
                    gene_names <- names(private$regions)
                    exon_length_by_exon_by_gene <- width(private$regions)
                    exon_length_by_exon_by_gene_cum <- vapply(
                                exon_length_by_exon_by_gene, sum, numeric(1))

                    exon_count_by_gene <- vapply(private$regions, length, 
                                                                numeric(1))
                    
                    ## standard data table structure for one replicat
                        col_gene <- rep(gene_names, times=unlist(map(
                                1:gene_count, 
                                ~sum(exon_length_by_exon_by_gene[[
                                                            gene_names[.x]]]))))

                        exon_names <- unlist(map(exon_count_by_gene, ~ 1:.x))
                        col_exon <- as.vector(rep(exon_names, 
                                    times=unlist(exon_length_by_exon_by_gene)))
                                    
                        col_exon_size <- rep(as.vector(#useful for flip function
                                unlist(exon_length_by_exon_by_gene)), 
                                times=as.vector(unlist(
                                    exon_length_by_exon_by_gene)))
                                    
                        gene_size <- unlist(map(1:length(self$get_regions()), 
                                        ~sum(width(self$get_regions()[[.x]]))))
                        col_gene_size <- rep(gene_size, times = gene_size)
                        
                        gene_length_cum <- c(0,
                                    cumsum(gene_size)[-length(gene_size)])+1
                        col_gene_start_nuc <- rep(gene_length_cum, 
                                                            times = gene_size)
                        
                        col_nuc <- unlist(map(as.vector(unlist(
                                    exon_length_by_exon_by_gene_cum)), ~ 1:.x))
                        
                        exon_strand_by_exon_by_gene <- unlist(map(gene_names, 
                                    ~as.vector(strand(private$regions[[.x]]))))
                        col_strand <- as.vector(rep(exon_strand_by_exon_by_gene,
                                    times=unlist(exon_length_by_exon_by_gene)))
                    
                        
                        length_std_dt_struct = length(col_gene)

                        # the number of not empty cases in the design param
                        copies_count <- sum(replace(unlist(design[,-1]),
                                    which(unlist(design[,-1]) == 2),1))

                        #multiplication of standard data table structure
                        col_gene <- rep(col_gene, copies_count)
                        col_exon <- rep(col_exon, copies_count)
                        col_nuctot <- 1:length(col_nuc)
                        col_nuctot <- rep(col_nuctot, copies_count)
                        col_nuc <- rep(col_nuc, copies_count)
                        col_exon_size <- rep(col_exon_size, copies_count)
                        col_gene_size <- rep(col_gene_size, copies_count)
                        col_strand <- rep(col_strand, copies_count)
                        col_gene_start_nuc <- rep(col_gene_start_nuc, 
                                                            copies_count)
                        
                    ## other columns of data table
                        design_names <- colnames(design)[-1]
                        bam_names_in_design <- tools::file_path_sans_ext(
                                                            basename(design[,1]))

                        bfile_names_by_design <- tools::file_path_sans_ext(
                            unlist(map(design_names , 
                                ~design[which(design[,
                                    which(colnames(
                                        design) == .x)] > 0),1])))
                        col_bam <- rep(bfile_names_by_design,
                                each=length_std_dt_struct)
                        
                        nb_bfile_by_design <- unlist(map(design_names , 
                                ~length(which(design[,which(
                                colnames(design) == .x)] > 0))))
                        col_design <- rep(design_names,
                                    times=(nb_bfile_by_design *
                                        length_std_dt_struct))
                        
                        ## col_values
                        #NB : lapply(Views...) -> out of limits of view
                        grtot <- self$get_regions()
                        col_values <- list()
                        idx <- 1 #index for col_values list
                        idx_sd_loop <- 1 
                        for(bam in bam_names_in_design) {
                            for (i in 1:length(grtot)){
								gr <- grtot[[i]]
								sq <- unique(as.character(seqnames(gr)))
                                val <- Views(
                                    coverages[[bam]][[sq]], 
                                    start(gr), 
                                    end(gr))
                                col_values[[idx]] <- unlist(lapply(
                                val, as.numeric))
                                idx <- idx + 1
                            }
                        }
                        col_values <- unlist(col_values)
                    
                    if (!is.null(bin_count)) {
                        message('produce data table : RNA-Seq binned')
                        col_bins <- trunc(
                                    (col_nuc/(col_gene_size+1))*bin_count)+1
                        col_bins <- as.integer(col_bins)
                        
                        private$table <- data.table(
                                region = col_gene,
                                exon = col_exon,
                                bam = col_bam,
                                design = col_design,
                                bin = col_bins,
                                nuc = col_nuc,
                                nuctot = col_nuctot,
                                exonsize = col_exon_size,
                                regionsize = col_gene_size,
                                regionstartnuc = col_gene_start_nuc,
                                value = col_values,
                                strand = col_strand)
                    } else {
                        message('produce data table : RNA-Seq')
                        private$table <- data.table(
                                region = col_gene,
                                exon = col_exon,
                                bam = col_bam,
                                design = col_design,
                                nuc = col_nuc,
                                nuctot = col_nuctot,
                                exonsize = col_exon_size,
                                regionsize = col_gene_size,
                                regionstartnuc = col_gene_start_nuc,
                                value = col_values,
                                strand = col_strand)
                    }    
                } else { # chipseq
                    if (!is.null(noise_removal)) {
                        coverages <- private$remove_controls(coverages, design)
                    } else {
                        coverages <- private$merge_chip(coverages, design)
                    }
                
                    message('produce data table : ChIP-Seq')
                    if (is.null(bin_count)) {
                        bin_count = 100
                    }
                    region_length <- vapply(self$get_regions(), length, 
                        numeric(1))
                    col_regions <- names(self$get_regions()) %>%
                        map(~ rep(.x, length(coverages) * bin_count * 
                                    region_length[.x])) %>% unlist()
                    col_designs <- map(region_length, ~ rep(names(coverages), 
                                        each = bin_count * .x)) %>% unlist
                    col_bins <- rep(1:bin_count,
                                    length(coverages) * sum(region_length))
                    pairs <- expand.grid(colnames(design)[-1], 
                                        names(self$get_regions()), 
                                        stringsAsFactors = FALSE)
                    col_values <- map2(pairs$Var1, pairs$Var2,
                        ~ private$get_subtable(coverages[[.x]], .y, 
                            bin_count)) %>% unlist
                    
                    #TODO : improve col_strand production
                    col_strand <- list()
                    idx <- 1
                    for (region_names in unique(col_regions)){
                        col_strand[[idx]] <- rep(rep(
                            as.vector(strand(private$regions)[[region_names]]),
                            each=bin_count),length(unique(col_designs)))
                        idx <- idx + 1
                    }
                    col_strand <- unlist(col_strand)
                    
                    private$table <- data.table(region = col_regions,
                                design = col_designs,
                                bin = col_bins,
                                value = col_values,
                                strand = col_strand)
                }

                private$params[["bin_size"]] <- bin_size
                private$params[["bin_count"]] <- bin_count
                private$params[["noise_removal"]] <- noise_removal
                private$params[["normalization"]] <- normalization
                private$params[["df_needs_update"]] <- TRUE
                private$params[["table_needs_update"]] <- FALSE
                private$design <- design
            } else {
                message(paste('WARNING : table is unchanged regarding',
                    'design, bin_count, noise_removal, normalization.'))
            }
            if (flip_regions == TRUE) {
                self$flip_regions()
            } else {
                self$unflip_regions()
            }
            invisible(self)
        },
        produce_data_frame = function(alpha = 0.05, sample_count = 1000, 
                                                    avoid_gaps = FALSE, 
                                                    bam_name = NULL, 
                                                    gaps_threshold = 0) {
            
            #arguments checking
            stopifnot(is.numeric(alpha))
            stopifnot(is.numeric(sample_count))
            stopifnot(alpha >= 0 & alpha <= 1)
            stopifnot(sample_count > 0)
            sample_count <- as.integer(sample_count)
            stopifnot(is.logical(avoid_gaps))
            if (!is.null(bam_name)){
                stopifnot(is.character(bam_name))
                bam_name <- tools::file_path_sans_ext(basename(bam_name))
                bam_names <- tools::file_path_sans_ext(basename(
                                                private$params[["bam_files"]]))
                if (!bam_name %in% bam_names){
                    stop(paste("bam_name argument is no one of bam_names",
                                        "provided to the metagene object"))
                }
            }
            stopifnot(gaps_threshold >= 0)
            
            #add checks
            list_of_arguments <- paste(alpha,
                                    sample_count,
                                    avoid_gaps,
                                    bam_name,
                                    gaps_threshold)
                                    
            if (private$params[["df_arguments"]] != list_of_arguments){
                private$params[["df_arguments"]] <- list_of_arguments
                private$params[["df_needs_update"]] <- TRUE
            }
            
            if (private$params[['df_needs_update']]){
            
                # 1. Get the correctly formatted table
                if (is.null(self$get_table())) {
                    self$produce_table()
                }
            
                # 2. Produce the data.frame 
                private$df <- data.table::copy(self$get_table())
            
                if (private$params[['assay']] == 'chipseq') {
                    message('produce data frame : ChIP-Seq')
                    private$data_frame_need_update(alpha, sample_count)
                    sample_size <- self$get_table()[bin == 1,][
                                            ,.N, by = .(region, design)][
                                            , .(min(N))]
                    sample_size <- as.integer(sample_size)

                    out_cols <- c("value", "qinf", "qsup")
                    bootstrap <- function(df) {
                        sampling <- matrix(df$value[sample(seq_along(
                                                df$value),
                                                sample_size * sample_count,
                                                replace = TRUE)],
                                        ncol = sample_size)
                        values <- colMeans(sampling)
                        res <- quantile(values, c(alpha/2, 1-(alpha/2)))
                        res <- c(mean(df$value), res)
                        names(res) <- out_cols
                        as.list(res)
                    }
                    private$df <- private$df[, 
                                        c(out_cols) := bootstrap(.SD), 
                                            by = .(region, design, bin)]
                    private$df <- unique(private$df)
                } else if (private$params[['assay']] == 'rnaseq' 
                                & !('bin' %in% colnames(private$df))){
                    message('produce data frame : RNA-Seq')
                    private$data_frame_need_update(alpha, sample_count)

                    
                    sample_size <- self$get_table()[nuc == 1,][
                                            ,.N, by = .(region, design)][
                                            , .(min(N))]
                    sample_size <- as.integer(sample_size)
                    
                    out_cols <- c("value", "qinf", "qsup")
                    bootstrap <- function(df) {
                        sampling <- matrix(df$value[sample(
                                                seq_along(df$value),
                                                sample_size * sample_count,
                                                replace = TRUE)],
                                        ncol = sample_size)
                        values <- colMeans(sampling)
                        res <- quantile(values, c(alpha/2, 1-(alpha/2)))
                        res <- c(mean(df$value), res)
                        names(res) <- out_cols
                        as.list(res)
                    }


                    if(avoid_gaps){
                        if (!is.null(bam_name)){
                            private$data_frame_avoid_gaps_updates(bam_name,
                                                        gaps_threshold)
                        } else {
                            private$data_frame_avoid_gaps_updates(
                                                        private$df$bam[1], 
                                                        gaps_threshold)
                        }
                    }

###################                    
                    if(all(rowSums(self$get_design()[,-1, drop=FALSE]) == 1) &
                        all(colSums(self$get_design()[,-1, drop=FALSE]) == 1)){
                        private$df$qinf <- private$df$value
                        private$df$qsup <- private$df$value
                        private$df$group <- paste0(private$df$design,'_',private$df$region)
                        private$df <- data.frame(private$df)
                        #return(private$df)
                    } else {
###################
                        private$df <- private$df[, 
                                c(out_cols) := bootstrap(.SD), 
                                by = .(region, design, nuctot)]
                    }
                    #filter to avoid duplicated ligne (to reduce df dims)
                    #it does not matter concerning the plot. Plot works !
                    private$df <- private$df[which(!duplicated(paste(
                                    private$df$region,
                                    private$df$design,
                                    private$df$nuctot))),]
                } else if (private$params[['assay']] == 'rnaseq' 
                                        & ('bin' %in% colnames(private$df))){
                    message('produce data frame : RNA-Seq binned')
                    private$data_frame_need_update(alpha, sample_count)
                        
                    sample_size <- self$get_table()[bin == 1,][
                                            ,.N, by = .(design)][
                                            , .(min(N))]
                    sample_size <- as.integer(sample_size)
                    
                    out_cols <- c("value", "qinf", "qsup")
                    bootstrap <- function(df) {
                        sampling <- matrix(df$value[sample(
                                                seq_along(df$value),
                                                sample_size * sample_count,
                                                replace = TRUE)],
                                        ncol = sample_size)
                        values <- colMeans(sampling)
                        res <- quantile(values, c(alpha/2, 1-(alpha/2)))
                        res <- c(mean(df$value), res)
                        names(res) <- out_cols
                        as.list(res)
                    }
                    
                    if(avoid_gaps){
                        message('Avoiding gaps')
                        if (!is.null(bam_name)){
                            private$data_frame_avoid_gaps_updates(bam_name,
                                                        gaps_threshold)
                        } else {
                            private$data_frame_avoid_gaps_updates(
                                                        private$df$bam[1], 
                                                        gaps_threshold)
                        }
                    }

                    private$df <- private$df[, 
                                c(out_cols) := bootstrap(.SD), 
                                by = .(design, bin)]

                    # checked : ok !
                    private$df <- private$df[which(!duplicated(paste(
                                    private$df$design,
                                    private$df$bin))),]
                }
                private$df <- as.data.frame(private$df)
                #private$df$design <- as.factor(private$df$design)
                private$df$group <- paste(private$df$design,
                                private$df$region,
                                sep="_")
                private$df$group <- as.factor(private$df$group)
                
                
                private$params[["df_needs_update"]] <- FALSE
                invisible(self)
            }
        },
        plot = function(region_names = NULL, design_names = NULL, title = NULL,
                        x_label = NULL) {
            # 1. Get the correctly formatted table
            if (length(private$table) == 0) {
                self$produce_table()
            }

            # 2. Produce the data frame
            if (nrow(private$df) == 0) {
                self$produce_data_frame()
            }
            df <- self$get_data_frame(region_names = region_names,
                                    design_names = design_names)
            # 3. Produce the graph
            if (is.null(title)) {
                title <- paste(unique(private$df[["group"]]), collapse=" vs ")
            }
            p <- private$plot_graphic(df = df, title = title, 
                                        x_label = x_label)
            print(p)
            private$graph <- p
            invisible(self)
        },
        export = function(bam_file, region, file) {
            region <- tools::file_path_sans_ext(basename(region))
            region <- private$regions[[region]]
            param <- Rsamtools::ScanBamParam(which = region)
            alignments <- GenomicAlignments::readGAlignments(bam_file,
                                                                param = param)
            weight <- 1 - private$bam_handler$get_rpm_coefficient(bam_file)
            seqlevels(alignments) <- seqlevels(region)
            # TODO: don't use the weight param of coverage
            coverage <- GenomicAlignments::coverage(alignments, weight=weight)
            rtracklayer::export(coverage, file, "BED")
            invisible(coverage)
        },
        flip_regions = function() {
            if (private$params[["flip_regions"]] == FALSE) {
                private$flip_table()
                private$params[["flip_regions"]] <- TRUE
            }
            invisible(self)
        },
        unflip_regions = function() {
            if (private$params[["flip_regions"]] == TRUE) {
                private$flip_table()
                private$params[["flip_regions"]] <- FALSE
            }
            invisible(self)
        }
    ),
        private = list(
        params = list(),
        regions = GRangesList(),
        table = data.table(),
        design = data.frame(),
        coverages = list(),
        df = data.frame(),
        graph = "",
        bam_handler = "",
        parallel_job = "",
        check_param = function(regions, bam_files, padding_size,
                                cores, verbose, force_seqlevels, assay) {
            # Check parameters validity
            if (!is.character(assay)) {
                stop("verbose must be a character value")
            }
            assayTypeAuthorized <- c('chipseq', 'rnaseq')
            if (!(tolower(assay) %in% assayTypeAuthorized)) {
                stop("assay values must be one of 'chipseq' or 'rnaseq'")
            }
            if (!is.logical(verbose)) {
                stop("verbose must be a logicial value (TRUE or FALSE)")
            }
            if (!is.logical(force_seqlevels)) {
                stop(paste("force_seqlevels must be a logicial ",
                            "value (TRUE or FALSE)",sep=""))
            }
            if (!(is.numeric(padding_size) || is.integer(padding_size)) ||
                padding_size < 0 || as.integer(padding_size) != padding_size) {
                stop("padding_size must be a non-negative integer")
            }
            isBiocParallel = is(cores, "BiocParallelParam")
            isInteger = ((is.numeric(cores) || is.integer(cores)) &&
                            cores > 0 &&as.integer(cores) == cores)
            if (!isBiocParallel && !isInteger) {
                stop(paste0("cores must be a positive numeric or ",
                            "BiocParallelParam instance"))
            }
            if (!is.vector(bam_files, "character")) {
                stop("bam_files must be a vector of BAM filenames")
            }
            if (!all(sapply(bam_files, file.exists))) {
                stop("At least one BAM file does not exist")
            }
            if (!is(regions, "GRangesList") && !is.character(regions)
                && !is(regions, "GRanges") && !is.list(regions)) {
                stop(paste0("regions must be either a vector of BED ",
                    "filenames, a GRanges object or a GrangesList object"))
            }
            # Validation specific to regions as a vector
            if (is.character(regions) && !all(sapply(regions, file.exists))) {
                stop("regions must be a list of existing BED files")
            }
        },
        check_design = function(design, check_bam_files = FALSE) {
            stopifnot(is.logical(check_bam_files))
            if(!is.null(design) && !is.data.frame(design) &&
                !identical(design, NA)) {
                stop("design must be a data.frame object, NULL or NA")
            }
            if (is.data.frame(design)) {
                if(ncol(design) < 2) {
                    stop("design must have at least 2 columns")
                }
                if (!(is.character(design[,1]) || is.factor(design[,1]))) {
                    stop("The first column of design must be BAM filenames")
                }
                if (!all(apply(design[, -1, drop=FALSE], MARGIN=2,
                            is.numeric))) {
                    stop(paste0("All design column, except the first one,",
                                    " must be in numeric format"))
                }
                if (check_bam_files == TRUE) {
                    samples <- as.character(design[,1])
                    if (!all(private$check_bam_files(samples))) {
                        stop("Design contains bam files absent from metagene.")
                    }
                }
            }
        },
        check_produce_table_params = function(bin_count, bin_size, design,
                                                noise_removal, normalization,
                                                flip_regions) {
            # At least one file must be used in the design
            if (!identical(design, NA)) {
                if (!is.null(design)) {
                    if (sum(rowSums(design[ , -1, drop=FALSE]) > 0) == 0) {
                        stop(paste("At least one BAM file must be ",
                                "used in the design",sep=""))
                    }
                }
            }
            # Test only BAM file used in the design
            if (!identical(design, NA)) {
                if(!is.null(design) &&
                    !all(apply(design[rowSums(design[, -1, drop=FALSE]) > 0, 1,
                                    drop=FALSE], MARGIN = 2,
                            FUN=private$check_bam_files))) {
                    stop("At least one BAM file does not exist")
                }
            }
            # bin_count should be a numeric value without digits
            if (!identical(bin_count, NA)) {
                if (!is.null(bin_count)) {
                    if (!is.numeric(bin_count) || bin_count < 0 ||
                        as.integer(bin_count) != bin_count) {
                        stop("bin_count must be NULL or a positive integer")
                    }
                }
            }
            if (!identical(noise_removal, NA)) {
                if (!is.null(noise_removal)) {
                    if (!noise_removal %in% c("NCIS", "RPM")) {
                        msg <- 'noise_removal must be NA, NULL, "NCIS" or '
                        msg <- paste0(msg, '"RPM".')
                        stop(msg)
                    }
                }
            }
            if (!identical(normalization, NA)) {
                if (!is.null(normalization)) {
                    if (!normalization == "RPM") {
                        msg <- "normalization must be NA, NULL or \"RPM\"."
                        stop(msg)
                    }
                }
            }
            if (!is.logical(flip_regions)) {
                msg <- "flip_regions must be a logical."
                stop(msg)
            }
        },
        table_need_update = function(design, bin_count, bin_size,
                                        noise_removal, normalization) {
            if (!identical(private$design, design)) {
                return(TRUE)
            }
            if (!identical(private$params[["bin_count"]], bin_count)) {
                return(TRUE)
            }
            if (!identical(private$params[["noise_removal"]], noise_removal)) {
                return(TRUE)
            }
            if (!identical(private$params[["normalization"]], normalization)) {
                return(TRUE)
            }
            return(FALSE)
        },
        data_frame_need_update = function(alpha = NA, sample_count = NA) {
            #need_update = FALSE
            # Fetch saved values
            alpha = private$get_param_value(alpha, "alpha")
            sample_count = private$get_param_value(sample_count, 
                                                    "sample_count")

            # Add default, if needed
            if (is.null(alpha)) {
                alpha <- 0.05
            }
            if (is.null(sample_count)) {
                sample_count <- 1000
            }
            if (nrow(private$df) == 0) {
                private$params[['df_needs_update']] <- TRUE
                #private$params[["alpha"]] <- alpha
                #private$params[["sample_count"]] <- sample_count
            } else {
                # Check if data frame need update
                if (!identical(private$params[["alpha"]], alpha)) {
                    private$params[['df_needs_update']] <- TRUE
                    #private$params[["alpha"]] <- alpha
                }
                if (!identical(private$params[["sample_count"]], 
                                                    sample_count)) {
                    private$params[['df_needs_update']] <- TRUE
                    #private$params[["sample_count"]] <- sample_count
                }
            }
            #need_update
        },
        get_param_value = function(param_value, param_name) {
            param_name <- as.character(param_name)
            if (identical(param_value, NA)) {
                if (! param_name %in% names(private$params)) {
                    param_value <- NULL
                } else {
                    param_value <- private$params[[param_name]]
                }
            }
            return(param_value)
        },
        print_verbose = function(to_print) {
            if (private$params[["verbose"]]) {
                cat(paste0(to_print, "\n"))
            }
        },
        get_subtable = function(coverages, region, bcount) {
            gr <- private$regions[[region]]
            grl <- split(gr, GenomeInfoDb::seqnames(gr))
            i <- vapply(grl, length, numeric(1)) > 0
            do.call("c", lapply(grl[i], private$get_view_means,
                                bcount = bcount, cov = coverages))
        },
        get_view_means = function(gr, bcount, cov) {
            chr <- unique(as.character(GenomeInfoDb::seqnames(gr)))
            gr <- intoNbins(gr, bcount)
            stopifnot(length(chr) == 1)
            views <- Views(cov[[chr]], start(gr), end(gr))
            viewMeans(views)
        },
        prepare_regions = function(regions) {
            if (class(regions) == "character") {
                names <- sapply(regions, function(x)
                    file_path_sans_ext(basename(x)))
                import_file <- function(region) {
                    ext <- tolower(tools::file_ext(region))
                    if (ext == "narrowpeak") {
                        extraCols <- c(signalValue = "numeric",
                                        pValue = "numeric", qValue = "numeric",
                                        peak = "integer")
                        rtracklayer::import(region, format = "BED",
                                            extraCols = extraCols)
                    } else if (ext == "broadpeak") {
                        extraCols <- c(signalValue = "numeric",
                                        pValue = "numeric", qValue = "numeric")
                        rtracklayer::import(region, format = "BED",
                                            extraCols = extraCols)
                    } else if (ext == "gtf" | ext == "gff") {
                        split(rtracklayer::import(region), 
                                rtracklayer::import(region)$gene_id)
                    } else {
                        rtracklayer::import(region)
                    }
                }
                regions <- private$parallel_job$launch_job(data = regions,
                                                        FUN = import_file)
                names(regions) <- names
            } else if (class(regions) == "GRanges") {
                regions <- GRangesList("regions" = regions)
            } else if (class(regions) == "list") {
                regions <- GRangesList(regions)
            }
            if (is.null(names(regions))) {
                names(regions) <- sapply(seq_along(regions), function(x) {
                    paste("region", x, sep = "_")
                    })        
            }
            
            if (private$params[['assay']] == "rnaseq"){
                stopifnot(all(sum(width(GenomicRanges::reduce(private$regions)))
                            == sum(width(private$regions))))
            }
            # TODO: Check if there is a id column in the mcols of every ranges.
            #    If not, add one by merging seqnames, start and end.

            GRangesList(lapply(regions, function(x) {
                # Add padding
                start(x) <- start(x) - private$params$padding_size
                start(x)[start(x) < 0] <- 1
                end(x) <- end(x) + private$params$padding_size
                # Clean seqlevels
                x <- sortSeqlevels(x)
                #seqlevels(x) <- unique(as.character(seqnames(x)))
                x
            }))
        },
        produce_coverages = function() {
        regions <- GenomicRanges::reduce(BiocGenerics::unlist(private$regions))
            res <- private$parallel_job$launch_job(
                        data = private$params[["bam_files"]],
                        FUN = private$bam_handler$get_coverage,
                        regions = regions,
                        force_seqlevels= private$params[["force_seqlevels"]])
            names(res) <- names(private$params[["bam_files"]])
            lapply(res, GenomeInfoDb::sortSeqlevels)
        },
        plot_graphic = function(df, title, x_label) {
            # Prepare x label
            
            if (is.null(x_label)) {
                if (private$params[['assay']] == "chipseq") {
                    x_label <- "Distance in bins"
                } else if (private$params[['assay']] == "rnaseq" &
                            !('bin' %in% colnames(private$table))) {
                    x_label <- "Distance in nucleotides"
                } else if (private$params[['assay']] == "rnaseq" &
                            ('bin' %in% colnames(private$table))) {
                    x_label <- "Distance in bins"
                }
            }

            # Prepare y label
            y_label <- "Mean coverage"
            if (is.null(private$params[["normalization"]])) {
                y_label <- paste(y_label, "(raw)")
            } else {
                y_label <- paste(y_label, "(RPM)")
            }
            
            # Produce plot
            p <- plot_metagene(df) +
                ylab(y_label) +
                xlab(x_label) +
                ggtitle(title)
            p
        },
        fetch_design = function(design, check_bam_files = FALSE) {
            private$check_design(design = design, check_bam_files)
            get_complete_design <- function() {
                bam_files <- names(private$params[["bam_files"]])
                design <- data.frame(bam_files = bam_files)
                for (bam_file in names(private$coverages)) {
                    colname <- file_path_sans_ext(basename(bam_file))
                    design[[colname]] <- rep(0, length(bam_files))
                    i <- bam_files == bam_file
                    design[[colname]][i] <- 1
                }
                design
            }
            if (is.null(design)) {
                return(get_complete_design())
            }
            if (identical(design, NA)) {
                if (all(dim(private$design) == c(0, 0))) {
                    return(get_complete_design())
                } else {
                    return(private$design)
                }
            }
            design[,1] <- as.character(design[,1])
            return(design)
        },
        remove_controls = function(coverages, design) {
            results <- list()
            for (design_name in colnames(design)[-1]) {
                i <- design[[design_name]] == 1
                j <- design[[design_name]] == 2
            chip_bam_files <- as.character(design[,1][i])
            chip_names <- private$get_bam_names(chip_bam_files)
                input_bam_files <- as.character(design[,1][j])
            input_names <- private$get_bam_names(input_bam_files)
                chip_coverages <- coverages[chip_names]
                chip_coverages <- Reduce("+", chip_coverages)
                if (length(input_bam_files) > 0) {
                    noise_ratio <-
                        private$bam_handler$get_noise_ratio(chip_names,
                                                            input_names)
                    input_coverages <- coverages[input_names]
                    input_coverages <- noise_ratio * Reduce("+",
                                                            input_coverages)
                    results[design_name] <- chip_coverages - input_coverages
                    i <- results[[design_name]] < 0
                    results[[design_name]][i] <- 0
                } else {
                    results[design_name] <- chip_coverages
                }
            }
            results
        },
        normalize_coverages = function(coverages) {
            bam_names <- unlist(lapply(private$bam_handler$get_bam_files()$bam, 
                                private$bam_handler$get_bam_name))
            for (bam_name in bam_names) {
                #which_rows = design[[design_name]]==1
                #bam_files <- as.character(design[,1][which_rows])
                count <- private$bam_handler$get_aligned_count(bam_name)
                #count <- sum(unlist(counts))
                weight <- 1 / (count / 1000000)
                coverages[[bam_name]] <-    coverages[[bam_name]] * weight
            }
            coverages
        },
        merge_chip = function(coverages, design) {
            result <- list()
            for (design_name in colnames(design)[-1]) {
                i <- design[[design_name]] == 1
                bam_files <- as.character(design[,1][i])
            bam_names <- private$get_bam_names(bam_files)
                cov <- coverages[bam_names]
                result[[design_name]] <- Reduce("+", cov)
            }
            result
        },
        flip_table = function() {
            if(!all(private$table[,length(levels(as.factor(strand))), 
                                    by=region][,2] == 1) &
                                    private$params[["assay"]] == 'rnaseq'){
                stop(paste('Strands of exons in one gene/region',
                                'must have the same sign to be flipped.'))
            }
            if (private$params[['assay']] == 'chipseq'){
                message('ChIP-Seq flip/unflip')
                i <- which(private$table$strand == '-')
                private$table$bin[i] <- (self$get_params()$bin_count + 1) - 
                                                        private$table$bin[i]
                private$table$bin <- as.integer(private$table$bin)
                private$params[["df_needs_update"]] <- TRUE
            } else if (private$params[['assay']] == 'rnaseq'){
                message('RNA-Seq flip/unflip')
                i <- which(private$table$strand == '-')
                #col_nuc
                private$table$nuc[i] <- (private$table$regionsize[i] + 1) - 
                                                private$table$nuc[i]
                private$table$nuc <- as.integer(private$table$nuc)
                #col_nuctot
                private$table$nuctot[i] <- (private$table$regionsize[i] + 1) - 
                                        private$table$nuctot[i] + 
                                        private$table$regionstartnuc[i] * 2 - 2
                private$table$nuctot <- as.integer(private$table$nuctot)
                #col_bin
                if(!is.null(private$params[["bin_count"]])){
                    private$table$bin[i] <- (self$get_params()$bin_count + 1) - 
                                                        private$table$bin[i]
                    private$table$bin <- as.integer(private$table$bin)
                }
                private$params[["df_needs_update"]] <- TRUE
            }
        },
        get_bam_names = function(filenames) {
            if (all(filenames %in% colnames(private$design)[-1])) {
                filenames
            } else {
                stopifnot(private$check_bam_files(filenames))
                vapply(filenames,
                    private$bam_handler$get_bam_name,
                    character(1))
            }
        },
        check_bam_files = function(bam_files) {
            all(vapply(bam_files,
            function(x) {
                !is.null((private$bam_handler$get_bam_name(x)))
            },
            logical(1)))
        },
        data_frame_avoid_gaps_updates = function(bam_name, gaps_threshold) {
            #bootstrap not executed at this point. Don't work on design !
            
            #how_namy_by_exon_by_design
            dfdt <- data.table::copy(private$df)
            nb_nuc_removed <- dfdt[value <= gaps_threshold 
                                    & bam == bam_name, length(value),
                                by=c('exon', 'region')]
            
            #assignment of new exonsize
            for (i in 1:length(nb_nuc_removed$V1)){
                #selected = lines of the ith region and exon of nb_nuc_removed
                selected <- which(
                    private$df$region == nb_nuc_removed$region[i] &
                    private$df$exon == nb_nuc_removed$exon[i])
                #retrieve the exonsize value of the ith region and exon
                original_exonsize <- unique(private$df$exonsize[selected])
                #replace former exonsixe
                new_exonsize <- original_exonsize-nb_nuc_removed$V1[i]
                private$df$exonsize[selected] <- new_exonsize
            }
            
            nb_nuc_removed_by_gene <- dfdt[value <= gaps_threshold 
                                    & bam == bam_name, length(value),
                                by=c('region')]
            #assignment of new region/genesize
            for (i in 1:length(unique(nb_nuc_removed_by_gene$region))){
                #selected = lines of the ith region of nb_nuc_removed
                selected <- which(
                    private$df$region == nb_nuc_removed_by_gene$region[i])
                #retrieve the regionsize value of the ith region and exon
                original_regionsize <- unique(private$df$regionsize[selected])
                #replace former regionsize
                new_regionsize <- (original_regionsize
                                    - nb_nuc_removed_by_gene$V1[i])
                private$df$regionsize[selected] <- new_regionsize
            }
            
            ### removal of zero values
            ## stop if all bam haven't the same amount of lines in table
            stopifnot(length(unique(private$table[, .N, by=bam]$N)) == 1)
            bam_line_count <- private$table[bam == bam_name, .N]
            #lines_to_remove for bam_name
            lines_to_remove <- which(private$df$bam == bam_name &
                                            private$df$value <= gaps_threshold)
            # %% provide the idx for the first bam
            lines_to_remove <- (lines_to_remove %% bam_line_count)
            #to avoid 0 if there is a x %% x = 0
            lines_to_remove <- replace(lines_to_remove, 
                                    which(lines_to_remove == 0), 
                                    bam_line_count)
            bam_count <- length(unique(private$df$bam))
            #lines_to_remove for all bam
            lines_to_remove <- unlist(map((0:(bam_count-1)), 
                                ~ lines_to_remove + bam_line_count * .x))
            private$df <- private$df[-lines_to_remove,]

            
            #reinitialization of nuctot before flip in next section to 
            # clear gaps in nuctot number seauence
            private$df$nuctot <- rep(1:length(which(
                            private$df$bam == bam_name)),
                            times = length(unique(private$df$bam)))
            
            #reorder the nuc and nuctot variables
            if(private$params[["flip_regions"]] == TRUE){
                flip_by_bam_n_region <- map2(rep(unique(private$df$bam), 
                            each=length(unique(private$df$region))), 
                    rep(unique(private$df$region),
                            times=length(unique(private$df$bam))), 
                    ~which(private$df$bam == .x & private$df$region == .y 
                                & private$df$strand == '-'))
                
                not_empty_idx <- which(map(flip_by_bam_n_region, 
                                                        ~length(.x)) > 0) 
                if (length(not_empty_idx) > 0){
                    map(flip_by_bam_n_region[not_empty_idx],
                                    ~ (private$df$nuc[.x] <- length(.x):1))
                    map(flip_by_bam_n_region[not_empty_idx],
                                    ~ (private$df$nuctot[.x] <- 
                                            max(private$df$nuctot[.x]):
                                                min(private$df$nuctot[.x])))
                }
                
                unflip_by_bam_n_region <- map2(rep(unique(private$df$bam), 
                            each=length(unique(private$df$region))), 
                    rep(unique(private$df$region), 
                            times=length(unique(private$df$bam))), 
                    ~which(private$df$bam == .x & private$df$region == .y 
                                & (private$df$strand == '+' | 
                                    private$df$strand == '*')))
                not_empty_idx <- which(map(unflip_by_bam_n_region, 
                                                        ~length(.x)) > 0) 
                if (length(not_empty_idx) > 0){
                    map(unflip_by_bam_n_region[not_empty_idx],
                                    ~ (private$df$nuc[.x] <- 1:length(.x)))
                    map(unflip_by_bam_n_region[not_empty_idx],
                                    ~ (private$df$nuctot[.x] <- 
                                            min(private$df$nuctot[.x]):
                                                max(private$df$nuctot[.x])))
                }
            } else { # if private$params[["flip_regions"]] == FALSE
                by_bam_n_region <- map2(rep(unique(private$df$bam), 
                            each=length(unique(private$df$region))), 
                    rep(unique(private$df$region), 
                            times=length(unique(private$df$bam))), 
                    ~which(private$df$bam == .x & private$df$region == .y))
                not_empty_idx <- which(map(by_bam_n_region, 
                                                        ~length(.x)) > 0) 
                if (length(not_empty_idx) > 0){
                    map(by_bam_n_region[not_empty_idx], 
                                    ~ (private$df$nuc[.x] <- 1:length(.x)))
                    map(by_bam_n_region[not_empty_idx], 
                                    ~ (private$df$nuctot[.x] <- 
                                            min(private$df$nuctot[.x]):
                                                max(private$df$nuctot[.x])))
                }
            }
            if(!is.null(private$params[["bin_count"]])){
                #reinitialization of region/gene_size to be able to rebuild 
                #bin column
                length_by_region_n_bam <- private$df[,length(nuc),
                                                    by=c('region','bam')]$V1
                private$df$regionsize <- rep(length_by_region_n_bam, 
                                            times=length_by_region_n_bam)
                #rebuild the correct bin column
                col_bins <- trunc((private$df$nuc/(private$df$regionsize+1))
                                        *private$params[["bin_count"]])+1
                private$df$bin <- as.integer(col_bins)
            }
        }
    )
)
