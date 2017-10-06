#' Produce a metagene plot
#' 
#' @param df a \code{data.frame} obtained with the \code{get_data_frame}
#' function. Must have the following columns: "region", "design", "bin", 
#' "value", "qinf" and "qsup".
#'
#' @return A `ggplot` object.
#'
#' @examples
#' region <- get_demo_regions()[1]
#' bam_file <- get_demo_bam_files()[1]
#' mg <- metagene$new(regions = region, bam_files = bam_file)
#' mg$produce_data_frame()
#' df <- mg$get_data_frame()
#' p <- plot_metagene(df)
plot_metagene <- function(df) {
    if (('bin' %in% colnames(df)) & !('nuc' %in% colnames(df))) { # if chipseq, for instance
    
        expected_cols <- c("bin", "value", "qinf", "qsup", "group")
        df<-df[,which(colnames(df) %in% expected_cols)]
        expected_class <- c("integer", "numeric", "factor",
                                                    rep("numeric", 2))
        stopifnot(all(expected_cols %in% colnames(df)))
        stopifnot(all(vapply(df, class, character(1)) == expected_class))

        ggplot(df, aes(x=bin, y=value, ymin=qinf, ymax=qsup)) +
            geom_ribbon(aes(fill=group), alpha=0.3) +
            geom_line(aes(color=group), size=1) +
            theme(panel.grid.major = element_line()) +
            theme(panel.grid.minor = element_line()) +
            theme(panel.background = element_blank()) +
            theme(panel.background = element_rect()) +
            theme_bw(base_size = 20)
            
    } else if (('nuc' %in% colnames(df)) & !('bin' %in% colnames(df))) { 
        # if rnaseq, for instance
        ascending = function(nuc) {
            nuc[1] < nuc[2]
        }
        
        unique_with_rep = function(vect){
            uniq <- list()
            uniq[[1]] <- vect[1]
            idx <- 2
            for (i in 2:length(vect)-1) {
                if (vect[i] != vect[i+1]) {
                    uniq[[idx]] <- vect[i+1]
                    idx <- idx + 1
                }
            }
            return(unlist(uniq))
        }
        
        #elbebg = exon_length_by_exon_by_gene
        rev_if_flipped = function(elbebg, is_genes_unflipped) {
            for (i in 1:length(elbebg)) {
                #if region/gene flipped
                if (!is_genes_unflipped[[i]]) {
                    elbebg[[i]] <- rev(elbebg[[i]])
                }
            }
            return(elbebg)
        }
        
        #if only one region/gene in the data_frame
        if (length(unique(df$region)) == 1) {
            print('One region/gene')
            
            is_genes_unflipped <- unlist(lapply(map(unique(df$region), 
                                ~ df[which(df$region == .x),]$nuc), ascending))
            if (all(is_genes_unflipped)){
                exon_separation_bars <- cumsum(unique_with_rep(df$exonsize))
            } else {
                exon_separation_bars <- cumsum(rev(unique_with_rep(
                                                                df$exonsize)))
            }
        
            expected_cols <- c("nuc", "value", "qinf", "qsup", "group")
            df<-df[,which(colnames(df) %in% expected_cols)]
            expected_class <- c("integer", "numeric", "factor",
                                                    rep("numeric", 2))
            stopifnot(all(expected_cols %in% colnames(df)))
            print(vapply(df, class, character(1)))
            stopifnot(all(vapply(df, class, character(1)) == expected_class))
            
            ggplot(df, aes(x=nuc, y=value, ymin=qinf, ymax=qsup)) +
                geom_ribbon(aes(fill = group), alpha=0.3) +
                geom_line(aes(color = group), size=1) +
                geom_vline(xintercept = exon_separation_bars, 
                                                    linetype = "dotted") +
                theme(panel.grid.major = element_line()) +
                theme(panel.grid.minor = element_line()) +
                theme(panel.background = element_blank()) +
                theme(panel.background = element_rect()) +
                theme_bw(base_size = 20)
                
        } else { #if multiple regions/genes in the data_frame
            print('Multiple regions/genes')
            
            is_genes_unflipped <- unlist(lapply(map(unique(df$region), 
                                ~ df[which(df$region == .x),]$nuc), ascending))
            if (all(is_genes_unflipped)){
                exon_separation_bars <- cumsum(unique_with_rep(df$exonsize))
            } else {
                exon_length_by_exon_by_gene <- map(unique(df$region),
                        ~unique_with_rep(df$exonsize[which(df$region == .x)]))
                exon_separation_bars <- cumsum(unlist(rev_if_flipped(
                                                    exon_length_by_exon_by_gene, 
                                                    is_genes_unflipped)))
            }
            
            expected_cols <- c("nuctot", "value", "qinf", "qsup", "group")
            df<-df[,which(colnames(df) %in% expected_cols)]
            expected_class <- c("integer", "numeric", "factor",
                                                    rep("numeric", 2))
            stopifnot(all(expected_cols %in% colnames(df)))
            print(vapply(df, class, character(1)))
            stopifnot(all(vapply(df, class, character(1)) == expected_class))
        
            #adjustment of nuctot in case of subset of original data frame
            df$nuctot <- df$nuctot - min(df$nuctot) +1
            
            ggplot(df, aes(x=nuctot, y=value, ymin=qinf, ymax=qsup)) +
                geom_ribbon(aes(fill = group), alpha = 0.3) +
                geom_line(aes(color = group), size = 1) +
                geom_vline(xintercept = exon_separation_bars, 
                                                    linetype = "dotted") +
                theme(panel.grid.major = element_line()) +
                theme(panel.grid.minor = element_line()) +
                theme(panel.background = element_blank()) +
                theme(panel.background = element_rect()) +
                theme_bw(base_size = 20)
        }
    } else if (('bin' %in% colnames(df)) & ('nuc' %in% colnames(df))) { 
        print('rnaseq with bins')
        
        expected_cols <- c("bin", "value", "qinf", "qsup", "design")
        df<-df[,which(colnames(df) %in% expected_cols)]
        print(vapply(df, class, character(1)))
        expected_class <- c("factor", "integer", rep("numeric", 3))
        stopifnot(all(expected_cols %in% colnames(df)))
        stopifnot(all(vapply(df, class, character(1)) == expected_class))

        ggplot(df, aes(x=bin, y=value, ymin=qinf, ymax=qsup)) +
            geom_ribbon(aes(fill=design), alpha=0.3) +
            geom_line(aes(color=design), size=1) +
            theme(panel.grid.major = element_line()) +
            theme(panel.grid.minor = element_line()) +
            theme(panel.background = element_blank()) +
            theme(panel.background = element_rect()) +
            theme_bw(base_size = 20)
    }
}