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
	if ('bin' %in% colnames(df)) { # if chipseq, for instance
	
		expected_cols <- c("bin", "value", "qinf", "qsup", "group")
		df<-df[,which(colnames(df) %in% expected_cols)]
		expected_class <- c("integer", rep("numeric", 3), "factor")
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
			
	} else if ('nuc' %in% colnames(df)) { # if rnaseq, for instance
		if (length(unique(df$region)) == 1) { #if only one region/gene in the data_frame
			
			#exon separation bars
			exon_separation_bars <- cumsum(unique(df[,which(colnames(df) %in% c('region','exonsize'))])[,2])
		
			expected_cols <- c("nuc", "value", "qinf", "qsup", "group")
			df<-df[,which(colnames(df) %in% expected_cols)]
			expected_class <- c("numeric", rep("numeric", 3), "factor")
			stopifnot(all(expected_cols %in% colnames(df)))
			stopifnot(all(vapply(df, class, character(1)) == expected_class))
			
			ggplot(df, aes(x=nuc, y=value, ymin=qinf, ymax=qsup)) +
				geom_ribbon(aes(fill = group), alpha=0.3) +
				geom_line(aes(color = group), size=1) +
				geom_vline(xintercept = exon_separation_bars, linetype = "dotted") +
				theme(panel.grid.major = element_line()) +
				theme(panel.grid.minor = element_line()) +
				theme(panel.background = element_blank()) +
				theme(panel.background = element_rect()) +
				theme_bw(base_size = 20)
				
		} else { #if multiple regions/genes in the data_frame
			
			#exon separation bars
			exon_separation_bars <- cumsum(unique(df[,which(colnames(df) %in% c('region','exonsize'))])[,2])
			
			expected_cols <- c("nuctot", "value", "qinf", "qsup", "group")
			df<-df[,which(colnames(df) %in% expected_cols)]
			expected_class <- c("numeric", rep("numeric", 3), "factor")
			stopifnot(all(expected_cols %in% colnames(df)))
			stopifnot(all(vapply(df, class, character(1)) == expected_class))
		
			#adjustment of nuctot in case of subset of original data frame
			df$nuctot = df$nuctot - min(df$nuctot) +1			
			
			ggplot(df, aes(x=nuctot, y=value, ymin=qinf, ymax=qsup)) +
				geom_ribbon(aes(fill = group), alpha = 0.3) +
				geom_line(aes(color = group), size = 1) +
				geom_vline(xintercept = exon_separation_bars, linetype = "dotted") +
				theme(panel.grid.major = element_line()) +
				theme(panel.grid.minor = element_line()) +
				theme(panel.background = element_blank()) +
				theme(panel.background = element_rect()) +
				theme_bw(base_size = 20)
		}
	}
}