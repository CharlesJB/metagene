#' Produce a metagene plot
#' 
#' @param df a \code{data.frame} obtained with the \code{get_data_frame}
#' function. Must have the following columns: "region", "design", "bin", "value",
#' "qinf" and "qsup".
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
	df$group <- paste(df$region,df$design,sep="_")
	df$group <- as.factor(df$group)
	df<-df[,3:7]
	expected_cols <- c("group", "bin", "value", "qinf", "qsup")
	expected_class <- c("integer", rep("numeric", 3), "factor")
	stopifnot(all(expected_cols %in% colnames(df)))
	stopifnot(all(vapply(df[df,on=expected_cols], class, character(1)) == expected_class))
	#why not : stopifnot(all(vapply(df, class, character(1)) == expected_class)) because df[df,on=expected_cols] is = to df ?

    ggplot(df, aes(x=bin, y=value, ymin=qinf, ymax=qsup)) +
        geom_ribbon(aes(fill=group), alpha=0.3) +
        geom_line(aes(color=group), size=1) +
        theme(panel.grid.major = element_line()) +
        theme(panel.grid.minor = element_line()) +
        theme(panel.background = element_blank()) +
        theme(panel.background = element_rect()) +
        theme_bw(base_size = 20)
}
