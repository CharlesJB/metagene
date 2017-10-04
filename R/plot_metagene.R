#' Produce a metagene plot
#' 
#' @param df a \code{data.frame} obtained with the \code{get_data_frame}
#' function. Must have the following columns: "group", "position", "value",
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
	expected_cols <- c("group", "position", "value", "qinf", "qsup")
	expected_class <- c("factor", rep("numeric", 4))
	stopifnot(all(expected_cols %in% colnames(df)))
	stopifnot(all(vapply(df[expected_cols], class, character(1)) == expected_class))

    ggplot(df, aes(x=position, y=value, ymin=qinf, ymax=qsup)) +
        geom_ribbon(aes(fill=group), alpha=0.3) +
        geom_line(aes(color=group), size=1) +
        theme(panel.grid.major = element_line()) +
        theme(panel.grid.minor = element_line()) +
        theme(panel.background = element_blank()) +
        theme(panel.background = element_rect()) +
        theme_bw(base_size = 20)
}
