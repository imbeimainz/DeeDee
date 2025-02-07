#' deedee_qqmult
#'
#' Q-Q plot for multiple comparisons of DE analyses
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param ref A numeric value, corresponding to the order of the element in the
#' `dde` object to be used as a reference.
#' @param pthresh Numeric value, corresponding to the p-value to use as a
#' threshold to subset the features to include.
#'
#' @return  A `ggplot` plot object.
#' @export
#'
#' @examples
#' data("de_named_list", package = "DeeDee")
#' library("SummarizedExperiment")
#'
#' rd_macrophage <- DataFrame(
#'   gene_id = rownames(de_named_list$ifng_vs_naive))
#' rownames(rd_macrophage) <- rownames(de_named_list$ifng_vs_naive)
#' se_macrophage_noassays <- SummarizedExperiment(
#'   assays = SimpleList(),
#'   rowData = rd_macrophage
#' )
#' dde <- DeeDeeExperiment(
#'   se_macrophage_noassays,
#'   de_results = de_named_list
#' )
#' dde
#'
#' deedee_qqmult(dde, pthresh = 0.05, ref = 1)
#'
deedee_qqmult <- function(dde,
                          ref = 1,
                          pthresh = 0.05) {
  if (length(dea(dde)) < 2)
    stop("Please provide a dde object including at least two DE analyses")

  if (!is.numeric(ref) | !(ref > 0 & ref <= length(dea(dde))))
    stop("Please specify a valid entry to use as reference")

  if (!is.numeric(pthresh) | !(pthresh > 0 & pthresh <= 1))
    stop("Please specify a valid p-value threshold")


  dea_list <- get_dea_list(dde)

  output <- list()
  nm <- c()

  for (i in seq_len(length(dea_list))) {
    if (i != ref) {
      output[[i]] <- data.frame(
        ggplot2::ggplot_build(deedee_qq(
          dde = dde,
          select1 = ref,
          select2 = i,
          as_line = TRUE
        ))$plot$data)

      nm[[i]] <- names(dea_list[i])
    }
  }

  names(output) <- nm

  res <- ggplot2::ggplot(
    dplyr::bind_rows(output, .id = "contrast"),
    ggplot2::aes(.data$x, .data$y, colour = .data$contrast)
  ) +
    ggplot2::xlab(names(dea_list)[ref]) +
    ggplot2::ylab("Contrasts") +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    viridis::scale_color_viridis(
      option = "magma", discrete = TRUE,
      begin = 0, end = 0.9
    )

  return(res)
}



