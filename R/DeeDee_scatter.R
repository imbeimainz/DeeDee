#' deedee_scatter
#'
#' Scatter plot for the log fold change in the different DE analyses
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param select1 Numeric value, corresponding to the order of the element in the
#' `dde` object to be selected first.
#' @param select2 Numeric value, corresponding to the order of the element in the
#' `dde` object to be selected as second.
#' @param color_by Character value, either "pval1" or "pval2" to color the
#' individual points mapping them to the p-value of either DE result set.
#' @param pthresh Numeric value, corresponding to the p-value to use as a
#' threshold to subset the features to include.
#'
#' @return A `ggplot` plot object.
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
#' deedee_scatter(dde, 1, 3, color_by = "pval2")
#' deedee_scatter(dde, 1, 2)
#' deedee_scatter(dde, 1, 4)
#' deedee_scatter(dde, 3, 2)
#'
deedee_scatter <- function(dde,
                           select1 = 1,
                           select2 = 2,
                           color_by = "pval1",
                           pthresh = 0.05) {

  if (length(dea(dde)) < 2)
    stop("Please provide a dde object including at least two DE analyses")

  if (!is.numeric(select1) | !(select1 > 0 & select1 <= length(dea(dde))))
    stop("Please specify a valid entry to use as 1st selection")
  if (!is.numeric(select2) | !(select2 > 0 & select2 <= length(dea(dde))))
    stop("Please specify a valid entry to use as 2nd selection")

  if (!is.numeric(pthresh) | !(pthresh > 0 & pthresh <= 1))
    stop("Please specify a valid p-value threshold")

  color_by <- match.arg(color_by, c("pval1", "pval2"))


  dea_names <- names(dea(dde))
  names1 <- dea_names[select1]
  names2 <- dea_names[select2]

  de_1 <- get_dea_df(dde, names1)
  de_2 <- get_dea_df(dde, names2)
  data_red <- list(
    tibble::rownames_to_column(as.data.frame(de_1)),
    tibble::rownames_to_column(as.data.frame(de_2))) # selected samples

  # 928
  # 2692
  # 185

  comp <- dplyr::inner_join(data_red[[1]],
                            data_red[[2]],
                            by = "rowname",
                            copy = FALSE
  )

  # for (i in 1:length(data_red)) {
  #   data_red[i][[1]] <- subset(
  #     data_red[i][[1]],
  #     data_red[i][[1]]$pval < pthresh
  #   )
  #   data_red[i][[1]] <- tibble::rownames_to_column(data_red[i][[1]])
  # }
  # comp <- dplyr::inner_join(data_red[1][[1]],
  #                           data_red[2][[1]],
  #                           by = "rowname",
  #                           copy = FALSE
  # )
  # comp <- comp[stats::complete.cases(comp[colnames(comp)]), ]

  names(comp) <- c("gene_id",
                   "logFC1", "punadj1", "pval1",
                   "logFC2", "punadj2", "pval2")


  comp_sel <- comp[(comp$pval1 <= pthresh) & (comp$pval2 <= pthresh), ]

  comp_sel$info <- paste0(
    comp_sel$gene_id,": ", comp_sel$logFC1, ", ", comp_sel$logFC2, "\n something else\nanother line"
  )

  if (length(comp[, 1]) == 0) {
    return(NULL)
  }

  res <- ggplot2::ggplot(
    data = comp_sel,
    aes(.data[["logFC1"]], .data[["logFC2"]],
        col = .data[[color_by]],
        label = .data[["gene_id"]],
        info = .data[["info"]]
    )) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis(option = "magma") +
    ggplot2::xlab(names1) +
    ggplot2::ylab(names2) +
    ggplot2::labs(color = color_by) +
    ggplot2::theme_light()

  return(res)
}
