#' Title
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param p_thresh Numeric value, corresponding to the threshold used to call a
#' feature as differentially expressed.
#' @param show_DEnumbers Logical value;
#'
#' @return A `ggplot` plot object, summarizing the amount of DE genes in each
#' comparison performed.
#'
#' @export
#'
#' @examples
#' data("de_named_list", package = "DeeDee")
#' library(SummarizedExperiment)
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
#' deedee_bars(dde)
deedee_bars <- function(dde,
                        p_thresh = 0.05,
                        show_DEnumbers = TRUE) {
  # checks and all
  # TODO


  # extract the summary infos
  dea_list <- get_dea_list(dde)

  dea_summary <- lapply(names(dea_list), function(dea) {
    this_dea <- dea_list[[dea]]
    dea_tally_up <- sum(this_dea$log2FoldChange >= 0 & this_dea$padj <= p_thresh)
    dea_tally_down <- sum(this_dea$log2FoldChange <= 0 & this_dea$padj <= p_thresh)

    dea_df <- data.frame(
      dea = dea,
      direction = c("Up", "Down"),
      value = c(dea_tally_up, dea_tally_down)
    )
  }) |> bind_rows()
  dea_summary$direction <- factor(dea_summary$direction, levels = c("Up", "Down"))

  res <- ggplot(dea_summary,
                aes(x = .data$dea,
                    y = .data$value,
                    fill = .data$direction)) +
    geom_bar(position = "dodge",
             stat="identity") +
    scale_fill_manual(values = c("#D62728FF", "#1F77B4FF")) +
    labs(x = "",
         y = "Number of DE genes",
         fill = "Direction of\nregulation",
         title = "DeeDeeExperiment summary",
         subtitle = paste0("Including ", length(dea(dde)), " DE analyses")
         ) +
    theme_bw()

  if (show_DEnumbers) {
    res <- res +
      geom_text(
        aes(label = .data$value),
        vjust = -0.7, position = position_dodge(width = 0.9)
      )
  }

  # alternatively: return the df itself

  return(res)
}

