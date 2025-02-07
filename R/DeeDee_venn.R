#' deedee_venn
#'
#' A Venn diagram on the DE features from different analyses
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param mode Character value, one of "both", "up", or "down". Specifies which
#' set of features to focus on for the overlap calculations.
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
#' library("dplyr") ## to have the inner_join, needed by ggvenn
#' deedee_venn(dde, pthresh = 0.05, mode = "both")
#'
deedee_venn <- function(dde,
                        mode = "both",
                        pthresh = 0.05) {

  if (length(dea(dde)) < 2)
    stop("Please provide a dde object including at least two DE analyses")

  if (!is.numeric(pthresh) | !(pthresh > 0 & pthresh <= 1))
    stop("Please specify a valid p-value threshold")

  mode <- match.arg(mode, c("up", "down", "both"))


  dea_list <- get_dea_list(dde)
  for (i in seq_len(length(dea_list))) {
    dea_list[[i]] <- subset(
      dea_list[[i]],
      dea_list[[i]]$padj < pthresh
    )

    if (length(dea_list[[i]][[1]]) == 0) {
      return(NULL)
    }

    if (mode == "up") {
      dea_list[[i]] <- subset(
        dea_list[[i]],
        dea_list[[i]]$log2FoldChange > 0
      )
    }
    if (mode == "down") {
      dea_list[[i]] <- subset(
        dea_list[[i]],
        dea_list[[i]]$log2FoldChange < 0
      )
    }
    dea_list[[i]] <- dea_list[[i]]["log2FoldChange"] # removing p-value column
    colnames(dea_list[[i]]) <- c(paste("log2FoldChange", i, sep = ""))
    dea_list[[i]] <- as.matrix(dea_list[[i]]) # conversion to matrix
    names <- rownames(dea_list[[i]])
    dea_list[[i]] <- as.vector(dea_list[[i]]) # conversion to vector
    names(dea_list[[i]]) <- names
    dea_list[[i]] <- sort(dea_list[[i]], decreasing = TRUE)
    dea_list[[i]] <- names(dea_list[[i]])
  }

  pal <- c(viridis::viridis(length(dea_list), option = "magma"))

  res <- ggvenn::ggvenn(dea_list,
                        fill_alpha = 0.2,
                        fill_color = pal,
                        show_percentage = FALSE,
                        stroke_color = "grey80",
                        set_name_size = 4
  )

  return(res)
}

# deedee_venn(dde, pthresh = 0.05, mode = "both")











