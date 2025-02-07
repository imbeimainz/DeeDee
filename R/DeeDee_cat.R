#' deedee_cat
#'
#' CAT plot - Concordance At the Top
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param ref A numeric value, corresponding to the order of the element in the
#' `dde` object to be used as a reference.
#' @param maxrank A numeric value. Indicates the maximum ranked feature to
#' include when computing the Concordance At the Top. Defaults to 1000.
#' @param mode A character value, could be one of "up", "down", or "both". Defines
#' which features to include in the computations. Defaults to "up".
#' @param pthresh Numeric value, corresponding to the p-value to use as a
#' threshold to subset the features to include. Defaults sensibly to 0.05.
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
#' deedee_cat(dde, ref = 1, maxrank = 1000, mode = "up", pthresh = 0.05)
#'
deedee_cat <- function(dde,
                       ref = 1,
                       maxrank = 1000,
                       mode = "up",
                       pthresh = 0.05) {

  if (length(dea(dde)) < 2)
    stop("Please provide a dde object including at least two DE analyses")

  if (!is.numeric(ref) | !(ref > 0 & ref <= length(dea(dde))))
    stop("Please specify a valid entry to use as reference")

  if (!is.numeric(maxrank) | !(maxrank > 1))
    stop("Please specify a positive integer (larger than 1) for the maximum rank")

  if (!is.numeric(pthresh) | !(pthresh > 0 & pthresh <= 1))
    stop("Please specify a valid p-value threshold")

  mode <- match.arg(mode, c("up", "down", "both"))


  dea_list <- get_dea_list(dde)

  for (i in seq_along(dea_list)) {
    dea_list[i][[1]] <- subset(
      dea_list[i][[1]],
      dea_list[i][[1]]$padj < pthresh
    )

    if (length(dea_list[i][[1]][[1]]) == 0) {
      return(NULL)
    }

    dea_list[i][[1]] <- dea_list[i][[1]]["log2FoldChange"] # remove p-value column
    colnames(dea_list[i][[1]]) <- c(paste("log2FoldChange", i, sep = ""))
    dea_list[i][[1]] <- as.matrix(dea_list[i][[1]]) # conversion to matrix
    names <- rownames(dea_list[i][[1]])
    dea_list[i][[1]] <- as.vector(dea_list[i][[1]]) # conversion to vector
    names(dea_list[i][[1]]) <- names
    if (mode == "up") {
      dea_list[i][[1]] <- sort(dea_list[i][[1]], decreasing = TRUE)
    } else if (mode == "down") {
      dea_list[i][[1]] <- sort(dea_list[i][[1]], decreasing = FALSE)
    } else if (mode == "both") {
      dea_list[i][[1]] <- abs(dea_list[i][[1]])
      dea_list[i][[1]] <- sort(dea_list[i][[1]], decreasing = TRUE)
    }

    dea_list[i][[1]] <- names(dea_list[i][[1]])
  }

  output <- list()
  nm <- c()

  for (i in seq_len(length(dea_list))) {
    if (i != ref) {
      output[[i]] <- data.frame(
        rank = seq_len(min(maxrank, length(dea_list[i][[1]]))),
        concordance = NA
      )

      for (j in seq_len(nrow(output[[i]]))) {
        intsec <- intersect(
          dea_list[ref][[1]][seq_len(j)],
          dea_list[i][[1]][seq_len(j)]
        )
        output[[i]][[j, "concordance"]] <- length(intsec) / j
      }
      nm[[i]] <- names(dea_list[i])
    }
  }

  names(output) <- nm

  out_aggr <- dplyr::bind_rows(output, .id = "contrast")

  res <- ggplot2::ggplot(
    out_aggr,
    ggplot2::aes(.data$rank, .data$concordance, colour = .data$contrast)
  ) +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    viridis::scale_color_viridis(
      option = "magma", discrete = TRUE,
      begin = 0, end = 0.9
    ) +
    ggplot2::annotate("text",
                      label = paste("reference: ",
                                    names(dea_list)[ref],
                                    sep = ""
                      ),
                      x = maxrank * 0.8,
                      y = max(out_aggr[["concordance"]]) * 1.1
    )

  return(res)
}
