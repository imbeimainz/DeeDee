#' deedee_upset
#'
#' Upset plot on the sets of DE features for the different analyses included
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param mode Character value, specifies which subset of features to include
#' in the overlap computations. Can be either of the following: "both", "up",
#' "down", or "both_colored".
#' @param min_setsize Numeric value, specifying the minimal number of
#' observations in an intersection for it to be included in the upset plot.
#' @param pthresh Numeric value, corresponding to the p-value to use as a
#' threshold to subset the features to include.
#'
#' @return A plot object, drawn with the `ComplexUpset` package.
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
#' deedee_upset(dde, pthresh = 0.05, mode = "both_colored", min_setsize = 10)
#'
deedee_upset <- function(dde,
                         mode = "both_colored",
                         min_setsize = 10,
                         pthresh = 0.05) {

  if (length(dea(dde)) < 2)
    stop("Please provide a dde object including at least two DE analyses")

  if (!is.numeric(min_setsize) | !(min_setsize > 0))
    stop("Please specify a positive integer for the min_setsize parameter")

  if (!is.numeric(pthresh) | !(pthresh > 0 & pthresh <= 1))
    stop("Please specify a valid p-value threshold")

  mode <- match.arg(mode, c("up", "down", "both", "both_colored"))


  dea_list <- get_dea_list(dde)
  for (i in seq_len(length(dea_list))) {
    dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$padj < pthresh) # pthresh

    if (length(dea_list[[i]][[1]]) == 0) {
      return(NULL)
    }

    dea_list[[i]] <- dea_list[[i]]["log2FoldChange"] # removing p-value column
    dea_list[[i]]$log2FoldChange[dea_list[[i]]$log2FoldChange < 0] <- -1 # 1/0/-1 distinction
    dea_list[[i]]$log2FoldChange[dea_list[[i]]$log2FoldChange > 0] <- 1

    if (mode == "both" || mode == "both_colored") {
      dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$log2FoldChange != 0)
    }

    if (mode == "up") {
      dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$log2FoldChange == 1)
    }

    if (mode == "down") {
      dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$log2FoldChange == -1)
    }
    dea_list[[i]] <- tibble::rownames_to_column(dea_list[[i]])
    names(dea_list[[i]]) <- c("rowname", names(dea_list[i]))
  }

  comp <- dplyr::full_join(dea_list[[1]], dea_list[[2]],
                           by = "rowname",
                           copy = FALSE
  )

  if (length(dea_list) > 2) {
    for (i in 3:length(dea_list)) {
      comp <- dplyr::full_join(comp, dea_list[[i]], by = "rowname", copy = FALSE)
    }
  }
  comp <- tibble::column_to_rownames(comp)
  contrasts <- colnames(comp)

  if (mode == "both_colored") {
    count <- vector(mode = "numeric", length = length(comp[[1]]))
    comp <- cbind(comp, count)
    comp["count"] <- rowSums(comp, na.rm = TRUE)

    na_count <- vector(mode = "numeric", length = length(comp[[1]]))
    comp <- cbind(comp, na_count)
    comp["na_count"] <- rowSums(is.na(comp))

    dr <- vector(mode = "character", length = length(comp[[1]]))
    comp <- cbind(comp, dr)

    nsets <- length(dea_list)

    for (i in seq_len(length(comp[[1]]))) {
      if ((comp[i, "count"] + comp[i, "na_count"]) == nsets) {
        comp[i, "dr"] <- "up"
      } else if ((comp[i, "count"] - comp[i, "na_count"]) == -nsets) {
        comp[i, "dr"] <- "down"
      } else {
        comp[i, "dr"] <- "different"
      }
    }

    comp <- subset(comp, select = -c(count, na_count))
    comp[is.na(comp)] <- FALSE
    comp[comp == -1 | comp == 1] <- TRUE
  }

  if (mode == "both_colored") {
    col_up <- viridis::viridis(n = 1, begin = 0.4, option = "magma")
    col_down <- viridis::viridis(n = 1, begin = 0.9, option = "magma")
    res <- ComplexUpset::upset(comp,
                               contrasts,
                               guide = NULL,
                               base_annotations = list(
                                 "Intersection of genes" = ComplexUpset::intersection_size(
                                   counts = FALSE,
                                   mapping = ggplot2::aes(fill = dr)
                                 )
                                 + ggplot2::scale_fill_manual(
                                   values = c(
                                     "different" = "dark grey",
                                     "down" = col_down,
                                     "up" = col_up
                                   ),
                                   name = "Regulation direction"
                                 )
                               ),
                               width_ratio = 0.1,
                               min_size = min_setsize
    )
  } else {
    comp <- (comp != 0 & !is.na(comp))
    comp <- as.data.frame(comp)
    res <- ComplexUpset::upset(comp,
                               contrasts,
                               guide = NULL,
                               width_ratio = 0.1,
                               min_size = min_setsize
    )
  }

  return(res)
}


# deedee_upset(dde,
#              pthresh = 0.05,
#              mode = "both_colored",
#              min_setsize = 10)
