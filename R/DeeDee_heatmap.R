#' deedee_heatmap
#'
#' Construct a logFC heatmap on a subset of features
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param show_first Numeric value, specifying the number of features to include.
#' The ranking is done according to ... TODO
#' @param show_gene_names Logical value, whether to display the gene names on the
#' heatmap's side.
#' @param dist Character value, specifying the distance type to use in the call
#' to `ComplexHeatmap`. Accordingly, it can be a pre-defined character
#' ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",
#' "pearson", "spearman", "kendall").
#' @param clust Character value. Defines the method to perform hierarchical
#' clustering, passed to `hclust`, as used in `ComplexHeatmap`.
#' @param show_na Logical value, whether to include features that have `NA` value
#' for the log fold change.
#' @param pthresh Numeric value, corresponding to the p-value to use as a
#' threshold to subset the features to include.
#'
#' @return  A plot object generated via the `ComplexHeatmap` package.
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
#' deedee_heatmap(dde,
#'                pthresh = 0.05, show_first = 25,
#'                show_gene_names = FALSE, dist = "euclidean",
#'                clust = "average", show_na = FALSE
#' )
#'
deedee_heatmap <- function(dde,
                           show_first = 25,
                           show_gene_names = FALSE,
                           dist = "euclidean",
                           clust = "average",
                           show_na = FALSE,
                           pthresh = 0.05) {
  if (length(dea(dde)) < 2)
    stop("Please provide a dde object including at least two DE analyses")

  if (!is.numeric(show_first) | !(show_first > 0))
    stop("Please specify a positive integer for the show_first parameter")

  if (!is.numeric(pthresh) | !(pthresh > 0 & pthresh <= 1))
    stop("Please specify a valid p-value threshold")

  stopifnot(is.logical(show_gene_names))
  stopifnot(is.logical(show_na))

  stopifnot(is.character(dist))
  stopifnot(is.character(clust))


  dea_list <- get_dea_list(dde)
  # and then proceed as for the old implementation

  for (i in seq_len(length(dea_list))) {
    dea_list[[i]] <- subset(dea_list[[i]], dea_list[[i]]$padj < pthresh) # pthresh
    dea_list[[i]] <- dea_list[[i]]["log2FoldChange"] # removing p-value column
    colnames(dea_list[[i]]) <- names(dea_list)[[i]] # creating unique colnames
    dea_list[[i]] <- tibble::rownames_to_column(dea_list[[i]])
  }

  comp <- dplyr::full_join(dea_list[[1]],
                           dea_list[[2]],
                           by = "rowname",
                           copy = FALSE
  )
  if (length(dea_list) > 2) {
    for (i in 3:length(dea_list)) {
      comp <- dplyr::full_join(comp, dea_list[[i]], by = "rowname", copy = FALSE)
    }
  }


  row.names(comp) <- comp$rowname

  # comp <- subset(comp, select = -c(rowname)) # removing column with rownames
  comp <- comp[, !(names(comp) %in% "rowname")]

  comp <- comp[rowSums(!is.na(comp)) >= floor(length(comp) / 2) + 1, ]
  if (show_na == FALSE) {
    comp <- comp[stats::complete.cases(comp[colnames(comp)]), ]
  }
  comp <- as.matrix(comp)

  if (length(comp[, 1]) == 0) {
    return(NULL)
  }

  if (show_gene_names == FALSE) {
    rownames(comp) <- c()
  }

  col <- viridis::viridis(n = 15, option = "magma")

  res <- ComplexHeatmap::Heatmap(
    comp[seq_len(min(show_first, length(comp[, 1]))), ],
  name = "logFC",
  col = col,
  clustering_distance_rows = dist,
  clustering_method_rows = clust
  )

  # ComplexHeatmap::draw(res)

  return(res)

}



#' deedee_heatmap(dd_list_original,
#'                pthresh = 0.05, show_first = 25,
#'                show_gene_names = FALSE, dist = "euclidean",
#'                clust = "average", show_na = FALSE
#' )

# deedee_heatmap(dde = dde,
#                pthresh = 0.05, show_first = 25,
#                show_gene_names = FALSE, dist = "euclidean",
#                clust = "average", show_na = FALSE
# )
#
# deedee_heatmap(dde = dde,
#                pthresh = 0.05, show_first = 25,
#                show_gene_names = TRUE, dist = "euclidean",
#                clust = "average", show_na = TRUE
# )
# TODOs: better scale, blue to red
# TODOs: do this by overriding the selection of genes
# TODOs: add option to return the data frame instead of the plotted object

#' Title
#'
#' @param dde TODO
#' @param assay_name TODO
#' @param pvalue_threshold TODO
#' @param logfc_threshold TODO
#' @param custom_genelist TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
#' NULL ## TODO: this needs to be fully in
deedee_deheat <- function(dde,
                          assay_name = "vst",
                          pvalue_threshold = 0.05,
                          logfc_threshold = 5,
                          custom_genelist = NULL) {

  # checks: need an assay to plot the expression values

  # extract DE genes in any contrast
  dea_list <- get_dea_list(dde)

  dea_features <- lapply(names(dea_list), function(dea) {
    this_dea <- dea_list[[dea]]

    this_de_features <-
      rownames(this_dea[(abs(this_dea$log2FoldChange) >= logfc_threshold) &
                        (this_dea$padj <= pvalue_threshold), ])
  })
  names(dea_features) <- names(dea_list)

  # merge them
  de_features <- unique(unlist(dea_features))


  # also using the custom genelist
  if (!is.null(custom_genelist)) {
    custom_genelist_present <- intersect(rownames(dde), custom_genelist)
    de_features <- unique(c(de_features, custom_genelist_present))
  }

  # prepare annotation columns
  logfc_df <- data.frame(matrix(nrow = length(de_features),
                                ncol = length(dea_list)))
  rownames(logfc_df) <- de_features
  colnames(logfc_df) <- names(dea_list)

  # matching the logFC info into it
  for (i in names(dea_list)) {
    this_dea <- dea_list[[i]]
    logfc_vector <- this_dea[rownames(logfc_df), ]$log2FoldChange
    logfc_df[, i] <- logfc_vector
  }

  # annotation for the samples as well?

  # take out a sensible assay - and vst that?
  assay_to_use <- assays(dde)[[assay_name]]
  assay_to_use <- vst(assay_to_use)
  assay_to_use_de <- assay_to_use[de_features , ]

  # TODO if centered?
  assay_to_use_de <- assay_to_use_de - rowMeans(assay_to_use_de)

  # create the heatmap itself
  # pheatmap(assay_to_use_de, scale = "row", annotation_row = logfc_df)

  max_logfc <- max(abs(logfc_df))
  col_fun <- circlize::colorRamp2(c(-max_logfc, 0, max_logfc), c("blue", "white", "red"))

  col_fun_list <- vector("list", length = length(dea_list))
  names(col_fun_list) <- names(dea_list)
  for (i in names(dea_list)) col_fun_list[[i]] <- col_fun

  ha <- Heatmap(assay_to_use_de,
                name = "Z-score\nexpression\nvalues",
                left_annotation = rowAnnotation(df = logfc_df, col = col_fun_list),
                show_row_names = FALSE)

  return(ha)
}
