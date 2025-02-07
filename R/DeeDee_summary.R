#' deedee_summary
#'
#' Generates a comprehensive summary on the DE analyses, in one command
#'
#' Most parameters in this function call are related to a specific function of the
#' `DeeDee` package. Please refer to the individual help pages for more details,
#' reminding that the parameter naming scheme is mirroring the pattern
#' `function_parameter`, e.g. `scatter_color_by` controlling the `color_by`
#' value in the [deedee_scatter()] function.
#'
#' @param dde A [DeeDeeExperiment] object.
#' @param template Character value. Defines the location of the Rmd template to use
#' for generating the summary
#' @param output_path Character value, name of the file the report will be saved
#' to. Defaults to `DeeDeeSummary.html`, in the current working directory.
#' @param overwrite Logical value, whether to overwrite an existing file if
#' choosing the exact same name.
#' @param pthresh Numeric value, corresponding to the p-value to use as a
#' threshold to subset the features to include. This parameter is shared across
#' a number of individual `DeeDee` functions.
#' @param scatter_select1 From [deedee_scatter()] - Numeric value, corresponding
#' to the order of the element in the dde object to be selected first.
#' @param scatter_select2 From [deedee_scatter()] - Numeric value, corresponding
#' to the order of the element in the dde object to be selected as second.
#' @param scatter_color_by From [deedee_scatter()] - Character value, either
#' "pval1" or "pval2" to color the individual points mapping them to the p-value
#' of either DE result set.
#' @param heatmap_show_first From [deedee_heatmap()] - Numeric value, specifying
#' the number of features to include.
#' @param heatmap_show_gene_names From [deedee_heatmap()] - Logical value,
#' whether to display the gene names on the heatmap's side.
#' @param heatmap_dist From [deedee_heatmap()] - Character value, specifying the
#' distance type to use in the call to ComplexHeatmap.
#' @param heatmap_clust From [deedee_heatmap()] - Character value. Defines the
#' method to perform hierarchical clustering, passed to `hclust`, as used in
#' `ComplexHeatmap`.
#' @param heatmap_show_na From [deedee_heatmap()] - Logical value, whether to
#' include features that have `NA` value for the log fold change.
#' @param venn_mode From [deedee_venn()] - Character value, one of "both", "up",
#' or "down". Specifies which set of features to focus on for the overlap
#' calculations.
#' @param upset_mode From [deedee_upset()] - Character value, specifies which
#' subset of features to include in the overlap computations.
#' Can be either of the following: "both", "up", "down", or "both_colored".
#' @param upset_min_setsize From [deedee_upset()] - Numeric value, specifying
#' the minimal number of observations in an intersection for it to be included
#' in the upset plot.
#' @param qqmult_ref From [deedee_qqmult()] - A numeric value, corresponding to
#' the order of the element in the `dde` object to be used as a reference for the
#' multi Q-Q plot.
#' @param cat_ref From [deedee_cat()] - A numeric value, corresponding to the
#' order of the element in the dde object to be used as a reference in the CAT
#' plot.
#' @param cat_maxrank From [deedee_cat()] - A numeric value. Indicates the
#' maximum ranked feature to include when computing the Concordance At the Top.
#' @param cat_mode From [deedee_cat()] - A character value, could be one of "up",
#' "down", or "both". Defines which features to include in the computations.
#' @param render_quiet Logical value. Whether to render the report quietly when
#' calling `rmarkdown::render()`.
#' @param silent Logical value, prints a message once completed if set to TRUE.
#' @param open_file Logical value, if set to TRUE opens up the report once
#' generated into a browser window.
#'
#' @seealso
#' [deedee_cat()], [deedee_heatmap()], [deedee_qq()], [deedee_qqmult()],
#' [deedee_scatter()], [deedee_upset()], and [deedee_venn()]
#'
#' @return Generates a fully fledged report in the location specified by
#' `output_path`, and returns (invisibly) the name of the generated report.
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
#' # deedee_summary(dde)
deedee_summary <- function(dde,
                           template = system.file("extdata",
                                                  "summary_template.Rmd",
                                                  package = "DeeDee"),
                           output_path = "DeeDee_Summary.html",
                           overwrite = FALSE,
                           pthresh = 0.05,
                           scatter_select1 = 1,
                           scatter_select2 = 2,
                           scatter_color_by = "pval1",
                           heatmap_show_first = 25,
                           heatmap_show_gene_names = FALSE,
                           heatmap_dist = "euclidean",
                           heatmap_clust = "average",
                           heatmap_show_na = FALSE,
                           venn_mode = "both",
                           upset_mode = "both_colored",
                           upset_min_setsize = 10,
                           qqmult_ref = 1,
                           cat_ref = 1,
                           cat_maxrank = 1000,
                           cat_mode = "up",
                           render_quiet = TRUE,
                           silent = FALSE,
                           open_file = TRUE) {

  # checkmate::assert_list(deedee_list, type = "data.frame", min.len = 2)
  # for (i in 1:length(deedee_list)) {
  #   checkmate::assert_data_frame(deedee_list[[i]], type = "numeric")
  # }
  # checkmate::assert_number(pthresh, lower = 0, upper = 1)
  #
  # checkmate::assert_logical(overwrite)
  #
  # checkmate::assert_number(scatter_select1, lower = 1, upper = length(deedee_list))
  # checkmate::assert_number(scatter_select2, lower = 1, upper = length(deedee_list))
  # choices <- c("pval1", "pval2")
  # checkmate::assert_choice(scatter_color_by, choices)
  #
  # checkmate::assert_number(heatmap_show_first, lower = 1)
  # checkmate::assert_logical(heatmap_show_gene_names)
  # choices1 <- c("euclidean", "manhattan", "pearson", "spearman")
  # checkmate::assert_choice(heatmap_dist, choices1)
  # choices2 <- c("single", "complete", "average", "centroid")
  # checkmate::assert_choice(heatmap_clust, choices2)
  # checkmate::assert_logical(heatmap_show_na)
  #
  # choices <- c("up", "down", "both")
  # checkmate::assert_choice(venn_mode, choices)
  #
  # choices <- c("up", "down", "both", "both_colored")
  # checkmate::assert_choice(upset_mode, choices)
  # checkmate::assert_number(upset_min_setsize, lower = 0)
  #
  # checkmate::assert_number(qqmult_ref, lower = 1, upper = length(deedee_list))
  #
  # checkmate::assert_number(cat_ref, lower = 1, upper = length(deedee_list))
  # checkmate::assert_number(cat_maxrank, lower = 1)
  # choices <- c("up", "down", "both")
  # checkmate::assert_choice(cat_mode, choices)

  # template <- system.file("extdata",
  #                         "summary_template.Rmd",
  #                         package = "DeeDee"
  # )



  file.exists(template)

  if (file.exists(output_path)) {
    if (!overwrite) {
      stop("Your declared output file already exists. ",
           "Set overwrite = TRUE to proceed anyway.",
           call. = FALSE
      )
    }
  }

  output_rmd <- paste(unlist(strsplit(output_path,
                                      split = ".",
                                      fixed = TRUE
  ))[1], ".Rmd", sep = "")

  file.copy(from = template, to = output_rmd, overwrite = TRUE)

  args <- list()
  args$input <- output_rmd
  args$output_format <- "html_document"
  args$output_path <- output_path

  output_path <- rmarkdown::render(output_rmd,
                                   params = args,
                                   quiet = render_quiet
  )
  if (open_file == TRUE) {
    utils::browseURL(output_path)
  }

  file.remove(output_rmd)

  if (silent == FALSE) {
    message("Your summary has been generated in ", output_path)
  }

  invisible(output_path)
}

# deedee_summary(dde,
#                template = "inst/extdata/summary_template.Rmd", output_path = "v2.html")


