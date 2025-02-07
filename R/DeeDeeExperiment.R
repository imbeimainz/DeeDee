#' @name DeeDeeExperiment
#'
#' @title The DeeDeeExperiment class
#'
#' @aliases
#' DeeDeeExperiment
#' DeeDeeExperiment-class
#'
#' @description
#' The `DeeDeeExperiment` class is designed to represent
#' It inherits from the SummarizedExperiment class, and additionally stores
#' DE-related information via dedicated slots and `colData`.
#'
#' @param se A `SummarizedExperiment` object, that will be used as a scaffold to
#' store the DE related information.
#' @param de_results A named list of DE results, in any of the formats supported by
#' the `DeeDee` package (currently: results from DESeq2, edgeR, limma).
#'
#' @details
#' The `se` parameter can be optionally left unspecified. If this is the case,
#' the resulting `DeeDeeExperiment` object will contain as features the ones
#' specified by the provided components of the object supplied via the
#' `de_results` parameter.
#'
#' The conversion of the components of the `de_results` list will be handled via
#' conversion functions to uniform the names and set of information which will
#' be stored in the returned `DeeDeeExperiment` object.
#' The names of the list will be used to define the `contrasts` for the different
#' DE analyses included, which will determine the way to access the information
#' stored in the `dea` slot of the `DeeDeeExperiment` object
#'
#' Since a `DeeDeeExperiment` is also a `SummarizedExperiment` object, it can be
#' seamlessly provided downstream for visualization and in-depth exploration to
#' packages such as `iSEE` or similar.
#'
#' @return A `DeeDeeExperiment` object.
#' @export
#'
#' @author Lea Rothörl and Federico Marini
#'
#' @examples
#' data("de_named_list", package = "DeeDee")
#'
#' dde_onlyde <- DeeDeeExperiment(
#'   de_results = de_named_list
#' )
#'
#' # or, with a SE object as support - even without assay data available
#' library("SummarizedExperiment")
#'
#' rd_macrophage <- DataFrame(
#'   gene_id = rownames(de_named_list$ifng_vs_naive))
#' rownames(rd_macrophage) <- rownames(de_named_list$ifng_vs_naive)
#' se_macrophage_noassays <- SummarizedExperiment(
#'   assays = SimpleList(),
#'   rowData = rd_macrophage
#' )
#'
#' dde <- DeeDeeExperiment(
#'   se_macrophage_noassays,
#'   de_results = de_named_list
#' )
DeeDeeExperiment <- function(se = NULL,
                             de_results = NULL) {

  # old <- S4Vectors:::disableValidity()
  # if (!isTRUE(old)) {
  #   S4Vectors:::disableValidity(TRUE)
  #   on.exit(S4Vectors:::disableValidity(old))
  # }

  if (!is.null(de_results))
    .check_de_results(de_results)

  # se <- SummarizedExperiment(...)
  if(!is.null(se)) {
    if (!is(se, "RangedSummarizedExperiment")) {
      if (is(se, "SummarizedExperiment")) {
        se <- as(se, "RangedSummarizedExperiment")
      } else {
        stop("'se' must be a RangedSummarizedExperiment object")
      }
    }
  } else {
    if(is.null(de_results))
      stop("You have to provide at least an se object or a de_results object")

    message("creating a mock SE from the rows of the DE result objects")
    # mock up the se from the de_results
    first_de <- de_results[[1]]

    # independently of the class, the feature names are in the
    # rownames slot, TODO: check
    ids <- rownames(first_de)

    rd_mock <- DataFrame(
      gene_id = ids,
      row.names = ids)

    # way1
    se_mock <- SummarizedExperiment(
      assays = SimpleList(),
      rowData = rd_mock
    )
    # se_mock@NAMES <- NULL
    # rownames(se_mock) <- ids

    # no clue why this is strictly needed, but still it seems it is, if mocking up
    se <- as(se_mock, "RangedSummarizedExperiment")
  }

  # TODO: if no SE is really provided, instantiate some rownames, at least directly
  # from the rownames of the result objects
  # TODO: the row names are taken from the FIRST object in the de results then - or
  # from the union of all of them?

  if (is.null(de_results)) {
    object <- new("DeeDeeExperiment",
                  se,
                  dea = list()
    )

    # stash the package version
    metadata(object)[["version"]] <- packageVersion("DeeDee")

    return(object)
  }

  # TODO: does not have to relate to an SE which has all the slots and all
  # ...


  # TODO: additional checks
  se_out <- se

  # here is where I add the names in the rowData to make all info matched
  # checks on the names
  names(de_results)
  # if not there, "force add"
  # TODO

  dde_ids <- rownames(se_out)

  dea_contrasts <- list()

  for (i in names(de_results)){
    this_de <- de_results[[i]]

    # do different things according to what these objects are
    if(is(this_de, "DESeqResults")) {


      input_deseq2 <- .importDE_DESeq2(se_out, this_de, i)
      se_out <- input_deseq2$se
      dea_contrasts[[i]] <- input_deseq2$dea_contrast

      # matched_ids <- match(rownames(se_out), rownames(this_de))
      #
      # # if not tested, add NA - everywhere? -> pre-fill?
      # rowData(se_out)[[paste0(i,"_log2FoldChange")]] <- NA
      # rowData(se_out)[[paste0(i,"_pvalue")]] <- NA
      # rowData(se_out)[[paste0(i,"_padj")]] <- NA
      #
      # rowData(se_out)[[paste0(i,"_log2FoldChange")]][matched_ids] <- this_de$log2FoldChange
      # rowData(se_out)[[paste0(i,"_pvalue")]][matched_ids] <- this_de$pvalue
      # rowData(se_out)[[paste0(i,"_padj")]][matched_ids] <- this_de$padj
      #
      # dea_contrasts[[i]] <- list(
      #   alpha = metadata(this_de)$alpha,
      #   lfcThreshold = metadata(this_de)$lfcThreshold,
      #   metainfo_logFC = mcols(this_de)$description[colnames(this_de) == "log2FoldChange"],
      #   metainfo_pvalue = mcols(this_de)$description[colnames(this_de) == "pvalue"],
      #   original_object = this_de,
      #   package = "DESeq2"
      # ) is(res_de, "DGEExact") | is(res_de, "DGELRT")
    } else if (is(this_de, "DGEExact") | is(this_de, "DGELRT")) {
      input_edgeR <- .importDE_edgeR(se_out, this_de, i)
      se_out <- input_edgeR$se
      dea_contrasts[[i]] <- input_edgeR$dea_contrast
    } else if (is(this_de, "MArrayLM")) {
      input_limma <- .importDE_limma(se_out, this_de, i)
      se_out <- input_limma$se
      dea_contrasts[[i]] <- input_limma$dea_contrast
    }
  }

  # rowData(dde)[["new_rd"]] <- de_name

  object <- new("DeeDeeExperiment",
                se_out,
                dea = dea_contrasts
  )

  # stash the package version
  metadata(object)[["version"]] <- packageVersion("DeeDee")

  return(object)
}





# extends the rowData slot of the provided SE and returns also metadata

#' Import from `DESeq2` DE results
#'
#' @param se A `SummarizedExperiment` object
#' @param res_de A set of DE results, provided as `DESeqResults` as in the `DESeq2`
#' framework
#' @param de_name A character value, describing the contrast of interest. Will be
#' used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SummarizedExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
#'
#' @examples
#' # TODO example
.importDE_DESeq2 <- function(se, res_de, de_name) {
  # checks TODO:
  # correct object format
  stopifnot(
    is(res_de, "DESeqResults")
  )
  # contain the right columns
  stopifnot(
    all(c("log2FoldChange", "pvalue", "padj") %in% colnames(res_de))
  )

  # contain the feature ids

  # p value different from NA respect the 0-1 interval
  stopifnot(
    all(na.omit(res_de$pvalue <= 1)) & all(na.omit(res_de$pvalue > 0))
  )


  matched_ids <- match(rownames(res_de), rownames(se))

  # if not tested, add NA - everywhere? -> pre-fill?
  rowData(se)[[paste0(de_name,"_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name,"_pvalue")]] <- NA
  rowData(se)[[paste0(de_name,"_padj")]] <- NA

  rowData(se)[[paste0(de_name,"_log2FoldChange")]][matched_ids] <- res_de$log2FoldChange
  rowData(se)[[paste0(de_name,"_pvalue")]][matched_ids] <- res_de$pvalue
  rowData(se)[[paste0(de_name,"_padj")]][matched_ids] <- res_de$padj

  dea_contrast <- list(
    alpha = metadata(res_de)$alpha,
    lfcThreshold = metadata(res_de)$lfcThreshold,
    metainfo_logFC = mcols(res_de)$description[colnames(res_de) == "log2FoldChange"],
    metainfo_pvalue = mcols(res_de)$description[colnames(res_de) == "pvalue"],
    original_object = res_de,
    # object_name = deparse(substitute(res_de)),
    package = "DESeq2"
  )

  return(
    list(
      se = se,
      dea_contrast = dea_contrast
    )
  )
}


#' Import from edgeR DE results
#'
#' @param se A SummarizedExperiment object
#' @param res_de A set of DE results, provided by the `edgeR` framework (either a
#' `DGEExact` or a `DGELRT` object).
#' @param de_name A character value, describing the contrast of interest. Will be
#' used to compose the column names in the rowData slot.
#'
#' @return A list, containing the updated SummarizedExperiment object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' DeeDee framework.
#'
#' @noRd
#'
#' @examples
#' # TODO example
.importDE_edgeR <- function(se, res_de, de_name) {
  # checks object
  stopifnot(
    is(res_de, "DGEExact") | is(res_de, "DGELRT")
  )

  # extract columns
  res_tbl <- topTags(res_de, n = nrow(res_de), sort.by = "none")


  matched_ids <- match(rownames(res_tbl), rownames(se))

  # if not tested, add NA - everywhere? -> pre-fill?
  rowData(se)[[paste0(de_name,"_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name,"_pvalue")]] <- NA
  rowData(se)[[paste0(de_name,"_padj")]] <- NA

  rowData(se)[[paste0(de_name,"_log2FoldChange")]][matched_ids] <- res_tbl$table$logFC
  rowData(se)[[paste0(de_name,"_pvalue")]][matched_ids] <- res_tbl$table$PValue
  rowData(se)[[paste0(de_name,"_padj")]][matched_ids] <- res_tbl$table$FDR

  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = res_tbl$comparison,
    metainfo_pvalue = NA,
    original_object = res_de,
    # object_name = deparse(substitute(res_tbl)),
    package = "edgeR"
  )

  return(
    list(
      se = se,
      dea_contrast = dea_contrast
    )
  )
}



#' Import from `limma` DE results
#'
#' @param se A `SummarizedExperiment` object
#' @param res_de A set of DE results, provided in the `limma` framework (a `MArrayLM`
#' object).
#' @param de_name A character value, describing the contrast of interest. Will be
#' used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SummarizedExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
#'
#' @examples
#' # TODO example
#' # ... limma_de <- lmFit
#' # will provide the outout of lmFit - see its examples
.importDE_limma <- function(se, res_de, de_name) {
  # checks object
  stopifnot(is(res_de, "MArrayLM"))

  # extract columns
  res_tbl <- topTable(res_de,
                      coef = 2,
                      number = nrow(res_de),
                      sort.by = "none")


  matched_ids <- match(rownames(res_tbl), rownames(se))

  # if not tested, add NA - everywhere? -> pre-fill?
  rowData(se)[[paste0(de_name,"_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name,"_pvalue")]] <- NA
  rowData(se)[[paste0(de_name,"_padj")]] <- NA

  rowData(se)[[paste0(de_name,"_log2FoldChange")]][matched_ids] <- res_tbl$logFC
  rowData(se)[[paste0(de_name,"_pvalue")]][matched_ids] <- res_tbl$P.Value
  rowData(se)[[paste0(de_name,"_padj")]][matched_ids] <- res_tbl$adj.P.Val

  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = res_tbl$comparison,
    metainfo_pvalue = NA,
    original_object = res_de,
    # object_name = deparse(substitute(res_tbl)),
    package = "limma"
  )

  return(
    list(
      se = se,
      dea_contrast = dea_contrast
    )
  )


  # returns info (in the standardized manner)

}


.importDE_custom <- function(x) {

}



.check_de_results <- function(x) {
  # checks that:
  # it is a list
  # its component are either of the expected/accepted elements
  # checks their columns?

  stopifnot(is.list(x))
  stopifnot(length(x) > 0)

  ok_types <- unlist(lapply(x, function(arg) {
    is(arg, "DESeqResults") | is(arg, "DGEExact") | is(arg, "DGELRT") | is(arg, "MArrayLM")
  }))

  stopifnot(all(ok_types))

  return(TRUE)
}


deedee_import <- function(x) {

  # legacy code:

  # # ----------------------------- argument check ------------------------------
  # choices <- c("DESeq2", "edgeR", "limma")
  # checkmate::assertChoice(input_type, choices)
  #
  # if (input_type == "DESeq2") {
  #   checkmate::assertClass(data, "DESeqResults")
  #   logFC <- data$log2FoldChange
  #   pval <- data$padj
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- data@rownames
  # } else if (input_type == "edgeR") {
  #   checkmate::assertClass(data, "DGEExact")
  #   logFC <- data[["table"]][["logFC"]]
  #   pval <- data[["table"]][["PValue"]]
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- data[["genes"]][["genes"]]
  # } else if (input_type == "limma") {
  #   checkmate::assertDataFrame(data, types = "numeric")
  #   logFC <- data$logFC
  #   pval <- data$adj.P.Val
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- rownames(data)
  # }
  # return(input)
}


