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
#' @param ... Additional arguments related to the `DESeq2::results()` function
#' that extracts DE results from the `dds` object.
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
                             de_results = NULL,
                             ...) {

  # old <- S4Vectors:::disableValidity()
  # if (!isTRUE(old)) {
  #   S4Vectors:::disableValidity(TRUE)
  #   on.exit(S4Vectors:::disableValidity(old))
  # }
  
  # set up de results list
  extracted_de_results <- list()
  
  # set up functional enrichment results list
  extracted_enrich_results <- list()
  
  if (!is.null(de_results)) {
    # capture variable name as a character
    entry_name <- deparse(substitute(de_results))
    de_results <- .check_de_results(de_results, entry_name)
    # assigned to the final de results list that might be later populated
    extracted_de_results <- append(extracted_de_results, de_results)
    # check validity of de results format before proceeding just in case
    extracted_de_results <- .check_de_results(extracted_de_results)
  }
  
  if (!is.null(se)) {
    # check if se is a dds object to set up whether to extract results from it or not
    if (is(se, "DESeqDataSet")) {
      message("Checking available contrasts in DESeqDataSet...")
      # available contrasts
      contrast_names <- DESeq2::resultsNames(se)
      # remove intercept to avoid conflicts, as it is not informative for comparisons
      contrast_names <- contrast_names[contrast_names != "Intercept"]
      
      if (length(contrast_names) == 0) {
        stop("No results found in the DESeqDataSet. You must run `DESeq()` first.")
      } else {
        for (contrast in contrast_names) {
          if (contrast %in% names(extracted_de_results)) {
            # if the name of the contrast in dds exists already in extracted_de_results
            # (i.e not empty) then don't extract results to avoid extracting the same twice
            message(
              "Contrast '",
              contrast,
              "' is already in extracted_de_results. Skipping..."
            )
            next # to skip this contrast
          }
          # if extracted_de_results is empty or the name of the contrast in dds doesn't
          # exist in extracted_de_results, then extract results
          # capturing additional arguments for the extraction if provided
          ### gotta figure out how to handle all the msg printed on the console!
          add_res_arg <- c(...)
          if (length(add_res_arg) == 0) {
            message("Extracting results for '",
                    contrast,
                    "' with default parameters")
          } else {
            message(
              "Extracting results for '",
              contrast,
              "' with the following parameters:",
              paste(paste0(
                names(add_res_arg), " = ", add_res_arg
              ), collapse = " ,")
            )
          }
          extracted_de_results[[contrast]] <- DESeq2::results(se, name = contrast, ...)
          
          # perform shrinkage if requested
          if (shrink_lfc == TRUE) {
            message(
              "Performing LFC shrinkage using 'apeglm' algorithm for contrast: '",
              contrast,
              " ' ..."
            )
            extracted_de_results[[contrast]] <- DESeq2::lfcShrink(se,
                                                          coef = contrast,
                                                          res = extracted_de_results[[contrast]],
                                                          type = "apeglm") # maybe give flexibility for this???
          }
        }
      }
    }
    
    else {
      if (!is(se, "RangedSummarizedExperiment")) {
        # check if it is SE and convert it into a RangedSE
        if (is(se, "SummarizedExperiment")) {
          se <- as(se, "RangedSummarizedExperiment")
        } else {
          stop("'se' must be a RangedSummarizedExperiment object")
        }
      }
      
    }
    
  }
  else {
    # if nothing is passed, return error
    if (length(extracted_de_results) == 0) {
      stop("You have to provide at least an se object or a de_results object!")
    }
    # if no se passed but extracted_de_results is not empty, create a mock from it
    message("creating a mock SE from the rows of the DE result objects")
    # mock up the se from the extracted_de_results
    first_de <- extracted_de_results[[1]]
    
    # independently of the class, the feature names are in the
    # rownames slot, TODO: check
    ids <- rownames(first_de)
    
    rd_mock <- DataFrame(gene_id = ids, row.names = ids)
    
    # way1
    se_mock <- SummarizedExperiment(assays = SimpleList(), rowData = rd_mock)
    # se_mock@NAMES <- NULL
    # rownames(se_mock) <- ids
    
    # no clue why this is strictly needed, but still it seems it is, if mocking up
    se <- as(se_mock, "RangedSummarizedExperiment")
  }
  
  
  # TODO: if no SE is really provided, instantiate some rownames, at least directly
  # from the rownames of the result objects
  # TODO: the row names are taken from the FIRST object in the de results then - or
  # from the union of all of them?
  
  if (is.null(de_results) && is.null(dds)) {
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
  names(extracted_de_results)
  # if not there, "force add"
  # TODO

  dde_ids <- rownames(se_out)

  dea_contrasts <- list()

  for (i in names(extracted_de_results)){
    this_de <- extracted_de_results[[i]]

    # do different things according to what these objects are
    if(is(this_de, "DESeqResults")) {
      input_deseq2 <- .importDE_DESeq2(se_out, this_de, i)
      se_out <- input_deseq2$se
      dea_contrasts[[i]] <- input_deseq2$dea_contrast

      
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

  #matched_ids <- match(rownames(res_de), rownames(se))

  matched_ids <- match(rownames(se), rownames(res_de)) # we align de res with se
  # only valid indices
  valid_matches <- !is.na(matched_ids)


  # Pre-fill rowData with NA
  rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name, "_pvalue")]]         <- NA
  rowData(se)[[paste0(de_name, "_padj")]]           <- NA


  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <- res_de$log2FoldChange[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches]         <- res_de$pvalue[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches]           <- res_de$padj[matched_ids[valid_matches]]


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


  # identify the logFC colss
  logFC_cols <- grep("^logFC", colnames(res_tbl), value = TRUE)
  
  
  matched_ids <- match(rownames(se), rownames(res_tbl)) # we align de res with se
  # only valid indices
  valid_matches <- !is.na(matched_ids)
  
  # pre-fill rowData with NA the assign the corresponding values only for matched
  # indices for logFC, accounting for the fact that the logFC column name in edgeR
  # depends on whether we have 1 or multple contrasts
  for (i in logFC_cols) {
    rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
    # assign correspionding values
    rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <- res_tbl$table[[i]][matched_ids[valid_matches]]
  }
  
  # pre-fill rowData with NA the assign the corresponding values for matched indices for pval and padj
  rowData(se)[[paste0(de_name, "_pvalue")]]         <- NA
  rowData(se)[[paste0(de_name, "_padj")]]           <- NA
  
  
  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches]         <- res_tbl$table$PValue[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches]           <- res_tbl$table$FDR[matched_ids[valid_matches]]
  
  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = res_tbl$comparison,
    metainfo_pvalue = NA,
    original_object = res_de,
    # object_name = deparse(substitute(res_tbl)),
    package = "edgeR"
  )
  
  return(list(se = se,
              dea_contrast = dea_contrast))
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

  # make sure there are at least 2 coefficients
  if (ncol(res_de$coefficients) < 2) { # we still need to manage the handling of 1 contrast
    warning(
      "The provided MArrayLM object has only ",
      ncol(res_de$coefficients),
      " coefficient(s). At least 2 are required."
    )
  }

  # extract columns
  res_tbl <- topTable(res_de,
                      coef = 2,
                      number = nrow(res_de),
                      sort.by = "none")


  #matched_ids <- match(rownames(res_tbl), rownames(se))

  matched_ids <- match(rownames(se), rownames(res_tbl)) # we align de res with se
  # only valid indices
  valid_matches <- !is.na(matched_ids)


  # Pre-fill rowData with NA
  rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name, "_pvalue")]]         <- NA
  rowData(se)[[paste0(de_name, "_padj")]]           <- NA


  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <- res_tbl$logFC[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches]         <- res_tbl$P.Value[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches]           <- res_tbl$adj.P.Val[matched_ids[valid_matches]]

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


# if de_results is one element, this function will capture the variable name and
# convert it into a character string to be passed to .check_de_results

.check_de_results <- function(x, entry_name=NULL) {
  # checks that:
  # it is a list
  # its component are either of the expected/accepted elements
  # checks their columns?

  #### update : checks and processes the input
  # if one single element, i.e not a list, convert if into a list of length 1
  if (is(x, "DGEExact") ||
      is(x, "DGELRT") || is(x, "MArrayLM") ||
      is(x, "DESeqResults")) {

    if(is.null(entry_name)) {
      stop("You must provide a name for your de results!")
    }

    x <- list(x)
    names(x) <- entry_name
  } else {
    # if a list
    ok_types <- unlist(lapply(x, function(arg) {
      is(arg, "DESeqResults") || is(arg, "DGEExact") ||
        is(arg, "DGELRT") || is(arg, "MArrayLM")
    }))

    if (!all(ok_types)) {
      stop("All elements in the list must be of type DESeqResults, DGEExact, DGELRT, or MArrayLM.")
    }
  }
  return(x)
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


