#' @name DeeDeeExperiment-methods
#'
#' @title Methods for [DeeDeeExperiment] objects
#'
#' @aliases
#' dea
#' dea<-
#' add_dea
#' remove_dea
#' get_dea_df
#' get_dea_list
#'
#' @description
#' The [DeeDeeExperiment()] class provides a family of methods to get
#' and set DE-related information in [DeeDeeExperiment] objects.
#'
#' @param x A [DeeDeeExperiment()] object
#' @param value Replacement value for replacement methods.
#' @param dea A named list of DE results, in any of the formats supported by
#' the `DeeDee` package (currently: results from DESeq2, edgeR, limma).
#' @param dea_name Character value, specifying the name of the DE analysis to
#' get or remove
#'
#' @return Return value varies depending on the individual methods, as described
#' below.
#'
#' @details
#' * `dea` and `dea<-` are the methods to get and set the `dea` information as a
#' whole. These methods return `DeeDeeExperiment` objects.
#' * `add_dea` and `remove_dea` are used to respectively add or remove DE-results
#' items. These methods also return `DeeDeeExperiment` objects, with updated
#' content in the `dea` slot.
#' * `get_dea_df` and `get_dea_list` retrieve the `dea` information and provide
#' this as a `DataFrame` object (for a specific analysis) or as a list, with one
#' element for each reported analysis.
#' * `show` is the method to nicely print out the information of a `DeeDeeExperiment`
#' object.
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
#'
#' # creating a `DeeDeeExperiment`
#' dde <- DeeDeeExperiment(
#'   se_macrophage_noassays,
#'   de_results = de_named_list
#' )
#' dde
#'
#' new_del <- list(
#'   ifng2 = de_named_list$ifng_vs_naive,
#'   ifngsalmo2 = de_named_list$ifngsalmo_vs_naive
#' )
#'
#' # add a new (set of) DE result(s)
#' dde_new <- add_dea(dde, new_del)
#' dde_new
#'
#' # removing DEAs
#' dde_removed <- remove_dea(dde, "ifng_vs_naive")
#' dde_removed
NULL


# dea slot - get & set ---------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea",
          signature = "DeeDeeExperiment",
          definition = function(x) {
            x@dea
          })


#' @rdname DeeDeeExperiment-methods
#' @export
setReplaceMethod("dea",
                 c("DeeDeeExperiment", "list"),
                 definition = function(x, value) {
                   x@dea <- value
                   validObject(x)
                   x
                 })


# dea info - add, remove, get --------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("add_dea",
          signature = c("DeeDeeExperiment", "list"),
          definition = function(x, dea) {
            # dde must be a DeeDeeExp
            # dea must be named list

            # check that names are all unique, and do not overlap with the existing ones
            names(dea)
            names(dea(x))

            dea_contrasts <- dea(x)
            dde_ids <- rownames(x)

            # update rowData, naming them correctly
            for (i in names(dea)){
              this_de <- dea[[i]]

              # do different things according to what these objects are
              if(is(this_de, "DESeqResults")) {
                matched_ids <- match(rownames(x), dde_ids)

                # if not tested, add NA - everywhere? -> pre-fill?
                rowData(x)[[paste0(i,"_log2FoldChange")]] <- NA
                rowData(x)[[paste0(i,"_pvalue")]] <- NA
                rowData(x)[[paste0(i,"_padj")]] <- NA

                rowData(x)[[paste0(i,"_log2FoldChange")]][matched_ids] <- this_de$log2FoldChange
                rowData(x)[[paste0(i,"_pvalue")]][matched_ids] <- this_de$pvalue
                rowData(x)[[paste0(i,"_padj")]][matched_ids] <- this_de$padj

                dea_contrasts[[i]] <- list(
                  alpha = metadata(this_de)$alpha,
                  lfcThreshold = metadata(this_de)$lfcThreshold,
                  metainfo_logFC = mcols(this_de)$description[colnames(this_de) == "log2FoldChange"],
                  metainfo_pvalue = mcols(this_de)$description[colnames(this_de) == "pvalue"],
                  original_object = this_de,
                  package = "DESeq2"
                )
              }
            }

            # update the deslot
            dea(x) <- dea_contrasts

            # check here the validity
            validObject(x)

            # return the object
            return(x)
          }
)

# TODO: might need one where I also simply add ONE single DE object, and that gets autoconverted to a named list (of length 1)



#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("remove_dea",
          signature = c("DeeDeeExperiment", "character"),
          definition = function(x, dea_name) {
            # x must be a DeeDeeExp

            # dea must be char vector
            deas <- names(dea(x))

            deas_to_remove <- intersect(dea_name, deas)
            # warning() if nothing to remove

            for (i in deas_to_remove) {
              cols_to_remove <- c(paste0(i, c("_log2FoldChange", "_pvalue", "_padj")))
              rowData(x) <- rowData(x)[, !(colnames(rowData(x)) %in% cols_to_remove)]
              # update the deslot
              dea(x)[[i]] <- NULL
            }

            # here check some validity?
            validObject(x)

            # return the object
            return(x)
          }
)




#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("get_dea_df",
          signature = c("DeeDeeExperiment", "character"),
          definition = function(x,
                                dea_name) {
            deas <- dea(x)
            dea_names <- names(deas)

            if (!(dea_name %in% dea_names)) {
              stop("dea not found")
            }

            rd_info <- paste0(dea_name,
                              c("_log2FoldChange", "_pvalue", "_padj"))

            if (! all(rd_info %in% colnames(rowData(x)))) {
              stop("Columns not found")
            }

            out <- rowData(x)[, rd_info]

            return(out)
          }
)

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("get_dea_list",
          signature = c("DeeDeeExperiment"),
          definition = function(x) {
            deas <- dea(x)
            dea_names <- names(deas)

            dea_list <- list()

            for (i in dea_names) {
              dea_list[[i]] <- as.data.frame(get_dea_df(x, i))
              colnames(dea_list[[i]]) <- c("log2FoldChange", "pvalue", "padj")
            }

            return(dea_list)
          }
)




# misc - show & more ------------------------------------------------------

#' @name DeeDeeExperiment-misc
#'
#' @title Miscellaneous DeeDeeExperiment methods
#'
#' @description
#' Miscellaneous methods for the \code{\link{DeeDeeExperiment}} class and its
#' descendants that do not fit into any other documentation category such as,
#' for example, show methods.
#'
#' @param object a \code{\link{DeeDeeExperiment}} object
#'
#' @return Returns NULL
NULL


#' @rdname DeeDeeExperiment-misc
#' @export
setMethod("show",
          signature = signature(object = "DeeDeeExperiment"),
          definition = function(object) {

            callNextMethod()
            cat(
              "Including ", length(object@dea), " DE analyses:\n",
              paste(names(object@dea), collapse = ", "),
              sep=""
            )
          })
