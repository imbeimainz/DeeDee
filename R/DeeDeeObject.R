#' DeeDeeObject (S4 class)
#'
#' @slot DeeDeeList A DeeDeeList object
#'
#' @description DeeDeeObject is the format that the DeeDee App works on. It
#'              holds a slot for your list of DeeDee tables. More functionality
#'              to come.
#' @export
#'
#' @examples
#'
#' data(DE_results_IFNg_naive, package = "DeeDee")
#' IFNg_naive <- deedee_prepare(IFNg_naive, "DESeq2")
#'
#' data(DE_results_IFNg_both, package = "DeeDee")
#' IFNg_both <- deedee_prepare(IFNg_both, "DESeq2")
#'
#' data(DE_results_Salm_naive, package = "DeeDee")
#' Salm_naive <- deedee_prepare(Salm_naive, "DESeq2")
#'
#' data(DE_results_Salm_both, package = "DeeDee")
#' Salm_both <- deedee_prepare(Salm_both, "DESeq2")
#'
#' dd_list <- list(
#'   IFNg_naive = IFNg_naive, IFNg_both = IFNg_both,
#'   Salm_naive = Salm_naive, Salm_both = Salm_both
#' )
#'
#' # obj <- DeeDeeObject(DeeDeeList = dd_list)
#'

DeeDeeObject <- setClass(

  "DeeDeeObject",
  slots = c(
    DeeDeeList = "list"
  ),
  prototype = list(
    DeeDeeList = NULL
  ),
  validity = function(object) {
    res <- object@DeeDeeList
    if (class(res) != "list") {
      return("DeeDeeList must be a list of DeeDee data.frames.")
    }
    for (i in 1:length(res)) {
      if (class(res[[i]]) != "data.frame") {
        return("DeeDeeList must be a list of DeeDee data.frames.")
      }
      if (checkmate::test_subset(names(res[[i]]), c("logFC", "pval")) == FALSE) {
        return("DeeDeeList must be a list of DeeDee data.frames.")
      }
      if (class(res[[i]]$logFC) != "numeric") {
        return("LogFC values must be a vector of type 'numeric'.")
      }
      if (class(res[[i]]$pval) != "numeric") {
        return("P-values must be a vector of type 'numeric'.")
      }
      if (max(res[[i]]$pval) > 1) {
        return("A p-value cannot be greater than 1.")
      }
      if (min(res[[i]]$pval) < 0) {
        return("A p-value cannot be lower than 0.")
      }
    }
    return(TRUE)
  }
)
