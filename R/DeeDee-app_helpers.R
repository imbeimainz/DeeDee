# nocov start
rds_input <- function(obj,
                      nm) {

  if (is(obj, "DeeDeeObject")) {
    obj <- obj@DeeDeeList
  }
  if (is(obj, "DESeqResults")) {
    ###### to-replace ####### obj <- DeeDeeLegacy::deedee_prepare(obj, "DESeq2")
    obj <- list(obj)
    names(obj) <- unlist(strsplit(nm,
                                  split = ".",
                                  fixed = TRUE
    ))[1]
  } else if (is(obj, "DGEExact")) {
    ###### to-replace ####### obj <- DeeDeeLegacy::deedee_prepare(obj, "edgeR")
    obj <- list(obj)
    names(obj) <- unlist(strsplit(nm,
                                  split = ".",
                                  fixed = TRUE
    ))[1]
  } else if (is(obj, "list")) {
    for (j in length(obj)) {
      if (checkmate::test_subset(
        names(obj[[j]]),
        c("logFC", "pval")
      ) == FALSE) {
        return(NULL)
      }
    }
  } else if (is(obj, "data.frame")) {
    if (length(obj) == 2) {
      if (checkmate::test_subset(
        names(obj),
        c("logFC", "pval")
      ) == FALSE) {
        return(NULL)
      }
      obj <- list(obj)
      names(obj) <- unlist(strsplit(nm,
                                    split = ".",
                                    fixed = TRUE
      ))[1]
    } else if (length(obj) == 6) {
      if (checkmate::test_subset(names(obj), c(
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B"
      )) == FALSE) {
        return(NULL)
      }
      ###### to-replace ####### obj <- DeeDeeLegacy::deedee_prepare(obj, "limma")
      obj <- list(obj)
      names(obj) <- unlist(strsplit(nm,
                                    split = ".",
                                    fixed = TRUE
      ))[1]
    } else {
      return(NULL)
    }
  }

  return(obj)
}


xlsx_input <- function(obj,
                       path,
                       nm) {

  if (length(obj) > 1) {
    out <- lapply(obj,
                  readxl::read_excel,
                  path = path
    )
    names(out) <- obj
    for (j in seq_len(length(obj))) {
      out[[obj[j]]] <- as.data.frame(out[[obj[j]]])
      out[[obj[j]]] <- tibble::column_to_rownames(
        out[[obj[j]]], "rowname"
      )
      if (checkmate::test_subset(
        names(out[[j]]) == FALSE,
        c("logFC", "pval")
      )) {
        return(NULL)
      }
    }
  } else {
    out <- readxl::read_excel(path = path)
    out <- tibble::column_to_rownames(out, "rowname")
    out <- list(out)
    names(out) <- unlist(strsplit(nm,
                                  split = ".",
                                  fixed = TRUE
    ))[1]
  }
  return(out)
}





ora <- function(geneList,
                universe,
                orgDB,
                key_type) {

  # ---------------------------- data preparation -----------------------------
  genes <- geneList[1]
  genes <- row.names(genes)
  genes <- genes[!is.na(genes)]

  universe[[1]] <- tibble::rownames_to_column(universe[[1]])
  universe[[2]] <- tibble::rownames_to_column(universe[[2]])
  univ <- dplyr::inner_join(universe[[1]], universe[[2]], by = "rowname")
  univ <- univ[["rowname"]]

  # ---------------------------------- GSEA -----------------------------------
  res <- clusterProfiler::enrichGO(
    gene = genes,
    universe = univ,
    OrgDb = get(orgDB),
    ont = "BP",
    keyType = key_type,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  # --------------------------------- return ----------------------------------
  return(res)
}



input_infobox <- function(deedee_obj,
                          sets,
                          md) {

  ext <- c()
  filename <- c()
  res <- list()
  type <- c()
  contrast <- c()
  genes <- c()
  count <- 0

  if (!is.null(deedee_obj)) {
    for (j in seq_len(length(deedee_obj@DeeDeeList))) {
      count <- count + 1
      type[count] <- "DeeDee object"
      filename[count] <- "input as argument"
      contrast[count] <- names(deedee_obj@DeeDeeList)[j]
      genes[count] <- length(deedee_obj@DeeDeeList[j]
                             [[contrast[count]]][["logFC"]])
    }
  }

  # reading out input files
  if (length(sets[, 1] > 0)) {
    for (i in seq_len((length(sets[, 1])))) {
      ext[i] <- tools::file_ext(sets[i, "datapath"])

      if (ext[[i]] == "rds" || ext[[i]] == "RDS") {
        res[[i]] <- readRDS(sets[[i, "datapath"]])
      } else if (ext[[i]] == "xlsx") {
        sheets <- readxl::excel_sheets(sets[[i, "datapath"]])
        res[[i]] <- lapply(sheets,
                           readxl::read_excel,
                           path = sets[[i, "datapath"]]
        )
        names(res[[i]]) <- sheets
        for (j in seq_len(length(sheets))) {
          res[[i]][[sheets[j]]] <- as.data.frame(res[[i]][[sheets[j]]])
          res[[i]][[sheets[j]]] <- tibble::column_to_rownames(
            res[[i]][[sheets[j]]], "rowname"
          )
        }
      } else if (ext[[i]] == "txt") {
        temp <- utils::read.table(sets[[i, "datapath"]])
        res[[i]] <- list(temp)
        names(res[[i]]) <- unlist(strsplit(sets[i, "name"],
                                           split = ".",
                                           fixed = TRUE
        ))[1]
      }

      if (is(res[[i]], "DESeqResults") ||
          is(res[[i]], "DGEExact") ||
          length(names(res[[i]])) == 6 ||
          is(res[[i]], "DeeDeeObject")) {
        count <- count + 1
        type[count] <- class(res[[i]])
        filename[count] <- sets[i, "name"]
        contrast[count] <- unlist(strsplit(filename[count],
                                           split = ".",
                                           fixed = TRUE
        ))[1]
        genes[count] <- length(md@DeeDeeList
                               [[contrast[count]]][["logFC"]])
      } else {
        if (is.data.frame(res[[i]])) {
          res[[i]] <- list(res[[i]])
          names(res[[i]]) <- unlist(strsplit(sets[i, "name"],
                                             split = ".",
                                             fixed = TRUE
          ))[1]
        }
        for (j in seq_len(length(res[[i]]))) {
          count <- count + 1
          type[count] <- "DeeDee object"
          filename[count] <- sets[i, "name"]
          contrast[count] <- names(res[[i]])[j]
          genes[count] <- length(res[[i]][j]
                                 [[contrast[count]]][["logFC"]])
        }
      }
    }
  }

  df <- data.frame(filename, type, contrast, genes)
  return(df)
}
# nocov end
