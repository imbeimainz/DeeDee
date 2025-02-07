suppressPackageStartupMessages(
  library("SummarizedExperiment")
)

data("de_named_list", package = "DeeDee")

rd_macrophage <- DataFrame(
  gene_id = rownames(de_named_list$ifng_vs_naive))
rownames(rd_macrophage) <- rownames(de_named_list$ifng_vs_naive)
se_macrophage_noassays <- SummarizedExperiment(
  assays = SimpleList(),
  rowData = rd_macrophage
)

names(de_named_list)

