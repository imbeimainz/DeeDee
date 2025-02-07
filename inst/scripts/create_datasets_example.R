# round 1 -----------------------------------------------------------------

library("macrophage")
library("DESeq2")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("dplyr")
library("ggplot2")

data(gse, "macrophage")
dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
dds_macrophage

se_macrophage <- gse
rownames(se_macrophage) <- substr(rownames(se_macrophage), 1, 15)
se_macrophage

se_macrophage_noassays <- SummarizedExperiment(
  assays = SimpleList(),
  rowData = rownames(se_macrophage)
)
rownames(se_macrophage_noassays) <- rownames(se_macrophage)

anno_df <- data.frame(
  gene_id = rownames(dds_macrophage),
  gene_name = mapIds(org.Hs.eg.db,
                     keys = rownames(dds_macrophage),
                     column = "SYMBOL",
                     keytype = "ENSEMBL"
  ),
  stringsAsFactors = FALSE,
  row.names = rownames(dds_macrophage)
)

colData(se_macrophage)
rowData(se_macrophage)

# DE run
keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
dds_unnormalized <- dds_macrophage



dds_macrophage <- DESeq(dds_macrophage)
vst_macrophage <- vst(dds_macrophage)
res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
                                        contrast = c("condition", "IFNg", "naive"),
                                        lfcThreshold = 1, alpha = 0.05
)
summary(res_macrophage_IFNg_vs_naive)
res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL

se_macrophage <- se_macrophage[keep, ]
se_macrophage_noassays <- se_macrophage_noassays[keep, ]

colData(se_macrophage)
rowData(se_macrophage)

# DE results
IFNg_naive <- results(dds_macrophage,
                      contrast = c("condition", "IFNg", "naive"),
                      lfcThreshold = 1, alpha = 0.05
)

save(IFNg_naive, file = "data/DE_results_IFNg_naive.RData", compress = "xz")


IFNg_both <- results(dds_macrophage,
                     contrast = c("condition", "IFNg_SL1344", "IFNg"),
                     lfcThreshold = 1, alpha = 0.05
)

save(IFNg_both, file = "data/DE_results_IFNg_both.RData", compress = "xz")


Salm_naive <- results(dds_macrophage,
                      contrast = c("condition", "SL1344", "naive"),
                      lfcThreshold = 1, alpha = 0.05
)

save(Salm_naive, file = "data/DE_results_Salm_naive.RData", compress = "xz")


Salm_both <- results(dds_macrophage,
                     contrast = c("condition", "IFNg_SL1344", "SL1344"),
                     lfcThreshold = 1, alpha = 0.05
)

save(Salm_both, file = "data/DE_results_Salm_both.RData", compress = "xz")


res_de <- res_macrophage_IFNg_vs_naive
str(res_de)

data(DE_results_IFNg_naive, package = "DeeDee")
IFNg_naive
data(DE_results_IFNg_both, package = "DeeDee")
IFNg_both
data(DE_results_Salm_naive, package = "DeeDee")
Salm_naive
data(DE_results_Salm_both, package = "DeeDee")
Salm_both

# for the original object/implementation:
# dd_list_original <- list(
#   IFNg_naive = DeeDeeLegacy::deedee_prepare(IFNg_naive, "DESeq2"),
#   IFNg_both = DeeDeeLegacy::deedee_prepare(IFNg_both, "DESeq2"),
#   Salm_naive = DeeDeeLegacy::deedee_prepare(Salm_naive, "DESeq2"),
#   Salm_both = DeeDeeLegacy::deedee_prepare(Salm_both, "DESeq2")
# )

# save(dd_list_original, file = "data/dd_list_original.RData", compress = "xz")
data("dd_list_original", package = "DeeDee")

# for the new version
de_named_list <- list(
  ifng_vs_naive = IFNg_naive,
  ifngsalmo_vs_naive = IFNg_both,
  salmonella_vs_naive = Salm_naive,
  salmo_both = Salm_both
)

# save(de_named_list, file = "data/de_named_list.RData", compress = "xz")
data("de_named_list", package = "DeeDee")

dde <- DeeDeeExperiment(se_macrophage_noassays, de_results = de_named_list)

