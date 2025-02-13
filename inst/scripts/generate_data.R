library("macrophage")
library("DESeq2")
source("~/Development/DeeDee_wip/deedee_prepare.R")
library("SummarizedExperiment")
library("edgeR")
library("limma")

# load macrophage dataset
data(gse, "macrophage")

# DE with DESeq2 ----------------------------------------------------------------

dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)

keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
dds_macrophage <- DESeq(dds_macrophage)

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

# DE with limma ----------------------------------------------------------------

# extract counts and sample metadata
counts <- assays(gse)$counts
sample_info <- colData(gse)

# extract the conditions and cell line info
condition <- factor(sample_info$condition)
line <- factor(sample_info$line)

# create design for DE
design <- model.matrix(~ line + condition)

# create DGE list
dge <- DGEList(counts)

# normalize the counts
dge <- calcNormFactors(dge)

# filter by expression
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes=TRUE]

# transform counts into logCPM
v <- voom(dge, design)

# fitting linear models using weighted least squares for each gene
fit <- lmFit(v, design)

# setup comparisons
contrast_matrix <- makeContrasts(
  IFNgNaive     = conditionIFNg,
  IFNg_both = conditionIFNg_SL1344 - conditionIFNg,
  SalmNaive     = conditionSL1344,
  Salm_both = conditionIFNg_SL1344 - conditionSL1344,
  levels = colnames(design)
)

# apply contrast
fit2 <- contrasts.fit(fit, contrast_matrix)

# empirical Bayes smoothing of standard errors
fit2 <- eBayes(fit2)
de_limma <- fit2 # MArrayLM object

# save
save(de_limma, file = "data/de_limma.RData", compress = "xz")

# DE with edgeR ----------------------------------------------------------------

# extract counts and sample metadata
counts <- assays(gse)$counts
sample_info <- colData(gse)

# extract the conditions and cell line info
condition <- factor(sample_info$condition)
line <- factor(sample_info$line)

# create DGEList object
dge <- DGEList(counts=counts, group=condition)

# normalize the counts
dge <- calcNormFactors(dge)

# filter by expression
keep <- filterByExpr(dge, group=condition)
dge <- dge[keep, , keep.lib.sizes=TRUE]

# create design
design <- model.matrix(~ line + condition)

# estimate dispersion
dge <- estimateDisp(dge, design)

# perform likelihood ratio test
fit <- glmFit(dge, design)
dge_lrt <- glmLRT(fit, coef=7:9) # DGELRT object

# perform exact test
# exact test doesn't handle multi factor models, so we have to subset
# IFNg vs naive
keep_samples <- sample_info$condition %in% c("naive", "IFNg")
dge_sub <- dge[, keep_samples]
# droplevel
dge_sub$samples$group <- factor(sample_info$condition[keep_samples])
# renormalizw
dge_sub <- calcNormFactors(dge_sub)
dge_exact_IFNg_naive <- exactTest(dge_sub, pair=c("naive", "IFNg")) # DGEExact object

# SL1344 vs naive
keep_samples <- sample_info$condition %in% c("naive", "SL1344")
dge_sub <- dge[, keep_samples]
dge_sub$samples$group <- factor(sample_info$condition[keep_samples])
dge_sub <- calcNormFactors(dge_sub)
dge_exact_Salm_naive <- exactTest(dge_sub, pair=c("naive", "SL1344"))

# IFNg_SL1344 vs IFNg
keep_samples <- sample_info$condition %in% c("IFNg", "IFNg_SL1344")
dge_sub <- dge[, keep_samples]
dge_sub$samples$group <- factor(sample_info$condition[keep_samples])
dge_sub <- calcNormFactors(dge_sub)
dge_exact_IFNg_both <- exactTest(dge_sub, pair=c("IFNg", "IFNg_SL1344"))

# IFNg_SL1344 vs SL1344
keep_samples <- sample_info$condition %in% c("SL1344", "IFNg_SL1344")
dge_sub <- dge[, keep_samples]
dge_sub$samples$group <- factor(sample_info$condition[keep_samples])
dge_sub <- calcNormFactors(dge_sub)
dge_exact_Salm_both <- exactTest(dge_sub, pair=c("SL1344", "IFNg_SL1344"))

# save
save(dge_lrt, file = "data/DGELRT_macrophage.RData", compress = "xz")
save(dge_exact_IFNg_naive, file = "data/DGEExact_IFNg_naive.RData", compress = "xz")
save(dge_exact_Salm_naive, file = "data/DGEExact_Salm_naive.RData", compress = "xz")
save(dge_exact_IFNg_both, file = "data/DGEExact_IFNg_both.RData", compress = "xz")
save(dge_exact_Salm_both, file = "data/DGEExact_Salm_both.RData", compress = "xz")

