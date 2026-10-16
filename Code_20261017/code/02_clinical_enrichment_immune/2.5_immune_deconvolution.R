# ==== 2.5 免疫细胞比例推断 ====
# Execution scope: 历史脚本；需原始数据及资源
# See CHAPTERS.md for order, dependencies and entry points.

# Historical implementation, reformatted and repathed for disclosure.
# Not part of the tested CSV reproduction entry point.
# Read docs/HISTORICAL_SOURCE.md before enabling execution.
if (Sys.getenv("PCA_ENABLE_HISTORICAL") != "1")
  stop("Historical replay is disabled by default. Read docs/HISTORICAL_SOURCE.md first.")
CODE_ROOT <- Sys.getenv("PCA_CODE_ROOT", unset=normalizePath(".",winslash="/",mustWork=TRUE))
if(!file.exists(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R")))stop("Set PCA_CODE_ROOT to the release directory")
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"))

################################################################################
# Immune Microenvironment Analysis — 4 Genes × 4 Cohorts
# ESTIMATE scores + CIBERSORT + Correlation Heatmaps
# Output: BCR_ML_Publication/Immune_Results/
################################################################################

library(Biobase)
library(estimate)
library(e1071)
library(parallel)
library(preprocessCore)
library(corrplot)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd(legacy_output())

# Output goes to NEW folder
out_dir <- "Immune_Results"
dir.create(out_dir, showWarnings=FALSE)

source(project_path("data", "external", "CIBERSORT", "CIBERSORT.R"))

# Helper: convert Ensembl rownames to gene symbols
# ---- 2.5.1 ensembl_to_symbol ----
ensembl_to_symbol <- function(expr_mat) {
  ensembl_ids <- rownames(expr_mat)
  symbol_map <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = ensembl_ids,
                                      columns = "SYMBOL",
                                      keytype = "ENSEMBL")
  symbol_map <- symbol_map[!is.na(symbol_map$SYMBOL), , drop = FALSE]
  symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), , drop = FALSE]
  symbol_map <- symbol_map[!duplicated(symbol_map$SYMBOL), , drop = FALSE]

  expr_mat <- expr_mat[symbol_map$ENSEMBL, , drop = FALSE]
  rownames(expr_mat) <- symbol_map$SYMBOL[match(rownames(expr_mat), symbol_map$ENSEMBL)]
  return(expr_mat)
}

core_ensembl <- c("ENSG00000124721"="DNAH8", "ENSG00000131747"="TOP2A",
                  "ENSG00000101938"="CHRDL1", "ENSG00000154330"="PGM5")
genes <- unname(core_ensembl)

datasets <- list("DKFZ"=readRDS(project_path("data", "external", "modelling", "DKFZ_eSet.RDS")),
                 "TCGA"=readRDS(project_path("data", "external", "modelling", "TCGA-PRAD_eSet.RDS")),
                 "Taylor"=readRDS(project_path("data", "external", "modelling", "Taylor_eSet.RDS")),
                 "CancerMap"=readRDS(project_path("data", "external", "modelling", "CancerMap_eSet.RDS")))

# ===========================================================================
# ---- 2.5.2 Extract expression matrices (gene symbol) for each cohort ----
# ===========================================================================
cat("=== 1. Extracting expression data ===\n")

expr_list <- list()
for(ds in names(datasets)) {
  eset <- datasets[[ds]]
  expr <- exprs(eset)
  pd   <- pData(eset)

  bcr <- pd$bcr_status
  keep <- which(!is.na(bcr) & bcr %in% c(0,1))

  common <- intersect(names(core_ensembl), rownames(expr))
  expr_sub <- expr[common, keep, drop=FALSE]
  rownames(expr_sub) <- core_ensembl[common]

  expr_list[[ds]] <- expr_sub
  cat(sprintf("  %-12s: %d genes x %d samples\n", ds, nrow(expr_sub), ncol(expr_sub)))
}

# ===========================================================================
# ---- 2.5.3 ESTIMATE scores for each cohort ----
# ===========================================================================
cat("\n=== 2. Computing ESTIMATE scores ===\n")

tmp_dir <- file.path(out_dir, "tmp")
dir.create(tmp_dir, showWarnings=FALSE, recursive=TRUE)

estimate_scores <- list()

for(ds in names(expr_list)) {
  expr_mat <- expr_list[[ds]]

  eset <- datasets[[ds]]
  full_expr <- exprs(eset)
  bcr <- pData(eset)$bcr_status
  keep <- which(!is.na(bcr) & bcr %in% c(0,1))
  full_expr <- full_expr[, keep]

  full_expr_symbol <- ensembl_to_symbol(full_expr)
  cat(sprintf("  %-12s: %d Ensembl IDs -> %d gene symbols\n", ds, nrow(full_expr), nrow(full_expr_symbol)))

  tmp_file <- file.path(tmp_dir, paste0(ds, "_expr.txt"))
  write.table(data.frame(Gene=rownames(full_expr_symbol), full_expr_symbol, check.names=FALSE),
              tmp_file, sep="\t", row.names=FALSE, quote=FALSE)

  gct_file <- file.path(tmp_dir, paste0(ds, ".gct"))
  est_file <- file.path(tmp_dir, paste0(ds, "_estimate.txt"))

  tryCatch({
    filterCommonGenes(input.f=tmp_file, output.f=gct_file, id="GeneSymbol")
    estimateScore(gct_file, est_file, platform="illumina")

    est_res <- read.table(est_file, sep="\t", skip=2, header=TRUE, row.names=1, check.names=FALSE)
    est_res <- est_res[, -1, drop=FALSE]
    est_res <- as.data.frame(t(est_res))

    for(cn in colnames(est_res)) est_res[[cn]] <- as.numeric(est_res[[cn]])

    estimate_scores[[ds]] <- est_res
    cat(sprintf("  %-12s: done\n", ds))
  }, error=function(e) {
    cat(sprintf("  %-12s: FAILED - %s\n", ds, e$message))
  })
}

saveRDS(estimate_scores, file.path(out_dir, "estimate_scores.RDS"))

# ===========================================================================
# ---- 2.5.4 CIBERSORT immune deconvolution ----
# ===========================================================================
cat("\n=== 3. Running CIBERSORT ===\n")
sig_matrix <- project_path("data", "external", "CIBERSORT", "LM22.txt")

cibersort_res <- list()

for(ds in names(expr_list)) {
  full_expr <- exprs(datasets[[ds]])
  bcr <- pData(datasets[[ds]])$bcr_status
  keep <- which(!is.na(bcr) & bcr %in% c(0,1))
  full_expr <- full_expr[, keep]

  full_expr_symbol <- ensembl_to_symbol(full_expr)

  mixture_file <- file.path(tmp_dir, paste0(ds, "_cibersort_input.txt"))
  write.table(data.frame(Gene=rownames(full_expr_symbol), full_expr_symbol, check.names=FALSE),
              mixture_file, sep="\t", row.names=FALSE, quote=FALSE)

  tryCatch({
    res <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
    res <- res[, 1:22]
    res <- res[, colSums(res) > 0, drop=FALSE]
    cibersort_res[[ds]] <- as.data.frame(res)
    cat(sprintf("  %-12s: %d cell types\n", ds, ncol(res)))
  }, error=function(e) {
    cat(sprintf("  %-12s: FAILED - %s\n", ds, e$message))
  })
}

saveRDS(cibersort_res, file.path(out_dir, "cibersort_res.RDS"))

# ===========================================================================
# ---- 2.5.5 Correlation: 4 genes vs ESTIMATE scores ----
# ===========================================================================
cat("\n=== 4. Gene-ESTIMATE correlation heatmaps ===\n")

for(ds in names(estimate_scores)) {
  est <- estimate_scores[[ds]]
  expr_mat <- expr_list[[ds]]

  common_samp <- intersect(colnames(expr_mat), rownames(est))
  if(length(common_samp) < 10) next

  available_genes <- intersect(genes, rownames(expr_mat))
  if(length(available_genes) < 2) next
  expr_4g <- t(expr_mat[available_genes, common_samp, drop=FALSE])
  est_sub <- est[common_samp, , drop=FALSE]

  cor_mat <- cor(expr_4g, est_sub, method="spearman", use="pairwise.complete.obs")

  pdf(file.path(out_dir, paste0("Fig13_Cor_ESTIMATE_", ds, ".pdf")), width=6, height=5)
  corrplot(cor_mat, method="color", type="full",
           col=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
           addCoef.col="black", number.cex=0.8,
           tl.col="black", tl.cex=1, cl.cex=0.8,
           title=paste0("Gene vs ESTIMATE — ", ds),
           mar=c(2,0,3,0))
  dev.off()

  png(file.path(out_dir, paste0("Fig13_Cor_ESTIMATE_", ds, ".png")), width=6*150, height=5*150, res=150)
  corrplot(cor_mat, method="color", type="full",
           col=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
           addCoef.col="black", number.cex=0.8,
           tl.col="black", tl.cex=1, cl.cex=0.8,
           title=paste0("Gene vs ESTIMATE — ", ds),
           mar=c(2,0,3,0))
  dev.off()

  write.csv(cor_mat, file.path(out_dir, paste0("Cor_ESTIMATE_", ds, ".csv")))
  cat(sprintf("  %-12s: saved\n", ds))
}

# ===========================================================================
# ---- 2.5.6 Correlation: 4 genes vs CIBERSORT immune cells ----
# ===========================================================================
cat("\n=== 5. Gene-CIBERSORT correlation heatmaps ===\n")

for(ds in names(cibersort_res)) {
  ciber <- cibersort_res[[ds]]
  expr_mat <- expr_list[[ds]]

  common_samp <- intersect(colnames(expr_mat), rownames(ciber))
  if(length(common_samp) < 10) next

  available_genes <- intersect(genes, rownames(expr_mat))
  if(length(available_genes) < 2) next
  expr_4g <- t(expr_mat[available_genes, common_samp, drop=FALSE])
  ciber_sub <- ciber[common_samp, , drop=FALSE]

  cor_mat <- cor(expr_4g, ciber_sub, method="spearman", use="pairwise.complete.obs")

  pdf(file.path(out_dir, paste0("Fig14_Cor_CIBERSORT_", ds, ".pdf")), width=10, height=5)
  corrplot(cor_mat, method="color", type="full",
           col=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
           addCoef.col="black", number.cex=0.55,
           tl.col="black", tl.cex=0.7, tl.srt=45, cl.cex=0.8,
           title=paste0("Gene vs Immune Cells — ", ds),
           mar=c(2,0,3,0))
  dev.off()

  png(file.path(out_dir, paste0("Fig14_Cor_CIBERSORT_", ds, ".png")), width=10*150, height=5*150, res=150)
  corrplot(cor_mat, method="color", type="full",
           col=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
           addCoef.col="black", number.cex=0.55,
           tl.col="black", tl.cex=0.7, tl.srt=45, cl.cex=0.8,
           title=paste0("Gene vs Immune Cells — ", ds),
           mar=c(2,0,3,0))
  dev.off()

  write.csv(cor_mat, file.path(out_dir, paste0("Cor_CIBERSORT_", ds, ".csv")))
  cat(sprintf("  %-12s: saved\n", ds))
}

# ===========================================================================
# ---- 2.5.7 Summary ----
# ===========================================================================
cat("\n=== Immune Analysis Complete ===\n")
cat("Output:", out_dir, "/\n")
cat("  Fig13_Cor_ESTIMATE_*.pdf  — 4 genes vs Stromal/Immune/ESTIMATE scores\n")
cat("  Fig14_Cor_CIBERSORT_*.pdf — 4 genes vs 22 immune cell types\n")
cat("  estimate_scores.RDS, cibersort_res.RDS\n")
cat("  Cor_ESTIMATE_*.csv, Cor_CIBERSORT_*.csv\n")
