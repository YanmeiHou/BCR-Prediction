# ==== 1.1 差异表达、LASSO与随机森林筛选 ====
# Execution scope: 历史脚本；需原始数据
# See CHAPTERS.md for order, dependencies and entry points.

# Historical implementation, reformatted and repathed for disclosure.
# Not part of the tested CSV reproduction entry point.
# Read docs/HISTORICAL_SOURCE.md before enabling execution.
if (Sys.getenv("PCA_ENABLE_HISTORICAL") != "1")
  stop("Historical replay is disabled by default. Read docs/HISTORICAL_SOURCE.md first.")
CODE_ROOT <- Sys.getenv("PCA_CODE_ROOT", unset=normalizePath(".",winslash="/",mustWork=TRUE))
if(!file.exists(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R")))stop("Set PCA_CODE_ROOT to the release directory")
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"))

# ============================================================
# TCGA BCR vs non-BCR 差异基因 + 随机森林(RF) + LASSO 特征筛选
# Output: outputs/historical_replay/Feature_discovery/
#
# 流程:
# ---- 1.1.1 limma 差异分析 (protein_coding, |log2FC|>=0.5, adj.P<0.05) ----
# ---- 1.1.2 randomForest 变量重要性排序 ----
# ---- 1.1.3 glmnet LASSO (10-fold CV, lambda.min) ----
# ---- 1.1.4 RF 与 LASSO 取交集 -> 关键基因 ----
# ============================================================

suppressMessages({
  library(Biobase); library(limma)
  library(glmnet); library(randomForest)
  library(ggplot2); library(ggrepel)
  library(VennDiagram); library(grid)
})

outdir <- file.path(legacy_output(), "Feature_discovery")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ---- 1.1.5 注释 ----
ann <- readRDS(project_path("data", "external", "discovery", "PCaDB_Gene_Annotation.RDS"))
ann_clean <- subset(ann, gene_biotype == "protein_coding" & !is.na(gene_name) & gene_name != "")
pc_ens <- unique(ann_clean$ensembl_id)
symbol_map <- setNames(ann_clean$gene_name, ann_clean$ensembl_id)

# ---- 1.1.6 载入 TCGA-PRAD ----
eset <- readRDS(project_path("data", "external", "discovery", "TCGA-PRAD_eSet.RDS"))
pd <- pData(eset)
keep <- which(pd$bcr_status %in% c("0", "1") & !is.na(pd$bcr_status))
esub <- eset[, keep]
group <- factor(pData(esub)$bcr_status, levels = c("0", "1"))
expr <- exprs(esub)
if (max(expr, na.rm = TRUE) > 100) expr <- log2(expr + 1)
y <- group

cat("=== 样本统计 ===\n")
cat("non-BCR:", sum(group == "0"), " BCR:", sum(group == "1"), "\n")

# ---- 1.1.7 差异表达分析 ----
design <- model.matrix(~ group)
fit <- lmFit(expr, design)
fit <- eBayes(fit, trend = TRUE)
res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
res$gene_symbol <- symbol_map[rownames(res)]
res$ensembl_id <- rownames(res)
res <- res[rownames(res) %in% pc_ens, , drop = FALSE]

keepdeg <- which(res$adj.P.Val < 0.05 & abs(res$logFC) >= 0.5)
degs <- res[keepdeg, , drop = FALSE]
degs <- degs[order(degs$adj.P.Val), , drop = FALSE]
degs$direction <- ifelse(degs$logFC > 0, "Up", "Down")
cat("=== 差异基因 (|log2FC|>=0.5 & adj.P<0.05) ===\n")
cat("DEGs:", nrow(degs), " | Up:", sum(degs$direction == "Up"),
    " | Down:", sum(degs$direction == "Down"), "\n")

write.csv(degs[, c("ensembl_id", "gene_symbol", "logFC", "AveExpr", "t",
                   "P.Value", "adj.P.Val", "B", "direction")],
          file.path(outdir, "TCGA_BCR_vs_NonBCR_DEGs.csv"), row.names = FALSE)

# ---- 1.1.8 表达矩阵 (样本 x 基因) ----
X <- t(expr[rownames(degs), , drop = FALSE])
sym <- degs$gene_symbol[match(colnames(X), rownames(degs))]
sym[is.na(sym)] <- colnames(X)[is.na(sym)]
colnames(X) <- make.unique(sym)

# 项目已确立的 4 个关键基因 (4 数据集共享 DEG)
key4 <- c("DNAH8", "TOP2A", "CHRDL1", "PGM5")
key4_ens <- c("ENSG00000124721", "ENSG00000131747",
              "ENSG00000101938", "ENSG00000154330")
names(key4_ens) <- key4

# ---- 1.1.9 LASSO ----
set.seed(123)
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10)
coef_min <- as.matrix(coef(cvfit, s = "lambda.min"))
coef_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
lasso_min <- setdiff(rownames(coef_min)[coef_min[, 1] != 0], "(Intercept)")
lasso_1se <- setdiff(rownames(coef_1se)[coef_1se[, 1] != 0], "(Intercept)")
cat("=== LASSO ===\n")
cat("lambda.min:", cvfit$lambda.min, " genes:", length(lasso_min), "\n")
cat("lambda.1se:", cvfit$lambda.1se, " genes:", length(lasso_1se), "\n")

write.csv(data.frame(gene = lasso_min), file.path(outdir, "LASSO_lambda.min_genes.csv"),
          row.names = FALSE)
write.csv(data.frame(gene = lasso_1se), file.path(outdir, "LASSO_lambda.1se_genes.csv"),
          row.names = FALSE)

# LASSO 系数表
coef_tab <- data.frame(
  gene = rownames(coef_min),
  coefficient = coef_min[, 1],
  stringsAsFactors = FALSE
)
coef_tab <- coef_tab[coef_tab$gene != "(Intercept)" & coef_tab$coefficient != 0, , drop = FALSE]
coef_tab <- coef_tab[order(-abs(coef_tab$coefficient)), , drop = FALSE]
write.csv(coef_tab, file.path(outdir, "LASSO_coefficients.csv"), row.names = FALSE)

# ---- 1.1.10 随机森林 ----
set.seed(123)
rf <- randomForest(x = X, y = y, ntree = 500, importance = TRUE)
mda <- importance(rf, type = 1)[, 1]
gini <- importance(rf, type = 2)[, 1]
imp_df <- data.frame(gene = names(mda), MeanDecreaseAccuracy = mda,
                     MeanDecreaseGini = gini, stringsAsFactors = FALSE)
imp_df <- imp_df[order(-imp_df$MeanDecreaseAccuracy), , drop = FALSE]
write.csv(imp_df, file.path(outdir, "RF_importance.csv"), row.names = FALSE)

# RF 筛选: 重要性 Top 20 (常用阈值)
rf_top20 <- imp_df$gene[1:20]
write.csv(data.frame(gene = rf_top20), file.path(outdir, "RF_top20_genes.csv"),
          row.names = FALSE)

# ---- 1.1.11 RF 与 LASSO 交集 ----
inter_min <- intersect(rf_top20, lasso_min)
inter_1se <- intersect(rf_top20, lasso_1se)
cat("=== RF(top20) 与 LASSO 交集 ===\n")
cat("RF_top20 ∩ LASSO.min =", length(inter_min), ":", paste(inter_min, collapse = ","), "\n")
cat("RF_top20 ∩ LASSO.1se =", length(inter_1se), ":", paste(inter_1se, collapse = ","), "\n")

write.csv(data.frame(gene = inter_min), file.path(outdir, "Key_genes_RF_LASSO_intersection.csv"),
          row.names = FALSE)

# ---- 1.1.12 火山图 ----
vol <- res
vol$logP <- -log10(vol$adj.P.Val)
vol$logP[is.infinite(vol$logP)] <- max(vol$logP[is.finite(vol$logP)], na.rm = TRUE)
vol$sig <- "NS"
vol$sig[vol$adj.P.Val < 0.05 & vol$logFC >= 0.5] <- "Up"
vol$sig[vol$adj.P.Val < 0.05 & vol$logFC <= -0.5] <- "Down"
# 只高亮并标注 4 个关键基因
vol_highlight <- vol[vol$gene_symbol %in% key4, , drop = FALSE]

cols <- c("Up" = "#C62828", "Down" = "#1565C0", "NS" = "grey80")
p <- ggplot(vol, aes(x = logFC, y = logP, color = sig)) +
  geom_point(size = 0.4, alpha = 0.55) +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  # 高亮 4 个关键基因的点
  geom_point(data = vol_highlight, aes(x = logFC, y = logP),
             color = "black", size = 2.5, shape = 21, fill = "red", stroke = 0.6) +
  # 只对 4 个关键基因标注；max.overlaps = Inf 保证 4 个标签全部显示
  geom_text_repel(data = vol_highlight, aes(x = logFC, y = logP, label = gene_symbol),
                  size = 4.5, max.overlaps = Inf, box.padding = 0.6,
                  segment.size = 0.4, fontface = "bold.italic", color = "black",
                  min.segment.length = 0) +
  labs(title = "TCGA-PRAD: BCR vs non-BCR (protein-coding DEGs)",
       x = expression(log[2] ~ "Fold Change"),
       y = expression(-log[10] ~ "(adjusted P)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(outdir, "plots", "Volcano_TCGA_BCR.pdf"), p, width = 8, height = 7)
ggsave(file.path(outdir, "plots", "Volcano_TCGA_BCR.png"), p, width = 8, height = 7, dpi = 300)

# ---- 1.1.13 RF 重要性图 (Top 20) ----
top20_df <- imp_df[1:20, , drop = FALSE]
top20_df <- top20_df[order(top20_df$MeanDecreaseAccuracy), , drop = FALSE]
bar_cols <- ifelse(top20_df$gene %in% key4, "#C62828",
            ifelse(top20_df$gene %in% inter_min, "#2E7D32", "#5B8FF9"))

pdf(file.path(outdir, "plots", "RF_importance_top20.pdf"), width = 9, height = 7)
par(mar = c(5, 10, 4, 2))
barplot(top20_df$MeanDecreaseAccuracy, horiz = TRUE, names.arg = top20_df$gene,
        las = 2, col = bar_cols, border = NA, cex.names = 0.85,
        xlab = "Mean Decrease Accuracy", main = "Random Forest: Top 20 Genes by Importance")
abline(v = 0, col = "grey40")
legend("bottomright", legend = c("Established 4-gene panel", "RF\u2229LASSO genes", "Other"),
       fill = c("#C62828", "#2E7D32", "#5B8FF9"), bty = "n", cex = 0.8)
dev.off()
png(file.path(outdir, "plots", "RF_importance_top20.png"), width = 9, height = 7,
    units = "in", res = 300)
par(mar = c(5, 10, 4, 2))
barplot(top20_df$MeanDecreaseAccuracy, horiz = TRUE, names.arg = top20_df$gene,
        las = 2, col = bar_cols, border = NA, cex.names = 0.85,
        xlab = "Mean Decrease Accuracy", main = "Random Forest: Top 20 Genes by Importance")
abline(v = 0, col = "grey40")
legend("bottomright", legend = c("Established 4-gene panel", "RF\u2229LASSO genes", "Other"),
       fill = c("#C62828", "#2E7D32", "#5B8FF9"), bty = "n", cex = 0.8)
dev.off()

# ---- 1.1.14 LASSO 图 ----
pdf(file.path(outdir, "plots", "LASSO_lambda_path.pdf"), width = 8, height = 6)
plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE,
     main = "LASSO Coefficient Path (TCGA BCR vs non-BCR)")
dev.off()
png(file.path(outdir, "plots", "LASSO_lambda_path.png"), width = 8, height = 6,
    units = "in", res = 300)
plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE,
     main = "LASSO Coefficient Path (TCGA BCR vs non-BCR)")
dev.off()

pdf(file.path(outdir, "plots", "LASSO_CV.pdf"), width = 7, height = 6)
plot(cvfit, main = "LASSO 10-fold Cross-Validation")
abline(v = log(cvfit$lambda.min), col = "red", lty = 2)
abline(v = log(cvfit$lambda.1se), col = "blue", lty = 2)
legend("topright", legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"), lty = 2, bty = "n")
dev.off()
png(file.path(outdir, "plots", "LASSO_CV.png"), width = 7, height = 6,
    units = "in", res = 300)
plot(cvfit, main = "LASSO 10-fold Cross-Validation")
abline(v = log(cvfit$lambda.min), col = "red", lty = 2)
abline(v = log(cvfit$lambda.1se), col = "blue", lty = 2)
legend("topright", legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"), lty = 2, bty = "n")
dev.off()

# ---- 1.1.15 Venn 图 ----
venn_list <- list(RF_top20 = rf_top20, LASSO = lasso_min)
pdf(file.path(outdir, "plots", "Venn_RF_LASSO.pdf"), width = 7, height = 7)
grid.draw(venn.diagram(
  x = venn_list, category.names = c("Random Forest\n(Top 20)", "LASSO\n(lambda.min)"),
  filename = NULL, col = c("#1F77B4", "#D62728"), fill = c("#1F77B4", "#D62728"),
  alpha = 0.4, cex = 1.4, cat.cex = 1.2, cat.col = c("#1F77B4", "#D62728"),
  margin = 0.1, main = "RF \u2229 LASSO"))
dev.off()
png(file.path(outdir, "plots", "Venn_RF_LASSO.png"), width = 7, height = 7,
    units = "in", res = 300)
grid.draw(venn.diagram(
  x = venn_list, category.names = c("Random Forest\n(Top 20)", "LASSO\n(lambda.min)"),
  filename = NULL, col = c("#1F77B4", "#D62728"), fill = c("#1F77B4", "#D62728"),
  alpha = 0.4, cex = 1.4, cat.cex = 1.2, cat.col = c("#1F77B4", "#D62728"),
  margin = 0.1, main = "RF \u2229 LASSO"))
dev.off()

# ---- 1.1.16 汇总 RDS ----
saveRDS(list(
  degs = degs, lasso_min = lasso_min, lasso_1se = lasso_1se,
  rf_importance = imp_df, rf_top20 = rf_top20,
  intersection_min = inter_min, intersection_1se = inter_1se,
  cvfit = cvfit, rf = rf, key4 = key4, key4_ens = key4_ens,
  n_degs = nrow(degs)
), file.path(outdir, "results_summary.RDS"))

cat("\n=== 完成，输出目录:", outdir, " ===\n")
