# ==== 1.3 模型训练、比较与验证 ====
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

################################################################################
# BCR 4-Gene ML Prediction — Historical Benchmark Analysis
#
# Training:   TCGA-PRAD (n=547)
# Validation: DKFZ (n=105), Taylor (n=140), CancerMap (n=233)
# Genes:      DNAH8, TOP2A, CHRDL1, PGM5
#
# Output:     BCR_ML_Publication/
################################################################################

# ===========================================================================
# ---- 1.3.1 Setup ----
# ===========================================================================
library(Biobase)
library(openxlsx)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(MASS)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(grid)
library(ggplot2)
library(ggpubr)
library(survival)

source(project_path("code/01_gene_selection_ml/1.2_ML_functions.R"))
setwd(legacy_output())

# ===========================================================================
# ---- 1.3.2 Data Preparation ----
# ===========================================================================
cat("═══════════════════════════════════════════\n")
cat("  BCR 4-Gene Prediction Model — Publication Analysis\n")
cat("═══════════════════════════════════════════\n\n")

core_ensembl <- c(
  "ENSG00000124721" = "DNAH8",
  "ENSG00000131747" = "TOP2A",
  "ENSG00000101938" = "CHRDL1",
  "ENSG00000154330" = "PGM5"
)

eset_files <- c(
  "DKFZ"      = project_path("data", "external", "modelling", "DKFZ_eSet.RDS"),
  "TCGA"      = project_path("data", "external", "modelling", "TCGA-PRAD_eSet.RDS"),
  "Taylor"    = project_path("data", "external", "modelling", "Taylor_eSet.RDS"),
  "CancerMap" = project_path("data", "external", "modelling", "CancerMap_eSet.RDS")
)

all_expr_list  <- list()
all_label_list <- list()

for (nm in names(eset_files)) {
  eset <- readRDS(eset_files[nm])
  expr <- exprs(eset)
  pd   <- pData(eset)

  common <- intersect(names(core_ensembl), rownames(expr))
  X <- t(expr[common, , drop = FALSE])
  colnames(X) <- core_ensembl[common]

  bcr  <- pd$bcr_status
  keep <- which(bcr %in% c(0, 1))

  all_expr_list[[nm]]  <- as.data.frame(X[keep, , drop = FALSE])
  all_label_list[[nm]] <- data.frame(
    outcome = factor(bcr[keep], levels = c("0", "1")),
    Cohort  = nm,
    row.names = rownames(X[keep, , drop = FALSE])
  )
}

all_expr  <- do.call(rbind, all_expr_list)
all_label <- do.call(rbind, all_label_list)

cat("--- Data Summary ---\n")
cat("Total samples:", nrow(all_expr), "\n")
cat("BCR events:   ", sum(all_label$outcome == "1"), "\n")
for(cn in names(eset_files)) {
  n <- sum(all_label$Cohort == cn)
  e <- sum(all_label$Cohort == cn & all_label$outcome == "1")
  cat(sprintf("  %-12s n=%d, BCR=%d\n", cn, n, e))
}

saveRDS(list(expr = all_expr, label = all_label, genes = core_ensembl),
        "prepared_data.RDS")

# ===========================================================================
# ---- 1.3.3 Train/Test Split ----
# ===========================================================================
train_cohort <- "TCGA"
test_cohorts <- c("DKFZ", "Taylor", "CancerMap")
all_cohorts  <- c(train_cohort, test_cohorts)

train_idx <- which(all_label$Cohort == train_cohort)
test_idx  <- which(all_label$Cohort %in% test_cohorts)

Train_expr  <- all_expr[train_idx, , drop = FALSE]
Train_class <- all_label[train_idx, , drop = FALSE]
Test_expr   <- all_expr[test_idx, , drop = FALSE]
Test_class  <- all_label[test_idx, , drop = FALSE]

Train_set <- scaleData(Train_expr, centerFlags = TRUE, scaleFlags = TRUE)
Test_set  <- scaleData(Test_expr, cohort = Test_class$Cohort,
                       centerFlags = TRUE, scaleFlags = TRUE)
all_combined <- rbind.data.frame(Train_set, Test_set)
all_label_combined <- rbind(Train_class, Test_class)

classVar <- "outcome"

cat(sprintf("\nTraining:   %s (n=%d, BCR=%d)\n",
            train_cohort, nrow(Train_set), sum(Train_class$outcome == "1")))
cat(sprintf("Validation:  %d samples (BCR=%d)\n\n",
            nrow(Test_set), sum(Test_class$outcome == "1")))

# ===========================================================================
# ---- 1.3.4 Load Methods & Train Models ----
# ===========================================================================
methods_df <- read.xlsx(project_path("data", "derived", "methods.xlsx"), startRow = 2)
methods <- na.omit(methods_df$Model)
methods <- gsub("-| ", "", methods)
methods <- methods[methods != ""]
cat("Total ML algorithms:", length(methods), "\n")

if (file.exists("model.rds")) {
  model <- readRDS("model.rds")
  cat("Loaded cached models:", length(model), "\n")
} else {
  preTrain.var <- list()
  preTrain.var[["simple"]] <- colnames(Train_set)
  model <- list()
  set.seed(777)
  Train_set_bk <- Train_set

  for (method in methods) {
    method_name <- method
    method_split <- strsplit(method, "\\+")[[1]]

    if (length(method_split) == 1) {
      var_method <- "simple"; fit_method <- method_split[1]
    } else {
      var_method <- method_split[1]; fit_method <- method_split[2]
    }

    if (!var_method %in% names(preTrain.var)) {
      tryCatch({
        preTrain.var[[var_method]] <- RunML(var_method, Train_set, Train_class,
                                            mode = "Variable", classVar = classVar)
      }, error = function(e) { preTrain.var[[var_method]] <<- colnames(Train_set) })
    }
    Variable <- preTrain.var[[var_method]]
    if (length(Variable) < 2) Variable <- colnames(Train_set)

    tryCatch({
      model[[method_name]] <- RunML(fit_method,
                                    Train_set_bk[, Variable, drop = FALSE],
                                    Train_class, mode = "Model", classVar = classVar)
      if (length(ExtractVar(model[[method_name]])) < 2) model[[method_name]] <- NULL
    }, error = function(e) { model[[method_name]] <<- NULL })
  }
  Train_set <- Train_set_bk
  saveRDS(model, "model.rds")
  cat("Trained", length(model), "models\n")
}

# ===========================================================================
# ---- 1.3.5 Compute Risk Scores & AUC ----
# ===========================================================================
cat("Computing risk scores & AUCs...\n")

RS_list <- list()
AUC_list <- list()
for (method in names(model)) {
  tryCatch({
    RS_list[[method]] <- CalPredictScore(model[[method]], all_combined)
    AUC_list[[method]] <- RunEval(model[[method]], Test_set, Test_class,
                                  Train_name = train_cohort,
                                  Train_set = Train_set, Train_label = Train_class,
                                  cohortVar = "Cohort", classVar = classVar)
  }, error = function(e) {})
}

AUC_mat <- do.call(rbind, AUC_list)
write.table(data.frame(ID = rownames(AUC_mat), AUC_mat), "AUC_mat.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

avg_AUC <- sort(apply(AUC_mat, 1, mean), decreasing = TRUE)
AUC_hm  <- AUC_mat[names(avg_AUC), , drop = FALSE]
avg_vals <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

# Features
fea_list <- lapply(model, ExtractVar)
fea_df <- do.call(rbind, lapply(names(fea_list), function(m) {
  if (length(fea_list[[m]]) > 0)
    data.frame(features = fea_list[[m]], algorithm = m, stringsAsFactors = FALSE)
}))
write.table(fea_df, "fea_df.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ===========================================================================
# ---- 1.3.6 Publication Figures ----
# ===========================================================================

# Color palette (Nature/Science style)
pal_train  <- "#D62728"  # red for training
pal_val    <- "#1F77B4"  # blue for validation
pal_cohort <- c("TCGA" = "#D62728", "DKFZ" = "#1F77B4",
                "Taylor" = "#2CA02C", "CancerMap" = "#FF7F0E")
pal_gene   <- c("DNAH8" = "#E41A1C", "TOP2A" = "#377EB8",
                "CHRDL1" = "#4DAF4A", "PGM5" = "#984EA3")
pal_heatmap <- c("#2166AC", "#F7F7F7", "#B2182B")

# ─────────────────────────────────────────────────────────────────
# Figure 1: AUC Heatmap (Using original SimpleHeatmap style)
# ─────────────────────────────────────────────────────────────────
cat("Figure 1: AUC Heatmap\n")

# Add "Training" / "Validation" to column names
col_labels <- colnames(AUC_hm)
col_roles  <- ifelse(col_labels == "TCGA", "Training", "Validation")
col_display <- paste0(col_labels, "\n[", col_roles, "]")
colnames(AUC_hm) <- col_display

CohortCol <- pal_cohort
names(CohortCol) <- col_display

cw <- 1.1; ch <- 0.38
hm1 <- SimpleHeatmap(AUC_hm, avg_vals, CohortCol, barCol = "#2166AC",
                     cellwidth = cw, cellheight = ch,
                     cluster_columns = FALSE, cluster_rows = FALSE)

pdf("Figure1_AUC_Heatmap.pdf",
    width  = cw * ncol(AUC_hm) + 7,
    height = max(ch * nrow(AUC_hm) * 0.42, 8))
draw(hm1)
dev.off()
png("Figure1_AUC_Heatmap.png",
    width  = (cw * ncol(AUC_hm) + 7) * 150,
    height = max(ch * nrow(AUC_hm) * 0.42, 8) * 150, res = 150)
draw(hm1)
dev.off()

# Restore original column names for later use
colnames(AUC_hm) <- col_labels

# ─────────────────────────────────────────────────────────────────
# Figure 2: Top Models Barplot
# ─────────────────────────────────────────────────────────────────
cat("Figure 2: Top Models Barplot\n")
top_n <- min(20, length(avg_AUC))
top_auc <- avg_AUC[1:top_n]

pdf("Figure2_TopModels.pdf", width = 13, height = 7.5)
par(mar = c(5, 24, 4, 2), family = "serif")
cols_bar <- colorRampPalette(c("#2166AC", "#D62728"))(top_n)
bp <- barplot(rev(top_auc), horiz = TRUE, las = 2, col = rev(cols_bar),
              border = NA, cex.names = 0.75, cex.axis = 1, cex.lab = 1.1,
              main = "", xlab = "Mean AUC (Training + 3 Validation Cohorts)",
              xlim = c(0, max(top_auc) * 1.14))
text(rev(top_auc) + 0.015, bp, format(rev(top_auc), digits = 3),
     pos = 4, cex = 0.7, col = "grey20")
abline(v = 0.7, lty = 2, col = "grey60", lwd = 1.2)
mtext("Top 20 Machine Learning Models for BCR Prediction",
      side = 3, line = 2.5, cex = 1.3, font = 2)
mtext("Training: TCGA-PRAD | Validation: DKFZ, Taylor, CancerMap",
      side = 3, line = 1, cex = 0.85, col = "grey40")
dev.off()

png("Figure2_TopModels.png", width = 13 * 150, height = 7.5 * 150, res = 150)
par(mar = c(5, 24, 4, 2), family = "serif")
bp <- barplot(rev(top_auc), horiz = TRUE, las = 2, col = rev(cols_bar),
              border = NA, cex.names = 0.75, cex.axis = 1, cex.lab = 1.1,
              main = "", xlab = "Mean AUC (Training + 3 Validation Cohorts)",
              xlim = c(0, max(top_auc) * 1.14))
text(rev(top_auc) + 0.015, bp, format(rev(top_auc), digits = 3),
     pos = 4, cex = 0.7, col = "grey20")
abline(v = 0.7, lty = 2, col = "grey60", lwd = 1.2)
mtext("Top 20 Machine Learning Models for BCR Prediction",
      side = 3, line = 2.5, cex = 1.3, font = 2)
mtext("Training: TCGA-PRAD | Validation: DKFZ, Taylor, CancerMap",
      side = 3, line = 1, cex = 0.85, col = "grey40")
dev.off()

# ─────────────────────────────────────────────────────────────────
# Figure 3: ROC Curves
# ─────────────────────────────────────────────────────────────────
cat("Figure 3: ROC Curves\n")
dir.create("ROC", showWarnings = FALSE)
top5 <- names(avg_AUC)[1:min(5, length(avg_AUC))]
cols_roc <- c("#D62728", "#1F77B4", "#2CA02C", "#FF7F0E", "#9467BD")
ltys_roc <- c(1, 1, 1, 1, 2)

# A. Combined Validation Set
val_y <- as.numeric(as.character(all_label_combined$outcome[test_idx]))

pdf("Figure3A_ROC_Validation.pdf", width = 7, height = 6.5)
par(mar = c(5, 5, 3, 1), family = "serif")
roc1 <- roc(val_y, RS_list[[top5[1]]][test_idx], quiet = TRUE)
plot(roc1, col = cols_roc[1], lwd = 3, cex.lab = 1.2, cex.axis = 1,
     main = "", legacy.axes = TRUE)
leg_txt <- sprintf("%s (AUC=%.3f)", top5,
                   sapply(top5, function(m) auc(roc(val_y, RS_list[[m]][test_idx], quiet = TRUE))))
for (i in 2:length(top5)) {
  roc_i <- roc(val_y, RS_list[[top5[i]]][test_idx], quiet = TRUE)
  lines(roc_i, col = cols_roc[i], lwd = 3)
}
legend("bottomright", legend = leg_txt, col = cols_roc, lwd = 3,
       cex = 0.8, bty = "n")
mtext("ROC Curves — Validation Cohorts", side = 3, line = 1.5, cex = 1.3, font = 2)
mtext("DKFZ + Taylor + CancerMap (n = 478, BCR = 124)", side = 3, line = 0.2, cex = 0.9, col = "grey40")
dev.off()

png("Figure3A_ROC_Validation.png", width = 7 * 150, height = 6.5 * 150, res = 150)
par(mar = c(5, 5, 3, 1), family = "serif")
roc1 <- roc(val_y, RS_list[[top5[1]]][test_idx], quiet = TRUE)
plot(roc1, col = cols_roc[1], lwd = 3, cex.lab = 1.2, cex.axis = 1,
     main = "", legacy.axes = TRUE)
for (i in 2:length(top5)) {
  roc_i <- roc(val_y, RS_list[[top5[i]]][test_idx], quiet = TRUE)
  lines(roc_i, col = cols_roc[i], lwd = 3)
}
legend("bottomright", legend = leg_txt, col = cols_roc, lwd = 3, cex = 0.8, bty = "n")
mtext("ROC Curves — Validation Cohorts", side = 3, line = 1.5, cex = 1.3, font = 2)
mtext("DKFZ + Taylor + CancerMap (n = 478, BCR = 124)", side = 3, line = 0.2, cex = 0.9, col = "grey40")
dev.off()

# B. Per-cohort ROC
for (cn in all_cohorts) {
  idx <- which(all_label_combined$Cohort == cn)
  y_true <- as.numeric(as.character(all_label_combined$outcome[idx]))
  role <- ifelse(cn == "TCGA", "Training", "Validation")

  pdf(paste0("ROC/Figure3_ROC_", cn, ".pdf"), width = 6.5, height = 6)
  par(mar = c(5, 5, 3, 1), family = "serif")
  roc1 <- roc(y_true, RS_list[[top5[1]]][idx], quiet = TRUE)
  plot(roc1, col = cols_roc[1], lwd = 3, cex.lab = 1.2, cex.axis = 1,
       main = "", legacy.axes = TRUE)
  for (i in 2:length(top5)) {
    lines(roc(y_true, RS_list[[top5[i]]][idx], quiet = TRUE),
          col = cols_roc[i], lwd = 3)
  }
  legend("bottomright", legend = leg_txt, col = cols_roc, lwd = 3,
         cex = 0.7, bty = "n")
  mtext(paste0("ROC Curve — ", cn, " (", role, ", n=", length(idx), ")"),
        side = 3, line = 1.2, cex = 1.2, font = 2)
  dev.off()

  png(paste0("ROC/Figure3_ROC_", cn, ".png"), width = 6.5 * 150, height = 6 * 150, res = 150)
  par(mar = c(5, 5, 3, 1), family = "serif")
  roc1 <- roc(y_true, RS_list[[top5[1]]][idx], quiet = TRUE)
  plot(roc1, col = cols_roc[1], lwd = 3, cex.lab = 1.2, cex.axis = 1,
       main = "", legacy.axes = TRUE)
  for (i in 2:length(top5)) {
    lines(roc(y_true, RS_list[[top5[i]]][idx], quiet = TRUE),
          col = cols_roc[i], lwd = 3)
  }
  legend("bottomright", legend = leg_txt, col = cols_roc, lwd = 3,
         cex = 0.7, bty = "n")
  mtext(paste0("ROC Curve — ", cn, " (", role, ", n=", length(idx), ")"),
        side = 3, line = 1.2, cex = 1.2, font = 2)
  dev.off()
}

# ─────────────────────────────────────────────────────────────────
# Figure 4: Gene Expression Boxplots
# ─────────────────────────────────────────────────────────────────
cat("Figure 4: Gene Expression Boxplots\n")
dir.create("Boxplots", showWarnings = FALSE)

for (gene_name in colnames(all_expr)) {
  pdf(paste0("Boxplots/Figure4_Boxplot_", gene_name, ".pdf"), width = 11, height = 9.5)
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2), family = "serif")
  for (cn in all_cohorts) {
    role <- ifelse(cn == "TCGA", "[Training]", "[Validation]")
    idx <- which(all_label$Cohort == cn)
    vals <- all_expr[idx, gene_name]
    labs <- factor(all_label$outcome[idx], levels = c("0", "1"))
    boxplot(vals ~ labs, main = paste(gene_name, role, cn), las = 1,
            col = c("#6BAED6", "#FB6A4A"), border = c("#2171B5", "#CB181D"),
            names = c("non-BCR", "BCR"), ylab = "log₂ Expression",
            cex.lab = 1.1, cex.axis = 1, cex.main = 1.2, lwd = 1.5)
    tryCatch({
      tt <- t.test(vals[labs == "1"], vals[labs == "0"])
      mtext(paste("t-test P =", format(tt$p.value, digits = 3, scientific = TRUE)),
            side = 3, line = -1.8, cex = 0.8, col = "grey30")
    }, error = function(e) {})
    # Add stripchart for individual data points
    stripchart(vals ~ labs, vertical = TRUE, method = "jitter",
               pch = 16, col = adjustcolor("grey30", 0.4), cex = 0.5, add = TRUE)
  }
  dev.off()

  png(paste0("Boxplots/Figure4_Boxplot_", gene_name, ".png"),
      width = 11 * 150, height = 9.5 * 150, res = 150)
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2), family = "serif")
  for (cn in all_cohorts) {
    role <- ifelse(cn == "TCGA", "[Training]", "[Validation]")
    idx <- which(all_label$Cohort == cn)
    vals <- all_expr[idx, gene_name]
    labs <- factor(all_label$outcome[idx], levels = c("0", "1"))
    boxplot(vals ~ labs, main = paste(gene_name, role, cn), las = 1,
            col = c("#6BAED6", "#FB6A4A"), border = c("#2171B5", "#CB181D"),
            names = c("non-BCR", "BCR"), ylab = "log₂ Expression",
            cex.lab = 1.1, cex.axis = 1, cex.main = 1.2, lwd = 1.5)
    tryCatch({
      tt <- t.test(vals[labs == "1"], vals[labs == "0"])
      mtext(paste("t-test P =", format(tt$p.value, digits = 3, scientific = TRUE)),
            side = 3, line = -1.8, cex = 0.8, col = "grey30")
    }, error = function(e) {})
    stripchart(vals ~ labs, vertical = TRUE, method = "jitter",
               pch = 16, col = adjustcolor("grey30", 0.4), cex = 0.5, add = TRUE)
  }
  dev.off()
}

# ─────────────────────────────────────────────────────────────────
# Figure 5: LODO Heatmap (using SimpleHeatmap)
# ─────────────────────────────────────────────────────────────────
cat("Figure 5: LODO Cross-Validation Heatmap\n")

top_lodo <- names(avg_AUC)[1:min(20, length(avg_AUC))]
lodo_auc_list <- list()

for (test_cn in all_cohorts) {
  cat("  LODO left-out:", test_cn, "\n")
  tr_idx <- which(all_label$Cohort != test_cn)
  te_idx <- which(all_label$Cohort == test_cn)
  TrX <- scaleData(all_expr[tr_idx, , drop = FALSE], centerFlags = TRUE, scaleFlags = TRUE)
  TeX <- scaleData(all_expr[te_idx, , drop = FALSE],
                   cohort = all_label$Cohort[te_idx], centerFlags = TRUE, scaleFlags = TRUE)
  TrY <- all_label[tr_idx, , drop = FALSE]
  TeY <- all_label[te_idx, , drop = FALSE]

  for (mn in top_lodo) {
    tryCatch({
      sp <- strsplit(mn, "\\+")[[1]]
      fm <- if (length(sp) == 1) sp[1] else sp[2]
      fit <- RunML(fm, TrX, TrY, mode = "Model", classVar = classVar)
      aucs <- RunEval(fit, TeX, TeY, cohortVar = "Cohort", classVar = classVar)
      lodo_auc_list[[paste(mn, test_cn, sep = "@")]] <- aucs
    }, error = function(e) {})
  }
}

keys <- names(lodo_auc_list)
lodo_methods <- unique(sapply(strsplit(keys, "@"), `[`, 1))
lodo_cohorts <- unique(sapply(strsplit(keys, "@"), `[`, 2))
lodo_mat <- matrix(NA, length(lodo_methods), length(lodo_cohorts),
                   dimnames = list(lodo_methods, lodo_cohorts))
for (k in keys) { p <- strsplit(k, "@")[[1]]; lodo_mat[p[1], p[2]] <- mean(lodo_auc_list[[k]], na.rm = TRUE) }

write.table(data.frame(ID = rownames(lodo_mat), lodo_mat), "LODO_AUC_mat.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

avg_lodo <- sort(apply(lodo_mat, 1, mean, na.rm = TRUE), decreasing = TRUE)
lodo_hm  <- lodo_mat[names(avg_lodo), , drop = FALSE]
avg_lv   <- as.numeric(format(avg_lodo, digits = 3, nsmall = 3))

# LODO display labels
lodo_labels <- paste0(colnames(lodo_hm), "\n[Left-Out]")
colnames(lodo_hm) <- lodo_labels
CohortCol_lodo <- pal_cohort
names(CohortCol_lodo) <- lodo_labels

hm_lodo <- SimpleHeatmap(lodo_hm, avg_lv, CohortCol_lodo, barCol = "#8B0000",
                         cellwidth = 1.2, cellheight = 0.5,
                         cluster_columns = FALSE, cluster_rows = FALSE)

pdf("Figure5_LODO_Heatmap.pdf",
    width  = 1.2 * ncol(lodo_hm) + 7,
    height = max(0.5 * nrow(lodo_hm) * 0.42, 5))
draw(hm_lodo)
dev.off()
png("Figure5_LODO_Heatmap.png",
    width  = (1.2 * ncol(lodo_hm) + 7) * 150,
    height = max(0.5 * nrow(lodo_hm) * 0.42, 5) * 150, res = 150)
draw(hm_lodo)
dev.off()

# Restore
colnames(lodo_hm) <- colnames(lodo_mat)

# ─────────────────────────────────────────────────────────────────
# Figure 6: Feature Frequency
# ─────────────────────────────────────────────────────────────────
cat("Figure 6: Gene Retention Frequency\n")
fea_freq <- sort(table(fea_df$features), decreasing = TRUE)
freq_pct <- round(fea_freq / length(model) * 100, 1)

pdf("Figure6_FeatureFrequency.pdf", width = 6, height = 5.5)
par(mar = c(5, 5, 4, 2), family = "serif")
bp <- barplot(fea_freq, col = pal_gene[names(fea_freq)], border = NA,
              las = 2, cex.names = 1.3, cex.axis = 1.1, cex.lab = 1.1,
              ylim = c(0, max(fea_freq) * 1.18),
              ylab = "Number of Models Retaining Gene",
              main = "")
text(bp, fea_freq + max(fea_freq) * 0.03,
     labels = paste0(fea_freq, " (", freq_pct, "%)"), cex = 1, font = 2)
mtext("Gene Retention Frequency Across 95 ML Models",
      side = 3, line = 2, cex = 1.2, font = 2)
dev.off()

png("Figure6_FeatureFrequency.png", width = 6 * 150, height = 5.5 * 150, res = 150)
par(mar = c(5, 5, 4, 2), family = "serif")
bp <- barplot(fea_freq, col = pal_gene[names(fea_freq)], border = NA,
              las = 2, cex.names = 1.3, cex.axis = 1.1, cex.lab = 1.1,
              ylim = c(0, max(fea_freq) * 1.18),
              ylab = "Number of Models Retaining Gene",
              main = "")
text(bp, fea_freq + max(fea_freq) * 0.03,
     labels = paste0(fea_freq, " (", freq_pct, "%)"), cex = 1, font = 2)
mtext("Gene Retention Frequency Across 95 ML Models",
      side = 3, line = 2, cex = 1.2, font = 2)
dev.off()

# ─────────────────────────────────────────────────────────────────
# Figure 7: Best Model Variable Importance
# ─────────────────────────────────────────────────────────────────
cat("Figure 7: Best Model Variable Importance\n")
best_name <- names(avg_AUC)[1]
best_fit  <- model[[best_name]]
best_feat <- fea_list[[best_name]]
fc <- class(best_fit)[1]

pdf("Figure7_VariableImportance.pdf", width = 7, height = 5.5)
par(mar = c(6, 5, 4, 2), family = "serif")
if (fc == "rfsrc") {
  imp <- sort(best_fit$importance, decreasing = TRUE)
  barplot(imp, las = 2, col = pal_gene[names(imp)], border = NA,
          cex.names = 1.2, cex.axis = 1.1, cex.lab = 1.1,
          ylab = "Random Forest Variable Importance (VIMP)", main = "")
  mtext(paste0("Best Model: ", best_name, "  |  Mean Validation AUC = ",
               format(avg_AUC[1], digits = 3)),
        side = 3, line = 2, cex = 1.1, font = 2, col = "grey20")
} else if (fc %in% c("lognet", "glm")) {
  cf <- coef(best_fit)
  if (fc == "lognet") { cf <- as.numeric(cf[, 1]); names(cf) <- rownames(coef(best_fit)) }
  cf <- cf[!names(cf) %in% "(Intercept)"]
  cols_cf <- ifelse(cf > 0, "#FB6A4A", "#6BAED6")
  barplot(cf, las = 2, col = cols_cf, border = NA,
          cex.names = 1.2, cex.axis = 1.1, cex.lab = 1.1,
          ylab = "Coefficient", main = "")
  mtext(paste0("Best Model: ", best_name, "  |  Mean Validation AUC = ",
               format(avg_AUC[1], digits = 3)),
        side = 3, line = 2, cex = 1.1, font = 2, col = "grey20")
} else {
  plot.new()
  text(0.5, 0.5, paste0("Best Model: ", best_name, "\nFeatures: ",
                         paste(best_feat, collapse = ", "),
                         "\nMean Validation AUC = ", round(avg_AUC[1], 3)),
       cex = 1.5)
}
dev.off()

png("Figure7_VariableImportance.png", width = 7 * 150, height = 5.5 * 150, res = 150)
par(mar = c(6, 5, 4, 2), family = "serif")
if (fc == "rfsrc") {
  imp <- sort(best_fit$importance, decreasing = TRUE)
  barplot(imp, las = 2, col = pal_gene[names(imp)], border = NA,
          cex.names = 1.2, cex.axis = 1.1, cex.lab = 1.1,
          ylab = "Random Forest Variable Importance (VIMP)", main = "")
  mtext(paste0("Best Model: ", best_name, "  |  Mean Validation AUC = ",
               format(avg_AUC[1], digits = 3)),
        side = 3, line = 2, cex = 1.1, font = 2, col = "grey20")
} else if (fc %in% c("lognet", "glm")) {
  cf <- coef(best_fit)
  if (fc == "lognet") { cf <- as.numeric(cf[, 1]); names(cf) <- rownames(coef(best_fit)) }
  cf <- cf[!names(cf) %in% "(Intercept)"]
  cols_cf <- ifelse(cf > 0, "#FB6A4A", "#6BAED6")
  barplot(cf, las = 2, col = cols_cf, border = NA,
          cex.names = 1.2, cex.axis = 1.1, cex.lab = 1.1,
          ylab = "Coefficient", main = "")
  mtext(paste0("Best Model: ", best_name, "  |  Mean Validation AUC = ",
               format(avg_AUC[1], digits = 3)),
        side = 3, line = 2, cex = 1.1, font = 2, col = "grey20")
} else {
  plot.new()
  text(0.5, 0.5, paste0("Best Model: ", best_name, "\nFeatures: ",
                         paste(best_feat, collapse = ", "),
                         "\nMean Validation AUC = ", round(avg_AUC[1], 3)),
       cex = 1.5)
}
dev.off()

# ===========================================================================
# ---- 1.3.7 Save Summary ----
# ===========================================================================
saveRDS(list(
  avg_AUC = avg_AUC, AUC_mat = AUC_mat, top_models = top5,
  best_model = best_name, best_features = best_feat,
  lodo_auc = lodo_mat, fea_freq = fea_freq,
  core_genes = core_ensembl, model_count = length(model)
), "summary_results.RDS")

# ===========================================================================
# ---- 1.3.8 Print Final Results (for paper) ----
# ===========================================================================
cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("               FINAL RESULTS FOR MANUSCRIPT\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("--- Study Cohorts ---\n")
cat(sprintf("  Training:   TCGA-PRAD     n=%d  BCR=%d\n",
            sum(all_label$Cohort == "TCGA"),
            sum(all_label$Cohort == "TCGA" & all_label$outcome == "1")))
for(cn in test_cohorts) {
  cat(sprintf("  Validation: %-12s n=%d  BCR=%d\n", cn,
              sum(all_label$Cohort == cn),
              sum(all_label$Cohort == cn & all_label$outcome == "1")))
}

cat(sprintf("\n--- ML Models: %d trained from 107 candidate algorithms ---\n", length(model)))

cat("\n--- Top 5 Models ---\n")
for(i in 1:min(5, length(avg_AUC))) {
  m <- names(avg_AUC)[i]
  cat(sprintf("  %d. %-30s  Mean AUC = %.3f\n", i, m, avg_AUC[i]))
}

cat(sprintf("\n--- Best Model: %s ---\n", best_name))
cat(sprintf("  Features: %s\n", paste(best_feat, collapse = ", ")))
cat("  Per-cohort AUC:\n")
print(round(AUC_mat[best_name, , drop = FALSE], 3))

cat("\n--- LODO Top 3 ---\n")
for(i in 1:min(3, length(avg_lodo))) {
  m <- names(avg_lodo)[i]
  cat(sprintf("  %d. %-30s  Mean AUC = %.3f\n", i, m, avg_lodo[i]))
}

cat("\n--- Gene Retention ---\n")
for(g in names(fea_freq)) {
  cat(sprintf("  %s: %d/%d models (%.1f%%)\n",
              g, fea_freq[g], length(model), 100*fea_freq[g]/length(model)))
}

cat("\n═══════════════════════════════════════════════════════════\n")
cat("Historical outputs:", legacy_output(), "\n")
cat("═══════════════════════════════════════════════════════════\n")
