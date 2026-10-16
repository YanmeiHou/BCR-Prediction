# ==== 2.3 单基因ROC ====
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
# Single-Gene BCR ROC Curves
#
# For each of the 4 genes (DNAH8, TOP2A, CHRDL1, PGM5):
#   One figure with 4 ROC curves (TCGA, DKFZ, Taylor, CancerMap)
#
# Output: SingleGene_ROC/
################################################################################

library(ggplot2)
library(pROC)
library(dplyr)

setwd(legacy_output())
dir.create("SingleGene_ROC", showWarnings = FALSE)

# Load prepared data
pd <- readRDS("prepared_data.RDS")
expr <- pd$expr
label <- pd$label

core_genes <- colnames(expr)  # DNAH8, TOP2A, CHRDL1, PGM5
cohorts <- c("TCGA", "DKFZ", "Taylor", "CancerMap")

# Cohort colors (consistent with publication)
pal_cohort <- c("TCGA" = "#D62728", "DKFZ" = "#1F77B4",
                "Taylor" = "#2CA02C", "CancerMap" = "#FF7F0E")

# Use sans = Arial on Windows
if (.Platform$OS.type == "windows") windowsFonts(sans = windowsFont("Arial"))

# ═══════════════════════════════════════════
# Generate 4 figures (one per gene)
# ═══════════════════════════════════════════
for (gene in core_genes) {
  cat(sprintf("─── %s ───\n", gene))

  # Open PDF
  pdf(sprintf("SingleGene_ROC/ROC_%s.pdf", gene),
      width = 7, height = 6.5, family = "sans")

  # Plot settings
  par(mar = c(5, 5, 3, 1), family = "sans")

  first_cohort <- TRUE
  auc_values <- c()

  for (cn in cohorts) {
    idx <- which(label$Cohort == cn)
    y_true <- as.numeric(as.character(label$outcome[idx]))
    x_pred <- expr[idx, gene]

    # Compute ROC
    roc_obj <- roc(y_true, x_pred, quiet = TRUE)
    auc_val <- auc(roc_obj)
    auc_values[cn] <- auc_val

    # Plot
    if (first_cohort) {
      plot(roc_obj, col = pal_cohort[cn], lwd = 3, cex.lab = 1.2, cex.axis = 1.1,
           main = "", legacy.axes = TRUE, family = "sans")
      first_cohort <- FALSE
    } else {
      lines(roc_obj, col = pal_cohort[cn], lwd = 3)
    }
  }

  # Legend
  leg_txt <- sprintf("%s (AUC=%.3f, n=%d)", names(auc_values), auc_values,
                     table(label$Cohort)[names(auc_values)])
  legend("bottomright", legend = leg_txt, col = pal_cohort[names(auc_values)],
         lwd = 3, cex = 0.9, bty = "n", text.font = 1)

  # Title
  title(main = paste(gene, "— BCR Prediction ROC"),
        cex.main = 1.4, font.main = 2, family = "sans")
  mtext("TCGA (Training) + DKFZ, Taylor, CancerMap (Validation)",
        side = 3, line = 0.2, cex = 0.9, col = "grey40", family = "sans")

  dev.off()

  # Also save PNG
  png(sprintf("SingleGene_ROC/ROC_%s.png", gene),
      width = 7 * 150, height = 6.5 * 150, res = 150)
  par(mar = c(5, 5, 3, 1), family = "sans")

  first_cohort <- TRUE
  for (cn in cohorts) {
    idx <- which(label$Cohort == cn)
    y_true <- as.numeric(as.character(label$outcome[idx]))
    x_pred <- expr[idx, gene]
    roc_obj <- roc(y_true, x_pred, quiet = TRUE)

    if (first_cohort) {
      plot(roc_obj, col = pal_cohort[cn], lwd = 3, cex.lab = 1.2, cex.axis = 1.1,
           main = "", legacy.axes = TRUE, family = "sans")
      first_cohort <- FALSE
    } else {
      lines(roc_obj, col = pal_cohort[cn], lwd = 3)
    }
  }

  legend("bottomright", legend = leg_txt, col = pal_cohort[names(auc_values)],
         lwd = 3, cex = 0.9, bty = "n", text.font = 1)
  title(main = paste(gene, "— BCR Prediction ROC"),
        cex.main = 1.4, font.main = 2, family = "sans")
  mtext("TCGA (Training) + DKFZ, Taylor, CancerMap (Validation)",
        side = 3, line = 0.2, cex = 0.9, col = "grey40", family = "sans")

  dev.off()

  cat(sprintf("  AUC: %s\n", paste(sprintf("%s=%.3f", names(auc_values), auc_values), collapse = ", ")))
}

cat("\n═══ 4 ROC figures saved to SingleGene_ROC/ ═══\n")
