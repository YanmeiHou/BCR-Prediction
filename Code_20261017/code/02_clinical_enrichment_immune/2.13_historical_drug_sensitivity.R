# 2.13 Historical oncoPredict implementation; manuscript linkage unconfirmed.
if(Sys.getenv("PCA_ENABLE_HISTORICAL")!="1") stop("Read docs/HISTORICAL_SOURCE.md before historical replay")
###############################################################################
# oncoPredict 前列腺癌四队列独立药物敏感性预测 + Meta 分析 (Fixed)
#
# 策略：各队列独立运行 oncoPredict (GDSC2) → BCR vs Non-BCR 比较 → Meta 汇总
#
# 修复要点：
#   0. 加入诊断步骤，检查 eSet 结构、BCR 编码、基因 ID 格式
#   1. 检查训练集与测试集基因重叠
#   2. robase BCR 状态解析（支持 "0"/"1", 0/1, "Yes"/"No", TRUE/FALSE）
#   3. calcPhenotype 多策略回退（eb → standard → none）
#   4. 详细的错误诊断输出
#   5. 最小样本量自适应
###############################################################################

library(oncoPredict)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(officer)
library(metafor)
library(reshape2)
library(grid)
library(Biobase)
library(ggrepel)

# ---- Arial font ----
theme_arial <- theme_minimal(base_family = "sans") +
  theme(text = element_text(family = "sans", color = "black"),
    plot.title   = element_text(family = "sans", color = "black", face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(family = "sans", color = "black", hjust = 0.5, size = 9),
    axis.title   = element_text(family = "sans", color = "black", size = 11),
    axis.text    = element_text(family = "sans", color = "black", size = 9),
    legend.text  = element_text(family = "sans", color = "black"),
    legend.title = element_text(family = "sans", color = "black"),
    strip.text   = element_text(family = "sans", color = "black", face = "bold", size = 9))

# ---- Paths ----
root       <- project_path('data','external','bulk')
data_dir   <- root
drug_dir   <- project_path('data','external','drug')
train_dir  <- file.path(drug_dir, "DataFiles", "Training Data")
out_dir    <- project_path('outputs','historical_drug_sensitivity')
fig_dir    <- file.path(out_dir, "Figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Log file (clean up any leftover sinks from interrupted runs) ----
while (sink.number(type = "message") > 0) sink(type = "message")
while (sink.number() > 0) sink()
log_file <- file.path(out_dir, "run_log.txt")
log_con  <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, split = TRUE, type = "message")

cat("========================================\n")
cat("oncoPredict: Prostate Cancer Drug Sensitivity (BCR)\n")
cat("Run time:", as.character(Sys.time()), "\n")
cat("========================================\n\n")

# 4 prostate cancer cohorts with BCR data
cohorts       <- c("DKFZ", "TCGA", "Taylor", "CancerMap")
cohort_esets  <- c("DKFZ"     = "DKFZ_eSet.RDS",
                   "TCGA"     = "TCGA-PRAD_eSet.RDS",
                   "Taylor"   = "Taylor_eSet.RDS",
                   "CancerMap" = "CancerMap_eSet.RDS")

# ============================================================
# 0. DIAGNOSTICS: Inspect eSet objects before anything else
# ============================================================
cat("=== STEP 0: DIAGNOSTICS ===\n\n")

for (coh in cohorts) {
  eset_path <- file.path(data_dir, cohort_esets[[coh]])
  if (!file.exists(eset_path)) {
    cat(sprintf("  [ERROR] %s: file not found: %s\n", coh, eset_path))
    next
  }
  eset <- readRDS(eset_path)
  pd   <- pData(eset)
  cat(sprintf("  --- %s ---\n", coh))
  cat(sprintf("    Class of eSet: %s\n", paste(class(eset), collapse = ", ")))
  cat(sprintf("    dim(exprs): %d genes x %d samples\n", nrow(exprs(eset)), ncol(exprs(eset))))
  cat(sprintf("    colnames(pData): %s\n", paste(colnames(pd), collapse = ", ")))
  
  # Inspect BCR column
  bcr_cols <- grep("bcr|BCR|recur|prog", colnames(pd), value = TRUE, ignore.case = TRUE)
  if (length(bcr_cols) == 0) {
    cat("    [WARNING] No BCR-like column found in pData! Printing all columns:\n")
    cat(sprintf("      %s\n", paste(colnames(pd), collapse = ", ")))
  } else {
    for (bc in bcr_cols) {
      vals <- pd[[bc]]
      cat(sprintf("    Column '%s': class=%s, unique values (first 10): %s\n",
                  bc, class(vals),
                  paste(head(unique(vals), 10), collapse = ", ")))
      cat(sprintf("      NA count: %d / %d\n", sum(is.na(vals)), length(vals)))
    }
  }
  
  # Inspect gene ID format (first 10 rownames)
  rn <- head(rownames(exprs(eset)), 10)
  cat(sprintf("    First 10 gene IDs: %s\n", paste(rn, collapse = ", ")))
  cat("\n")
}

# ============================================================
# 0b. Load GDSC2 training data & inspect gene IDs
# ============================================================
cat("--- GDSC2 Training Data ---\n")
traDataExp <- readRDS(file.path(train_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
traDataRes <- readRDS(file.path(train_dir, "GDSC2_Res.rds"))
cat(sprintf("  Training genes: %d, cell lines: %d, drugs: %d\n",
            nrow(traDataExp), ncol(traDataExp), ncol(traDataRes)))
cat(sprintf("  First 10 training gene IDs: %s\n",
            paste(head(rownames(traDataExp), 10), collapse = ", ")))
cat("\n")

# ============================================================
# 1. Load expression & BCR status (robust parsing)
# ============================================================
cat("=== STEP 1: Loading expression & BCR data ===\n\n")

expr_list   <- list()
bcr_list    <- list()

# Helper: robust BCR parsing
parse_bcr <- function(pd) {
  # Find BCR column
  bcr_col <- grep("^bcr_status$|^BCR$|^bcr$", colnames(pd), value = TRUE, ignore.case = TRUE)
  if (length(bcr_col) == 0) {
    # Try fuzzy match
    bcr_col <- grep("bcr|recurrence|progression", colnames(pd),
                    value = TRUE, ignore.case = TRUE)
  }
  if (length(bcr_col) == 0) {
    stop("No BCR column found in pData. Available columns: ",
         paste(colnames(pd), collapse = ", "))
  }
  bcr_col <- bcr_col[1]
  cat(sprintf("    Using BCR column: '%s'\n", bcr_col))
  
  raw <- pd[[bcr_col]]
  cat(sprintf("    Raw unique values: %s\n",
              paste(unique(raw), collapse = ", ")))
  
  # Parse to 0/1 numeric
  if (is.numeric(raw)) {
    # Numeric: 0/1
    bcr_num <- as.numeric(raw)
  } else if (is.logical(raw)) {
    bcr_num <- as.numeric(raw)
  } else if (is.character(raw)) {
    raw_trim <- trimws(raw)
    # "0"/"1"
    if (all(raw_trim %in% c("0", "1", NA_character_))) {
      bcr_num <- as.numeric(raw_trim)
    } else if (all(toupper(raw_trim) %in% c("YES", "NO", NA_character_))) {
      bcr_num <- ifelse(toupper(raw_trim) == "YES", 1, 0)
    } else if (all(toupper(raw_trim) %in% c("TRUE", "FALSE", NA_character_))) {
      bcr_num <- ifelse(toupper(raw_trim) == "TRUE", 1, 0)
    } else {
      # Try to map: assume 1/positive is "BCR" or "Yes" or "Recurrence" etc.
      bcr_num <- ifelse(tolower(raw_trim) %in% c("1", "yes", "true", "bcr", "recurrence",
                                                  "recurrent", "progressed", "positive",
                                                  "pos", "event"), 1,
                 ifelse(tolower(raw_trim) %in% c("0", "no", "false", "non-bcr", "disease-free",
                                                  "free", "negative", "neg", "nonevent"), 0, NA))
    }
  } else if (is.factor(raw)) {
    # Factor: convert to character then parse
    raw_char <- as.character(raw)
    return(parse_bcr(data.frame(bcr_status = raw_char, stringsAsFactors = FALSE)))
  } else {
    stop("Unsupported BCR column type: ", class(raw))
  }
  
  cat(sprintf("    Parsed: BCR=1 (n=%d), Non-BCR=0 (n=%d), NA (n=%d)\n",
              sum(bcr_num == 1, na.rm = TRUE),
              sum(bcr_num == 0, na.rm = TRUE),
              sum(is.na(bcr_num))))
  return(bcr_num)
}

for (coh in cohorts) {
  eset_path <- file.path(data_dir, cohort_esets[[coh]])
  if (!file.exists(eset_path)) {
    cat(sprintf("  [SKIP] %s: file not found: %s\n", coh, eset_path))
    next
  }
  eset <- readRDS(eset_path)
  pd   <- pData(eset)
  
  # Get expression matrix
  expr_mat <- exprs(eset)
  cat(sprintf("  %s raw: %d genes x %d samples\n", coh, nrow(expr_mat), ncol(expr_mat)))
  
  # Log2 transform if needed
  if (max(expr_mat, na.rm = TRUE) > 100) {
    cat(sprintf("    Applying log2(x+1) transformation (max=%.1f)\n", max(expr_mat, na.rm = TRUE)))
    expr_mat <- log2(expr_mat + 1)
  }
  
  # Parse BCR status
  bcr_num <- tryCatch(parse_bcr(pd), error = function(e) {
    cat(sprintf("    [ERROR] BCR parsing failed: %s\n", e$message))
    return(NULL)
  })
  if (is.null(bcr_num)) next
  
  # Filter to samples with known BCR status
  keep <- which(!is.na(bcr_num))
  if (length(keep) == 0) {
    cat(sprintf("    [SKIP] %s: no samples with known BCR status\n", coh))
    next
  }
  expr_mat <- expr_mat[, keep, drop = FALSE]
  bcr_num  <- bcr_num[keep]
  
  bcr_vec <- ifelse(bcr_num == 1, "BCR", "Non-BCR")
  names(bcr_vec) <- colnames(expr_mat)
  
  expr_list[[coh]]   <- expr_mat
  bcr_list[[coh]]    <- bcr_vec
  
  cat(sprintf("    Final: %d genes x %d samples (BCR=%d, Non-BCR=%d)\n\n",
              nrow(expr_mat), ncol(expr_mat),
              sum(bcr_vec == "BCR"), sum(bcr_vec == "Non-BCR")))
}

cat(sprintf("Total cohorts loaded: %d / %d\n\n", length(expr_list), length(cohorts)))

if (length(expr_list) == 0) {
  stop("No cohorts loaded successfully. Check BCR column encoding in eSet objects.")
}

# ============================================================
# 2. Run oncoPredict per cohort (with fallback strategies)
# ============================================================
cat("=== STEP 2: Running oncoPredict per cohort (GDSC2) ===\n\n")

# Checkpoint directory
all_drug_results <- list()
checkpoint_dir <- file.path(out_dir, "Checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

# Resume from checkpoints
for (coh in names(expr_list)) {
  ckpt_file <- file.path(checkpoint_dir, paste0("oncoPredict_", coh, ".rds"))
  if (file.exists(ckpt_file)) {
    all_drug_results[[coh]] <- readRDS(ckpt_file)
    cat(sprintf("  [RESUME] %s loaded from checkpoint (%d samples)\n",
                coh, nrow(all_drug_results[[coh]])))
  }
}

# ---- calcPhenotype wrapper with fallback ----
run_calcPhenotype <- function(traDataExp, traDataRes, expr_mat, coh_name, tmp_wd) {
  
  # Pre-check: gene overlap
  train_genes <- rownames(traDataExp)
  test_genes  <- rownames(expr_mat)
  common <- intersect(train_genes, test_genes)
  cat(sprintf("    Gene overlap: %d / %d training genes, %d / %d test genes\n",
              length(common), length(train_genes),
              length(common), length(test_genes)))
  
  if (length(common) < 50) {
    cat("    [CRITICAL] Very few overlapping genes! Check gene ID format.\n")
    cat(sprintf("      Training gene examples: %s\n",
                paste(head(train_genes, 5), collapse = ", ")))
    cat(sprintf("      Test gene examples: %s\n",
                paste(head(test_genes, 5), collapse = ", ")))
    return(NULL)
  }
  
  # Determine minNumSamples
  n_samples <- ncol(expr_mat)
  min_samp  <- min(10, max(3, floor(n_samples * 0.3)))
  cat(sprintf("    n_samples=%d, minNumSamples set to %d\n", n_samples, min_samp))
  
  # Strategy 1: eb (default batch correction)
  strategies <- list(
    list(batchCorrect = "eb",     label = "eb (ComBat)"),
    list(batchCorrect = "standard", label = "standard"),
    list(batchCorrect = "none",   label = "no batch correction")
  )
  
  for (strat in strategies) {
    cat(sprintf("    Trying strategy: %s...\n", strat$label))
    
    result <- tryCatch({
      pred_mat <- calcPhenotype(
        trainingExprData       = traDataExp,
        trainingPtype          = traDataRes,
        testExprData           = expr_mat,
        batchCorrect           = strat$batchCorrect,
        powerTransformPhenotype = TRUE,
        removeLowVaryingGenes  = 0.2,
        minNumSamples          = min_samp,
        printOutput            = TRUE,
        removeLowVaringGenesFrom = "rawData"
      )
      pred_mat
    }, error = function(e) {
      cat(sprintf("      [FAIL] %s: %s\n", strat$label, e$message))
      NULL
    })
    
    if (!is.null(result)) {
      cat(sprintf("    [SUCCESS] with strategy: %s\n", strat$label))
      return(result)
    }
  }
  
  cat(sprintf("    [ALL STRATEGIES FAILED] for cohort %s\n", coh_name))
  return(NULL)
}

# ---- Process each cohort ----
for (coh in names(expr_list)) {
  ckpt_file <- file.path(checkpoint_dir, paste0("oncoPredict_", coh, ".rds"))
  if (!is.null(all_drug_results[[coh]])) {
    cat(sprintf("\n  === %s === [SKIP - already completed]\n", coh))
    next
  }
  
  cat(sprintf("\n  === %s ===\n", coh))
  t_start <- Sys.time()
  
  expr_sub <- expr_list[[coh]]
  bcr_sub  <- bcr_list[[coh]]
  
  # Remove genes with NA/NaN/Inf
  keep_genes <- complete.cases(expr_sub)
  if (sum(!keep_genes) > 0) {
    cat(sprintf("    Removed %d genes with NA values\n", sum(!keep_genes)))
    expr_sub <- expr_sub[keep_genes, , drop = FALSE]
  }
  
  expr_mat <- as.matrix(expr_sub)
  
  # Create temp working directory
  tmp_wd <- file.path(out_dir, paste0("tmp_", coh))
  unlink(tmp_wd, recursive = TRUE)  # clean previous attempts
  dir.create(tmp_wd, showWarnings = FALSE, recursive = TRUE)
  
  owd <- getwd()
  setwd(tmp_wd)
  cat(sprintf("    Working dir: %s\n", tmp_wd))
  
  # Run calcPhenotype with fallback
  pred_mat <- run_calcPhenotype(traDataExp, traDataRes, expr_mat, coh, tmp_wd)
  
  if (is.null(pred_mat)) {
    setwd(owd)
    cat(sprintf("    [SKIP] %s: calcPhenotype failed with all strategies\n", coh))
    # Save error log for debugging
    err_log <- file.path(checkpoint_dir, paste0("FAILED_", coh, ".txt"))
    writeLines(c(
      paste("Cohort:", coh),
      paste("Time:", as.character(Sys.time())),
      paste("Genes:", nrow(expr_mat)),
      paste("Samples:", ncol(expr_mat)),
      paste("Training genes:", nrow(traDataExp)),
      paste("Overlap genes:", length(intersect(rownames(traDataExp), rownames(expr_mat)))),
      paste("First 5 test genes:", paste(head(rownames(expr_mat), 5), collapse = ", ")),
      paste("First 5 train genes:", paste(head(rownames(traDataExp), 5), collapse = ", "))
    ), err_log)
    next
  }
  
  # ---- Parse results ----
  if (is.matrix(pred_mat) || is.data.frame(pred_mat)) {
    cat(sprintf("    Got predictions matrix: %d samples x %d drugs\n",
                nrow(pred_mat), ncol(pred_mat)))
    drug_pred <- as.data.frame(pred_mat)
    drug_pred$Sample <- rownames(drug_pred)
  } else {
    cat("    Checking for DrugPredictions.csv on disk...\n")
    pred_file <- file.path(tmp_wd, "DrugPredictions.csv")
    if (!file.exists(pred_file)) {
      all_csv <- list.files(tmp_wd, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
      if (length(all_csv) > 0) {
        cat(sprintf("    Using alternative CSV: %s\n", all_csv[1]))
        pred_file <- all_csv[1]
      } else {
        cat("    [SKIP] No prediction output found\n")
        setwd(owd)
        next
      }
    }
    drug_pred <- read.csv(pred_file, row.names = 1, check.names = FALSE)
    drug_pred <- as.data.frame(t(drug_pred))
    drug_pred$Sample <- rownames(drug_pred)
  }
  
  # Attach BCR and cohort info
  drug_pred$BCR    <- bcr_sub[drug_pred$Sample]
  drug_pred$Cohort <- coh
  
  # Save checkpoint
  all_drug_results[[coh]] <- drug_pred
  saveRDS(drug_pred, ckpt_file)
  
  # Clean up temp dir
  setwd(owd)
  unlink(tmp_wd, recursive = TRUE)
  
  t_end <- Sys.time()
  drug_cols <- setdiff(colnames(drug_pred), c("Sample", "BCR", "Cohort"))
  cat(sprintf("    Predicted %d drugs for %d samples | Saved: %s | Time: %.1f min\n",
              length(drug_cols), nrow(drug_pred), basename(ckpt_file),
              difftime(t_end, t_start, units = "mins")))
}

saveRDS(all_drug_results, file.path(out_dir, "oncoPredict_All_Cohorts.rds"))
cat(sprintf("\n  Combined RDS saved: %d/%d cohorts\n", length(all_drug_results), length(expr_list)))

# ============================================================
# 3. BCR vs Non-BCR differential drug sensitivity per cohort
# ============================================================
cat("\n=== STEP 3: Computing BCR vs Non-BCR drug sensitivity differences ===\n\n")

if (length(all_drug_results) == 0) {
  cat("\n[FATAL] No cohorts completed successfully.\n")
  cat("Please review the diagnostic output above and check:\n")
  cat("  1. BCR column encoding in eSet pData (see STEP 0 output)\n")
  cat("  2. Gene ID overlap between training and test data (see STEP 2 output)\n")
  cat("  3. Sample sizes per cohort after BCR filtering\n")
  cat("  4. Check failed cohort logs in:", file.path(checkpoint_dir, "FAILED_*.txt"), "\n")
  stop("No cohorts completed successfully.")
}

df_all <- do.call(rbind, all_drug_results)
drug_cols <- setdiff(colnames(df_all), c("Sample", "BCR", "Cohort"))

cat(sprintf("  Combined data: %d samples, %d drugs, %d cohorts\n",
            nrow(df_all), length(drug_cols), length(unique(df_all$Cohort))))

drug_stats <- data.frame()
for (coh in unique(df_all$Cohort)) {
  sub <- subset(df_all, Cohort == coh)
  cat(sprintf("  Processing %s: %d samples (BCR=%d, Non-BCR=%d)\n",
              coh, nrow(sub), sum(sub$BCR == "BCR"), sum(sub$BCR == "Non-BCR")))
  
  for (drug in drug_cols) {
    bcr      <- sub[[drug]][sub$BCR == "BCR"]
    nonbcr   <- sub[[drug]][sub$BCR == "Non-BCR"]
    
    # Need at least 2 in each group for meaningful stats
    if (length(bcr) < 2 || length(nonbcr) < 2) next
    
    m_bcr    <- mean(bcr, na.rm = TRUE)
    m_nonbcr <- mean(nonbcr, na.rm = TRUE)
    n_bcr    <- length(bcr)
    n_nonbcr <- length(nonbcr)
    sd_bcr   <- sd(bcr, na.rm = TRUE)
    sd_nonbcr<- sd(nonbcr, na.rm = TRUE)
    diff     <- m_bcr - m_nonbcr
    
    # Hedge's g style SE (more robust)
    se_diff  <- sqrt(sd_bcr^2/n_bcr + sd_nonbcr^2/n_nonbcr)
    if (is.na(se_diff) || se_diff == 0) se_diff <- NA
    
    pv <- tryCatch(wilcox.test(bcr, nonbcr)$p.value, error = function(e) NA)
    
    drug_stats <- rbind(drug_stats, data.frame(
      Cohort = coh, Drug = drug,
      BCR_mean = m_bcr, NonBCR_mean = m_nonbcr,
      Diff = diff, SE_diff = se_diff, pvalue = pv,
      n_BCR = n_bcr, n_NonBCR = n_nonbcr,
      Direction = ifelse(diff > 0, "BCR_higher_IC50", "BCR_lower_IC50"),
      stringsAsFactors = FALSE
    ))
  }
}

drug_stats$Sig <- ifelse(drug_stats$pvalue < 0.001, "***",
                  ifelse(drug_stats$pvalue < 0.01,  "**",
                  ifelse(drug_stats$pvalue < 0.05,  "*", "ns")))

cat(sprintf("  Drug stats: %d rows (%d drugs x %d cohorts)\n",
            nrow(drug_stats), length(drug_cols), length(unique(drug_stats$Cohort))))

# ============================================================
# 4. Meta-analysis per drug
# ============================================================
cat("\n=== STEP 4: Meta-analysis per drug (Random Effects) ===\n\n")

drug_meta <- data.frame()
for (drug in drug_cols) {
  sub <- subset(drug_stats, Drug == drug)
  if (nrow(sub) < 2) next  # need at least 2 cohorts
  sub <- sub[!is.na(sub$SE_diff) & sub$SE_diff > 0, ]
  if (nrow(sub) < 2) next
  
  res <- tryCatch({
    rma(yi = sub$Diff, sei = sub$SE_diff, method = "REML")
  }, error = function(e) NULL)
  
  if (!is.null(res)) {
    pi <- tryCatch(predict(res), error = function(e) NULL)
    drug_meta <- rbind(drug_meta, data.frame(
      Drug = drug,
      k = nrow(sub),
      Pooled_Diff = res$b[1,1],
      SE_pooled = res$se,
      CI_low = res$ci.lb,
      CI_high = res$ci.ub,
      pval = res$pval,
      I2 = res$I2,
      tau2 = res$tau2,
      Q_pval = res$QEp,
      PI_low  = if(!is.null(pi)) pi$cr.lb else NA,
      PI_high = if(!is.null(pi)) pi$cr.ub else NA,
      Direction = ifelse(res$b[1,1] > 0, "BCR_higher_IC50", "BCR_lower_IC50"),
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(drug_meta) == 0) {
  stop("Meta-analysis produced no results. Check drug statistics data.")
}

drug_meta <- drug_meta[order(drug_meta$Pooled_Diff)]
drug_meta$Sig <- ifelse(drug_meta$pval < 0.001, "***",
                 ifelse(drug_meta$pval < 0.01,  "**",
                 ifelse(drug_meta$pval < 0.05,  "*", "ns")))
drug_meta$pval_FDR <- p.adjust(drug_meta$pval, method = "BH")

n_sig   <- sum(drug_meta$pval < 0.05, na.rm = TRUE)
n_fdr   <- sum(drug_meta$pval_FDR < 0.05, na.rm = TRUE)
mean_i2 <- mean(drug_meta$I2, na.rm = TRUE)

cat(sprintf("\n  %d drugs tested\n", nrow(drug_meta)))
cat(sprintf("  Significant (p<0.05): %d | FDR<0.05: %d\n", n_sig, n_fdr))
cat(sprintf("  Mean I2 = %.1f%%\n", mean_i2))

cat("\n  Top 10 differential drugs:\n")
top10 <- drug_meta[order(drug_meta$pval)]
for (i in 1:min(10, nrow(top10))) {
  cat(sprintf("    %-45s Diff=%+.4f p=%.2e I2=%.0f%% %s\n",
              top10$Drug[i], top10$Pooled_Diff[i], top10$pval[i], top10$I2[i], top10$Sig[i]))
}

# ============================================================
# 5. Generate figures
# ============================================================
cat("\n=== STEP 5: Generating figures ===\n\n")

palette_dir <- c("BCR_higher_IC50" = "#B2182B", "BCR_lower_IC50" = "#2166AC")

# ---- Fig A: Meta forest plot (top 30) ----
cat("  [1/4] Meta forest plot...\n")

top_drugs <- drug_meta[order(drug_meta$pval)]
top_drugs <- top_drugs[1:min(30, nrow(top_drugs))]
top_drugs$Display <- factor(top_drugs$Drug, levels = rev(top_drugs$Drug))

# Handle edge case: all PI/CI are NA
if (all(is.na(top_drugs$PI_low))) {
  top_drugs$PI_low  <- top_drugs$CI_low
  top_drugs$PI_high <- top_drugs$CI_high
}

x_min <- min(top_drugs$PI_low, top_drugs$CI_low, na.rm = TRUE) - 0.05
x_max <- max(top_drugs$PI_high, top_drugs$CI_high, na.rm = TRUE) + 0.05
# Safety: if x_min == x_max (only 1 drug or all same), expand
if (x_max - x_min < 0.01) { x_min <- x_min - 0.1; x_max <- x_max + 0.1 }
text_zone <- (x_max - x_min) * 0.65
x_max_plot <- x_max + text_zone

col1_x <- x_max + text_zone * 0.08
col2_x <- x_max + text_zone * 0.38
col3_x <- x_max + text_zone * 0.72

p_forest <- ggplot(top_drugs, aes(x = Pooled_Diff, y = Display, color = Direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = PI_low, xmax = PI_high), height = 0.15, linewidth = 0.5, alpha = 0.35, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.35, linewidth = 1.0, na.rm = TRUE) +
  geom_point(size = 3.5, shape = 18) +
  geom_text(aes(label = sprintf("%+.3f", Pooled_Diff), x = col1_x),
            size = 2.6, family = "sans", color = "black", hjust = 0, fontface = "bold") +
  geom_text(aes(label = sprintf("[%.3f, %.3f]", CI_low, CI_high), x = col2_x),
            size = 2.4, family = "sans", color = "black", hjust = 0) +
  geom_text(aes(label = sprintf("%s  I2=%.0f%%", Sig, I2), x = col3_x),
            size = 2.4, family = "sans", color = "black", hjust = 0) +
  scale_color_manual(values = palette_dir, guide = "none") +
  scale_x_continuous(limits = c(x_min - 0.02, x_max_plot), expand = c(0, 0)) +
  labs(title = "Drug Sensitivity: BCR - Non-BCR IC50 Difference (Meta-Analysis, Top 30)",
       subtitle = paste0("GDSC2  |  ", length(unique(df_all$Cohort)), " cohorts pooled (REML)  |  ",
                         n_sig, " nominal sig  |  ", n_fdr, " FDR<0.05  |  Mean I2=",
                         sprintf("%.0f%%", mean_i2)),
       x = "Pooled BCR - Non-BCR IC50 Difference", y = NULL) +
  theme_arial +
  theme(panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3))

ggsave(file.path(fig_dir, "FigA_Drug_Meta_Forest.pdf"), p_forest, width = 15, height = 10)
ggsave(file.path(fig_dir, "FigA_Drug_Meta_Forest.png"), p_forest, width = 15, height = 10, dpi = 300)

# ---- Fig B: Volcano plot ----
cat("  [2/4] Volcano plot...\n")

drug_meta$logP <- -log10(drug_meta$pval)
drug_meta$label <- ifelse(drug_meta$pval < 0.01 & abs(drug_meta$Pooled_Diff) > 0.01,
                          drug_meta$Drug, "")

p_volcano <- ggplot(drug_meta, aes(x = Pooled_Diff, y = logP, color = Direction)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(aes(label = label), size = 2.2, family = "sans", color = "black",
                  max.overlaps = 20) +
  scale_color_manual(values = palette_dir, guide = guide_legend(title = "Sensitivity")) +
  labs(title = "Drug Sensitivity Meta-Analysis: Volcano Plot",
       subtitle = paste0(nrow(drug_meta), " drugs  |  ", n_sig, " nominal p<0.05  |  ",
                         n_fdr, " FDR<0.05"),
       x = "Pooled BCR - Non-BCR IC50 Difference", y = "-log10(p-value)") +
  theme_arial + theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "FigB_Drug_Volcano.pdf"), p_volcano, width = 8, height = 8)
ggsave(file.path(fig_dir, "FigB_Drug_Volcano.png"), p_volcano, width = 8, height = 8, dpi = 300)

# ---- Fig C: Heatmap of top drugs x cohorts ----
cat("  [3/4] Cohort concordance heatmap...\n")

# Only include cohorts that actually completed
completed_cohorts <- unique(drug_stats$Cohort)

topN <- min(20, nrow(drug_meta))
topN_drugs <- drug_meta$Drug[order(drug_meta$pval)][1:topN]
diff_mat <- tryCatch({
  mat <- dcast(subset(drug_stats, Drug %in% topN_drugs), Drug ~ Cohort, value.var = "Diff")
  rownames(mat) <- mat$Drug; mat$Drug <- NULL
  drug_order <- drug_meta$Drug[drug_meta$Drug %in% topN_drugs]
  mat <- mat[drug_order, completed_cohorts, drop = FALSE]
  as.matrix(mat)
}, error = function(e) {
  cat(sprintf("    [WARNING] dcast failed: %s. Generating simple heatmap.\n", e$message))
  NULL
})

if (!is.null(diff_mat) && nrow(diff_mat) > 0 && ncol(diff_mat) > 0) {
  sig_mat <- tryCatch({
    mat <- dcast(subset(drug_stats, Drug %in% topN_drugs), Drug ~ Cohort, value.var = "Sig")
    rownames(mat) <- mat$Drug; mat$Drug <- NULL
    mat <- mat[rownames(diff_mat), colnames(diff_mat), drop = FALSE]
    as.matrix(mat)
  }, error = function(e) matrix("", nrow = nrow(diff_mat), ncol = ncol(diff_mat),
                                 dimnames = dimnames(diff_mat)))
  
  rng <- max(abs(diff_mat), na.rm = TRUE)
  if (rng == 0) rng <- 1
  breaks <- seq(-rng, rng, length.out = 101)
  
  pdf(file.path(fig_dir, "FigC_Drug_Heatmap_Cohorts.pdf"), width = 9, height = 9)
  pheatmap(diff_mat,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = breaks,
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = sig_mat,
           fontsize = 11, fontsize_number = 13,
           fontfamily = "sans", number_color = "black",
           main = "", border_color = "white",
           cellwidth = 55, cellheight = 18, angle_col = 45,
           na_col = "grey90")
  grid.text(paste0("Drug Sensitivity: BCR - Non-BCR IC50 Difference Across Cohorts (Top ", topN, ")"),
            y = 0.87, x = 0.5, just = "center",
            gp = gpar(fontfamily = "sans", fontsize = 14, fontface = "bold", col = "black"))
  dev.off()
  
  png(file.path(fig_dir, "FigC_Drug_Heatmap_Cohorts.png"), width = 9, height = 9,
      units = "in", res = 300)
  pheatmap(diff_mat,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = breaks,
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = sig_mat,
           fontsize = 11, fontsize_number = 13,
           fontfamily = "sans", number_color = "black",
           main = "", border_color = "white",
           cellwidth = 55, cellheight = 18, angle_col = 45,
           na_col = "grey90")
  grid.text(paste0("Drug Sensitivity: BCR - Non-BCR IC50 Difference Across Cohorts (Top ", topN, ")"),
            y = 0.87, x = 0.5, just = "center",
            gp = gpar(fontfamily = "sans", fontsize = 14, fontface = "bold", col = "black"))
  dev.off()
} else {
  cat("  [SKIP] Heatmap: not enough data\n")
}

# ---- Fig D: Boxplots of top drugs ----
cat("  [4/4] Boxplots for top drugs...\n")

top6 <- drug_meta$Drug[order(drug_meta$pval)][1:min(6, nrow(drug_meta))]
available_top <- intersect(top6, colnames(df_all))
if (length(available_top) > 0) {
  df_top <- df_all[, c("Sample", "BCR", "Cohort", available_top), drop = FALSE]
  df_long <- pivot_longer(df_top, cols = all_of(available_top),
                          names_to = "Drug", values_to = "IC50")
  df_long$BCR <- factor(df_long$BCR, levels = c("Non-BCR", "BCR"))
  
  p_box <- ggplot(df_long, aes(x = BCR, y = IC50, fill = BCR)) +
    geom_boxplot(alpha = 0.8, width = 0.5, outlier.size = 0.3, outlier.alpha = 0.3) +
    facet_grid(Drug ~ Cohort, scales = "free_y") +
    scale_fill_manual(values = c("Non-BCR" = "#4DBBD5", "BCR" = "#E64B35"), guide = "none") +
    labs(title = paste0("Top ", length(available_top), " Differential Drugs: BCR vs Non-BCR Across Cohorts"),
         subtitle = "GDSC2 predicted IC50 | oncoPredict",
         x = NULL, y = "Predicted IC50") +
    theme_arial +
    theme(strip.text = element_text(size = 6.5), axis.text.x = element_text(size = 7))
  
  box_h <- min(12, 2.2 * length(available_top))
  ggsave(file.path(fig_dir, "FigD_TopDrugs_Boxplots.pdf"), p_box,
         width = 10, height = box_h, limitsize = FALSE)
  ggsave(file.path(fig_dir, "FigD_TopDrugs_Boxplots.png"), p_box,
         width = 10, height = box_h, dpi = 300, limitsize = FALSE)
} else {
  cat("  [SKIP] Boxplots: no top drugs available in data\n")
}

# ============================================================
# 6. Save + Report
# ============================================================
cat("\n=== Saving results ===\n\n")

write.csv(drug_stats, file.path(out_dir, "Drug_PerCohort_Stats.csv"), row.names = FALSE)
write.csv(drug_meta,  file.path(out_dir, "Drug_MetaAnalysis_Results.csv"), row.names = FALSE)

# Word report
doc <- read_docx()
doc <- body_add_par(doc, "oncoPredict Drug Sensitivity Prediction — Prostate Cancer Multi-Cohort Meta-Analysis (BCR vs Non-BCR)", "heading 1")
doc <- body_add_par(doc, paste("Date:", Sys.Date()), "Normal")
doc <- body_add_par(doc, "", "Normal")

doc <- body_add_par(doc, "1. Methods Overview", "heading 1")
doc <- body_add_par(doc, paste0(
  "oncoPredict (GDSC2 training set, ", nrow(traDataExp), " genes x ", ncol(traDataExp),
  " cell lines) was applied to ", length(expr_list), " independent prostate cancer cohorts (",
  paste(names(expr_list), collapse = ", "),
  ") to predict drug IC50. Within each cohort, BCR vs Non-BCR differences were compared ",
  "(Wilcoxon test), then integrated via Random-Effects Meta-Analysis (REML)."),
  "Normal")
doc <- body_add_par(doc, "", "Normal")

doc <- body_add_par(doc, "2. Overall Results", "heading 1")
doc <- body_add_par(doc, paste0("Total drugs tested: ", nrow(drug_meta)), "Normal")
doc <- body_add_par(doc, sprintf("Nominal p < 0.05: %d | FDR < 0.05: %d", n_sig, n_fdr), "Normal")
doc <- body_add_par(doc, sprintf("Mean I2 = %.1f%% (cross-cohort heterogeneity)", mean_i2), "Normal")
doc <- body_add_par(doc, "", "Normal")

doc <- body_add_par(doc, "3. Top Significant Drugs (FDR < 0.05)", "heading 1")
fdr_sig <- subset(drug_meta, pval_FDR < 0.05)
if (nrow(fdr_sig) > 0) {
  for (i in 1:nrow(fdr_sig)) {
    dir_label <- ifelse(fdr_sig$Direction[i] == "BCR_higher_IC50",
                        "BCR resistant (IC50 higher)", "BCR sensitive (IC50 lower)")
    doc <- body_add_par(doc,
      sprintf("  %s: %s  Diff=%+.4f [%.4f, %.4f] p=%.2e FDR=%.2e I2=%.0f%%",
              fdr_sig$Drug[i], dir_label, fdr_sig$Pooled_Diff[i],
              fdr_sig$CI_low[i], fdr_sig$CI_high[i],
              fdr_sig$pval[i], fdr_sig$pval_FDR[i], fdr_sig$I2[i]), "Normal")
  }
} else {
  doc <- body_add_par(doc, "  (No drugs passed FDR < 0.05 threshold)", "Normal")
}
doc <- body_add_par(doc, "", "Normal")

doc <- body_add_par(doc, "4. Figures", "heading 1")

if (file.exists(file.path(fig_dir, "FigA_Drug_Meta_Forest.png"))) {
  doc <- body_add_par(doc, "Fig A: Meta Forest Plot (Top 30)", "heading 2")
  doc <- body_add_img(doc, file.path(fig_dir, "FigA_Drug_Meta_Forest.png"), width = 5.5, height = 3.5)
}
if (file.exists(file.path(fig_dir, "FigB_Drug_Volcano.png"))) {
  doc <- body_add_par(doc, "Fig B: Volcano Plot", "heading 2")
  doc <- body_add_img(doc, file.path(fig_dir, "FigB_Drug_Volcano.png"), width = 4, height = 4)
}
if (file.exists(file.path(fig_dir, "FigC_Drug_Heatmap_Cohorts.png"))) {
  doc <- body_add_par(doc, "Fig C: Cohort Concordance Heatmap", "heading 2")
  doc <- body_add_img(doc, file.path(fig_dir, "FigC_Drug_Heatmap_Cohorts.png"), width = 4.5, height = 4.5)
}
if (file.exists(file.path(fig_dir, "FigD_TopDrugs_Boxplots.png"))) {
  doc <- body_add_par(doc, "Fig D: Top Drug Boxplots", "heading 2")
  doc <- body_add_img(doc, file.path(fig_dir, "FigD_TopDrugs_Boxplots.png"), width = 5, height = 3)
}

doc <- body_add_par(doc, "", "Normal")
doc <- body_add_par(doc, paste("=== Report End |", Sys.Date(), "==="), "Normal")
report_path <- file.path(out_dir, "oncoPredict_BCR_DrugReport.docx")
print(doc, target = report_path)
cat(sprintf("Report: %s\n", report_path))

# ---- Close log ----
cat("\n========================================\n")
cat("COMPLETE\n")
cat("========================================\n")
cat(sprintf("Drugs: %d | Nominal sig: %d | FDR sig: %d | Mean I2: %.1f%%\n",
            nrow(drug_meta), n_sig, n_fdr, mean_i2))
cat(sprintf("Output: %s\n", out_dir))
cat(sprintf("Report: %s\n", report_path))
cat(sprintf("Log:   %s\n", log_file))

sink(type = "message")
sink()
close(log_con)

