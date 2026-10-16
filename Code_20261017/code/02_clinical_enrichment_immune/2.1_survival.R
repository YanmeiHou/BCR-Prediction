# ==== 2.1 生存分析 ====
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
# Survival Analysis — 4 Genes × 4 Cohorts
# BCR-free survival (all 4 cohorts) + OS (TCGA only)
# Output: BCR_ML_Publication/Survival/
################################################################################

library(Biobase)
library(survival)
library(survminer)

setwd(legacy_output())

core_ensembl <- c(
  "ENSG00000124721" = "DNAH8",
  "ENSG00000131747" = "TOP2A",
  "ENSG00000101938" = "CHRDL1",
  "ENSG00000154330" = "PGM5"
)

datasets <- list(
  "DKFZ"      = readRDS(project_path("data", "external", "modelling", "DKFZ_eSet.RDS")),
  "TCGA"      = readRDS(project_path("data", "external", "modelling", "TCGA-PRAD_eSet.RDS")),
  "Taylor"    = readRDS(project_path("data", "external", "modelling", "Taylor_eSet.RDS")),
  "CancerMap" = readRDS(project_path("data", "external", "modelling", "CancerMap_eSet.RDS"))
)

pal_gene <- c("DNAH8"="#E41A1C", "TOP2A"="#377EB8", "CHRDL1"="#4DAF4A", "PGM5"="#984EA3")

# ===========================================================================
# Function: KM plot for one gene × one cohort
# ===========================================================================
# ---- 2.1.1 km_gene_cohort ----
km_gene_cohort <- function(eset, cohort_name, gene_ensembl, gene_symbol,
                           surv_type = "BCR") {

  expr <- exprs(eset)
  pd   <- pData(eset)

  # Check gene exists
  if (!gene_ensembl %in% rownames(expr)) return(NULL)

  gene_expr <- expr[gene_ensembl]

  # Build survival data
  if (surv_type == "BCR") {
    time_col <- "time_to_bcr"
    event_col <- "bcr_status"
  } else {
    time_col <- "time_to_death"
    event_col <- "os_status"
  }

  # Check columns exist and have data
  if (!all(c(time_col, event_col) %in% colnames(pd))) return(NULL)

  time_vals <- pd[[time_col]]
  event_vals <- pd[[event_col]]

  # Remove NAs
  valid <- which(!is.na(time_vals) & !is.na(event_vals))
  if (length(valid) < 10) return(NULL)

  time_vals <- as.numeric(time_vals[valid])
  event_vals <- as.numeric(event_vals[valid])
  gene_expr  <- gene_expr[valid]

  # Median split
  cutoff <- median(gene_expr, na.rm=TRUE)
  group <- ifelse(gene_expr > cutoff, "High", "Low")
  group <- factor(group, levels=c("Low", "High"))

  # Create surv object
  surv_obj <- Surv(time_vals, event_vals)

  # Fit KM
  fit <- survfit(surv_obj ~ group)

  # P-value
  diff_test <- survdiff(surv_obj ~ group)
  pval <- 1 - pchisq(diff_test$chisq, df=1)

  # Plot
  g <- ggsurvplot(fit, data=data.frame(time=time_vals, event=event_vals, group=group),
                  pval = TRUE, pval.method = TRUE,
                  conf.int = TRUE,
                  palette = c("#2E86AB", "#D64045"),
                  legend.labs = c("Low Expression", "High Expression"),
                  legend.title = gene_symbol,
                  title = paste0(gene_symbol, " — ", cohort_name, " (", surv_type, ")"),
                  xlab = ifelse(surv_type=="BCR", "Time to BCR (months)", "Overall Survival (months)"),
                  ylab = paste0(surv_type, "-free Survival Probability"),
                  risk.table = TRUE,
                  risk.table.height = 0.25,
                  ggtheme = theme_classic(base_family="ArialMT"),
                  font.title = c(14, "bold"),
                  font.x = c(12), font.y = c(12),
                  font.legend = c(11),
                  font.tickslab = c(10))

  return(list(plot=g, fit=fit, pval=pval, n=nrow(fit)))
}

# ===========================================================================
# Generate all KM curves
# ===========================================================================
cat("Generating survival curves...\n")
dir.create("Survival", showWarnings=FALSE)

for (ds_name in names(datasets)) {
  eset <- datasets[[ds_name]]

  for (gid in names(core_ensembl)) {
    gsym <- core_ensembl[gid]

    # BCR-free survival
    res_bcr <- km_gene_cohort(eset, ds_name, gid, gsym, "BCR")
    if (!is.null(res_bcr)) {
      pdf(paste0("Survival/KM_BCR_", gsym, "_", ds_name, ".pdf"),
          width=8, height=7, family="ArialMT")
      print(res_bcr$plot)
      dev.off()
      png(paste0("Survival/KM_BCR_", gsym, "_", ds_name, ".png"),
          width=8*150, height=7*150, res=150)
      print(res_bcr$plot)
      dev.off()
      cat(sprintf("  BCR %s %s: p=%.4f\n", gsym, ds_name, res_bcr$pval))
    }

    # OS (TCGA only)
    if (ds_name == "TCGA") {
      res_os <- km_gene_cohort(eset, ds_name, gid, gsym, "OS")
      if (!is.null(res_os)) {
        pdf(paste0("Survival/KM_OS_", gsym, "_TCGA.pdf"),
            width=8, height=7, family="ArialMT")
        print(res_os$plot)
        dev.off()
        png(paste0("Survival/KM_OS_", gsym, "_TCGA.png"),
            width=8*150, height=7*150, res=150)
        print(res_os$plot)
        dev.off()
        cat(sprintf("  OS  %s TCGA: p=%.4f\n", gsym, res_os$pval))
      }
    }
  }
}

# ===========================================================================
# Combined figure: 4 genes × 4 cohorts BCR KM
# ===========================================================================
cat("\nGenerating combined figure...\n")

pdf("Survival/Fig8_KM_BCR_AllGenes_AllCohorts.pdf", width=20, height=18, family="ArialMT")
par(mfrow=c(4,4), mar=c(4,4,2.5,1), oma=c(2,2,3,1))

for (gid in names(core_ensembl)) {
  gsym <- core_ensembl[gid]
  for (ds_name in names(datasets)) {
    eset <- datasets[[ds_name]]
    expr <- exprs(eset)
    pd <- pData(eset)

    if (!gid %in% rownames(expr)) { plot.new(); next }

    gene_expr <- expr[gid]
    time_vals <- as.numeric(pd$time_to_bcr)
    event_vals <- as.numeric(pd$bcr_status)

    valid <- which(!is.na(time_vals) & !is.na(event_vals))
    if (length(valid) < 10) { plot.new(); next }

    time_vals <- time_vals[valid]
    event_vals <- event_vals[valid]
    gene_expr  <- gene_expr[valid]

    cutoff <- median(gene_expr)
    group <- factor(ifelse(gene_expr > cutoff, "High", "Low"), levels=c("Low", "High"))

    fit <- survfit(Surv(time_vals, event_vals) ~ group)
    diff_test <- survdiff(Surv(time_vals, event_vals) ~ group)
    pval <- 1 - pchisq(diff_test$chisq, df=1)

    plot(fit, col=c("#2E86AB", "#D64045"), lwd=2.5,
         xlab="Time (months)", ylab="BCR-free Survival",
         main=paste(gsym, "-", ds_name),
         cex.main=1.1, cex.lab=1, cex.axis=0.9)
    legend("bottomleft", c("Low", "High"), col=c("#2E86AB", "#D64045"),
           lwd=2.5, bty="n", cex=0.7, title=gsym)
    mtext(paste0("P = ", format(pval, digits=3, scientific=TRUE)),
          side=3, line=-1.2, cex=0.65, col="grey30")
  }
}
mtext("BCR-Free Survival — 4 Genes × 4 Cohorts", outer=TRUE, side=3,
      cex=1.5, font=2, line=0.5)
dev.off()

png("Survival/Fig8_KM_BCR_AllGenes_AllCohorts.png", width=20*150, height=18*150, res=150)
par(mfrow=c(4,4), mar=c(4,4,2.5,1), oma=c(2,2,3,1))
for (gid in names(core_ensembl)) {
  gsym <- core_ensembl[gid]
  for (ds_name in names(datasets)) {
    eset <- datasets[[ds_name]]
    expr <- exprs(eset)
    pd <- pData(eset)

    if (!gid %in% rownames(expr)) { plot.new(); next }
    gene_expr <- expr[gid]
    time_vals <- as.numeric(pd$time_to_bcr)
    event_vals <- as.numeric(pd$bcr_status)
    valid <- which(!is.na(time_vals) & !is.na(event_vals))
    if (length(valid) < 10) { plot.new(); next }
    time_vals <- time_vals[valid]; event_vals <- event_vals[valid]
    gene_expr <- gene_expr[valid]
    cutoff <- median(gene_expr)
    group <- factor(ifelse(gene_expr > cutoff, "High", "Low"), levels=c("Low", "High"))
    fit <- survfit(Surv(time_vals, event_vals) ~ group)
    diff_test <- survdiff(Surv(time_vals, event_vals) ~ group)
    pval <- 1 - pchisq(diff_test$chisq, df=1)
    plot(fit, col=c("#2E86AB", "#D64045"), lwd=2.5,
         xlab="Time (months)", ylab="BCR-free Survival",
         main=paste(gsym, "-", ds_name),
         cex.main=1.1, cex.lab=1, cex.axis=0.9)
    legend("bottomleft", c("Low", "High"), col=c("#2E86AB", "#D64045"),
           lwd=2.5, bty="n", cex=0.7, title=gsym)
    mtext(paste0("P = ", format(pval, digits=3, scientific=TRUE)),
          side=3, line=-1.2, cex=0.65, col="grey30")
  }
}
mtext("BCR-Free Survival — 4 Genes × 4 Cohorts", outer=TRUE, side=3,
      cex=1.5, font=2, line=0.5)
dev.off()

# ===========================================================================
# Summary table
# ===========================================================================
cat("\n=== BCR-Free Survival P-values (Median Split) ===\n")
cat(sprintf("%-10s %-12s %-12s %-12s %-12s\n", "Gene", "DKFZ", "TCGA", "Taylor", "CancerMap"))
cat(strrep("-", 58), "\n")

for (gid in names(core_ensembl)) {
  gsym <- core_ensembl[gid]
  pvals <- c()
  for (ds_name in names(datasets)) {
    eset <- datasets[[ds_name]]
    expr <- exprs(eset)
    pd <- pData(eset)
    if (!gid %in% rownames(expr)) { pvals <- c(pvals, NA); next }
    gene_expr <- expr[gid]
    time_vals <- as.numeric(pd$time_to_bcr)
    event_vals <- as.numeric(pd$bcr_status)
    valid <- which(!is.na(time_vals) & !is.na(event_vals))
    if (length(valid) < 10) { pvals <- c(pvals, NA); next }
    group <- factor(ifelse(gene_expr[valid] > median(gene_expr[valid]), "High", "Low"))
    diff_test <- survdiff(Surv(time_vals[valid], event_vals[valid]) ~ group)
    pvals <- c(pvals, 1 - pchisq(diff_test$chisq, df=1))
  }
  cat(sprintf("%-10s %-12s %-12s %-12s %-12s\n", gsym,
              format(pvals[1], digits=3, scientific=TRUE),
              format(pvals[2], digits=3, scientific=TRUE),
              format(pvals[3], digits=3, scientific=TRUE),
              format(pvals[4], digits=3, scientific=TRUE)))
}

cat("\n=== OS P-values (TCGA only, 10 events) ===\n")
for (gid in names(core_ensembl)) {
  gsym <- core_ensembl[gid]
  eset <- datasets[["TCGA"]]
  expr <- exprs(eset)
  pd <- pData(eset)
  gene_expr <- expr[gid]
  time_vals <- as.numeric(pd$time_to_death)
  event_vals <- as.numeric(pd$os_status)
  valid <- which(!is.na(time_vals) & !is.na(event_vals))
  group <- factor(ifelse(gene_expr[valid] > median(gene_expr[valid]), "High", "Low"))
  diff_test <- survdiff(Surv(time_vals[valid], event_vals[valid]) ~ group)
  pval <- 1 - pchisq(diff_test$chisq, df=1)
  cat(sprintf("  %s: p = %s\n", gsym, format(pval, digits=3, scientific=TRUE)))
}

cat("\nAll survival curves saved to: Survival/\n")
cat("Done!\n")
