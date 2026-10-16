# ==== 2.2 Cox回归与森林图 ====
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
# Univariate Cox + Forest Plot (forestploter package)
# Output: Cox/Fig12_ForestPlot.pdf, Cox/univariate_cox_results.csv
################################################################################

library(Biobase)
library(survival)
library(forestploter)

setwd(legacy_output())
dir.create("Cox", showWarnings=FALSE)

# ===========================================================================
# ---- 2.2.1 Cox regression: 4 genes x 4 cohorts ----
# ===========================================================================
core_ensembl <- c("ENSG00000124721"="DNAH8", "ENSG00000131747"="TOP2A",
                  "ENSG00000101938"="CHRDL1", "ENSG00000154330"="PGM5")

datasets <- list("DKFZ"      = readRDS(project_path("data", "external", "modelling", "DKFZ_eSet.RDS")),
                 "TCGA"      = readRDS(project_path("data", "external", "modelling", "TCGA-PRAD_eSet.RDS")),
                 "Taylor"    = readRDS(project_path("data", "external", "modelling", "Taylor_eSet.RDS")),
                 "CancerMap" = readRDS(project_path("data", "external", "modelling", "CancerMap_eSet.RDS")))

genes   <- c("DNAH8","TOP2A","CHRDL1","PGM5")
cohorts <- c("DKFZ","TCGA","Taylor","CancerMap")

HR <- LO <- HI <- PV <- Nmat <- Emat <- matrix(NA,4,4,dimnames=list(genes,cohorts))

for(gid in names(core_ensembl)) {
  gsym <- core_ensembl[gid]
  for(ds in names(datasets)) {
    eset <- datasets[[ds]]; pd <- pData(eset)
    tb <- as.numeric(pd$time_to_bcr); eb <- as.numeric(pd$bcr_status)
    valid <- which(!is.na(tb) & !is.na(eb))
    if(!gid %in% rownames(exprs(eset))) next
    ge <- exprs(eset)[gid, valid]; tb <- tb[valid]; eb <- eb[valid]

    cox <- coxph(Surv(tb, eb) ~ ge)
    s <- summary(cox); ci <- suppressMessages(exp(confint(cox)))
    HR[gsym,ds]  <- exp(coef(cox))
    LO[gsym,ds]  <- ci[1]; HI[gsym,ds] <- ci[2]
    PV[gsym,ds]  <- s$coefficients[1,"Pr(>|z|)"]
    Nmat[gsym,ds] <- length(valid); Emat[gsym,ds] <- sum(eb)
  }
}

# Save CSV
res <- data.frame(stringsAsFactors=FALSE)
for(g in genes) for(cn in cohorts) if(!is.na(HR[g,cn]))
  res <- rbind(res, data.frame(Gene=g, Cohort=cn, N=Nmat[g,cn], Event=Emat[g,cn],
            HR=HR[g,cn], HR_lower=LO[g,cn], HR_upper=HI[g,cn], P=PV[g,cn], row.names=NULL))
write.csv(res, "Cox/univariate_cox_results.csv", row.names=FALSE)

# Print table
cat("=== Univariate Cox Regression ===\n")
cat(sprintf("%-8s %-12s %5s %5s  %-18s %s\n","Gene","Cohort","N","Event","HR (95% CI)","P"))
cat(strrep("-",70),"\n")
for(i in 1:nrow(res)) {
  stars <- if(res$P[i]<0.001)"***"else if(res$P[i]<0.01)"**"else if(res$P[i]<0.05)"*"else""
  cat(sprintf("%-8s %-12s %5d %5d  %.2f (%.2f-%.2f) %s   %s\n",
      res$Gene[i],res$Cohort[i],res$N[i],res$Event[i],
      res$HR[i],res$HR_lower[i],res$HR_upper[i],
      if(res$P[i]<0.001)"<0.001"else sprintf("%.3f",res$P[i]),stars))
}

# ===========================================================================
# ---- 2.2.2 Build forestploter data frame ----
# ===========================================================================
df <- data.frame(stringsAsFactors=FALSE)

for(g in genes) {
  for(cn in cohorts) {
    hr <- HR[g,cn]; lo <- LO[g,cn]; hi <- HI[g,cn]; pv <- PV[g,cn]
    if(is.na(hr)) next

    stars <- if(pv<0.001)"***"else if(pv<0.01)"**"else if(pv<0.05)"*"else""
    p_str <- if(pv<0.001)"<0.001"else sprintf("%.3f",pv)

    df <- rbind(df, data.frame(
      Gene    = if(cn==cohorts[1]) g else "",
      Cohort  = cn,
      N       = as.character(Nmat[g,cn]),
      Event   = as.character(Emat[g,cn]),
      HR_text = sprintf("%.2f (%.2f-%.2f) %s", hr, lo, hi, stars),
      Pval    = p_str,
      ` `     = paste(rep(" ", 20), collapse=" "),
      HR      = hr,
      Lower   = lo,
      Upper   = hi,
      stringsAsFactors = FALSE, check.names = FALSE))
  }
}

# ===========================================================================
# ---- 2.2.3 Forest plot ----
# ===========================================================================
p <- forest(
  df[, c("Gene","Cohort","N","Event","HR_text","Pval"," ")],
  est       = df$HR,
  lower     = df$Lower,
  upper     = df$Upper,
  ci_column = 7,
  ref_line  = 1,
  xlim      = c(0.08, 5.5),
  xlab      = "Hazard Ratio (95% CI)",
  title     = "Univariate Cox Regression \u2014 BCR-Free Survival",
  ticks_at  = c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5)
)

pdf("Cox/Fig12_ForestPlot.pdf", width=14, height=8.5)
print(p)
dev.off()

png("Cox/Fig12_ForestPlot.png", width=14*150, height=8.5*150, res=150)
print(p)
dev.off()

cat("\nDone: Cox/Fig12_ForestPlot.pdf (forestploter)\n")
