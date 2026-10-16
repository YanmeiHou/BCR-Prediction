# ==== 2.6 基因与免疫细胞相关性 ====
# Execution scope: 历史脚本；需对应输入
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
# 3 Figures: ESTIMATE 4-metrics + CIBERSORT Lymphoid + CIBERSORT Myeloid
# Arial, pure black, black border, 2x2 layout
################################################################################
library(Biobase)
library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)

out_dir <- file.path(legacy_output(), "Immune_Results")
FONT <- "sans"

# ---- 2.6.1 Load ----
cibersort_res <- readRDS(file.path(out_dir, "cibersort_res.RDS"))
estimate_scores <- readRDS(file.path(out_dir, "estimate_scores.RDS"))

# Compute TumorPurity
# ---- 2.6.2 compute_purity ----
compute_purity <- function(es) { p <- cos(0.6049872018+0.0001467884*es); p[p<0]<-NA; p }
for(ds in names(estimate_scores)) {
  estimate_scores[[ds]]$TumorPurity <- compute_purity(estimate_scores[[ds]]$ESTIMATEScore)
}

datasets <- list(
  DKFZ=readRDS(project_path("data", "external", "modelling", "DKFZ_eSet.RDS")),
  TCGA=readRDS(project_path("data", "external", "modelling", "TCGA-PRAD_eSet.RDS")),
  Taylor=readRDS(project_path("data", "external", "modelling", "Taylor_eSet.RDS")),
  CancerMap=readRDS(project_path("data", "external", "modelling", "CancerMap_eSet.RDS"))
)

core_ensembl <- c("ENSG00000124721"="DNAH8","ENSG00000131747"="TOP2A",
                  "ENSG00000101938"="CHRDL1","ENSG00000154330"="PGM5")
genes <- unname(core_ensembl)

# ---- 2.6.3 Helpers ----
# ---- 2.6.4 cor_with_pval ----
cor_with_pval <- function(x, y) {
  r_mat <- matrix(NA, ncol(x), ncol(y))
  p_mat <- matrix(NA, ncol(x), ncol(y))
  rownames(r_mat) <- colnames(x); colnames(r_mat) <- colnames(y)
  rownames(p_mat) <- colnames(x); colnames(p_mat) <- colnames(y)
  for(i in seq_len(ncol(x))) for(j in seq_len(ncol(y))) {
    ok <- complete.cases(x[,i], y[,j])
    if(sum(ok) >= 5) {
      ct <- cor.test(x[ok,i], y[ok,j], method="spearman")
      r_mat[i,j] <- ct$estimate
      p_mat[i,j] <- ct$p.value
    }
  }
  list(r=r_mat, p=p_mat)
}

# ---- 2.6.5 match_samples ----
match_samples <- function(s1, s2) {
  cm <- intersect(s1, s2)
  if(length(cm) > 0) return(cm)
  s1n <- gsub("\\.","-",s1); s2n <- gsub("\\.","-",s2)
  mp <- intersect(s1n, s2n)
  if(length(mp) > 0) return(list(idx1=match(mp,s1n), idx2=match(mp,s2n), names=mp))
  return(character(0))
}

# ---- 2.6.6 Compute correlations ----
all_ciber <- list()
all_est   <- list()

for(ds in names(datasets)) {
  eset <- datasets[[ds]]
  expr_full <- exprs(eset)
  pd <- pData(eset)
  bcr <- pd$bcr_status
  keep <- which(!is.na(bcr) & bcr %in% c(0,1))
  common_genes <- intersect(names(core_ensembl), rownames(expr_full))
  expr_4g <- t(expr_full[common_genes, keep, drop=FALSE])
  colnames(expr_4g) <- core_ensembl[common_genes]

  # CIBERSORT
  ciber <- cibersort_res[[ds]]
  res <- match_samples(rownames(expr_4g), rownames(ciber))
  if(is.list(res)) {
    cr <- cor_with_pval(expr_4g[res$idx1,,drop=FALSE], as.matrix(ciber)[res$idx2,,drop=FALSE])
  } else {
    cr <- cor_with_pval(expr_4g[res,,drop=FALSE], as.matrix(ciber)[res,,drop=FALSE])
  }
  all_ciber[[ds]] <- cr

  # ESTIMATE
  est <- estimate_scores[[ds]]
  res2 <- match_samples(rownames(expr_4g), rownames(est))
  if(is.list(res2)) {
    er <- cor_with_pval(expr_4g[res2$idx1,,drop=FALSE], as.matrix(est)[res2$idx2,,drop=FALSE])
  } else {
    er <- cor_with_pval(expr_4g[res2,,drop=FALSE], as.matrix(est)[res2,,drop=FALSE])
  }
  all_est[[ds]] <- er
}

# ---- 2.6.7 Build long-format data ----
# ---- 2.6.8 build_df ----
build_df <- function(cor_list) {
  df_list <- list()
  for(ds in names(cor_list)) {
    r_mat <- cor_list[[ds]]$r
    p_mat <- cor_list[[ds]]$p
    r_long <- melt(r_mat, varnames=c("Gene","Variable"), value.name="Correlation")
    p_long <- melt(p_mat, varnames=c("Gene","Variable"), value.name="p_value")
    df <- merge(r_long, p_long, by=c("Gene","Variable"))
    df$sig <- ifelse(is.na(df$p_value), "",
              ifelse(df$p_value < 0.001, "***",
              ifelse(df$p_value < 0.01, "**",
              ifelse(df$p_value < 0.05, "*", ""))))
    df$Dataset <- ds
    df$label <- sprintf("%.2f%s", df$Correlation, df$sig)
    df_list[[ds]] <- df
  }
  do.call(rbind, df_list)
}

df_ciber <- build_df(all_ciber)
df_est   <- build_df(all_est)

# Clean CIBERSORT cell names
# ---- 2.6.9 clean_cell ----
clean_cell <- function(x) {
  x <- gsub("T cells follicular helper", "Tfh", x)
  x <- gsub("T cells regulatory \\(Tregs\\)", "Treg", x)
  x <- gsub("T cells CD4 memory resting", "CD4(R)", x)
  x <- gsub("T cells CD4 memory activated", "CD4(A)", x)
  x <- gsub("T cells CD4 naive", "CD4(n)", x)
  x <- gsub("T cells CD8", "CD8", x)
  x <- gsub("T cells gamma delta", "T(gd)", x)
  x <- gsub("B cells naive", "B naive", x)
  x <- gsub("B cells memory", "B memory", x)
  x <- gsub("NK cells resting", "NK(R)", x)
  x <- gsub("NK cells activated", "NK(A)", x)
  x <- gsub("Plasma cells", "Plasma", x)
  x <- gsub("Macrophages M0", "M0", x)
  x <- gsub("Macrophages M1", "M1", x)
  x <- gsub("Macrophages M2", "M2", x)
  x <- gsub("Monocytes", "Mono", x)
  x <- gsub("Mast cells resting", "Mast(R)", x)
  x <- gsub("Mast cells activated", "Mast(A)", x)
  x <- gsub("Dendritic cells resting", "DC(R)", x)
  x <- gsub("Dendritic cells activated", "DC(A)", x)
  x <- gsub("Eosinophils", "Eosino", x)
  x <- gsub("Neutrophils", "Neutro", x)
  x
}

df_ciber$Variable <- clean_cell(as.character(df_ciber$Variable))

# Order
dataset_order <- c("TCGA","DKFZ","Taylor","CancerMap")
df_ciber$Gene <- factor(df_ciber$Gene, levels=genes)
df_ciber$Dataset <- factor(df_ciber$Dataset, levels=dataset_order)
df_est$Gene <- factor(df_est$Gene, levels=genes)
df_est$Dataset <- factor(df_est$Dataset, levels=dataset_order)
df_est$Variable <- factor(df_est$Variable, levels=c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"))

# ---- 2.6.10 Theme ----
base_theme <- theme_minimal(base_size=10) +
  theme(
    text = element_text(color="black", family=FONT),
    axis.text = element_text(color="black", family=FONT),
    axis.title = element_text(color="black", family=FONT),
    plot.title = element_text(color="black", family=FONT, face="bold"),
    strip.text = element_text(color="black", family=FONT, face="bold"),
    legend.text = element_text(color="black", family=FONT),
    legend.title = element_text(color="black", family=FONT),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# ---- 2.6.11 make_heatmap ----
make_heatmap <- function(df, title, tsz=2.8, xa=45, xsz=7.5) {
  ggplot(df, aes(x=Variable, y=Gene, fill=Correlation)) +
    geom_tile(color="black", linewidth=0.6) +
    geom_text(aes(label=label), size=tsz, fontface="bold", family=FONT, color="black") +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                         midpoint=0, limits=c(-1,1), name="Spearman R") +
    facet_wrap(~Dataset, nrow=2, ncol=2) +
    labs(title=title, x="", y="") +
    base_theme +
    theme(
      axis.text.x = element_text(angle=xa, hjust=if(xa==45) 1 else 0.5, size=xsz, face="bold"),
      axis.text.y = element_text(size=9, face="bold"),
      plot.title = element_text(size=13, hjust=0.5),
      strip.text = element_text(size=10),
      legend.key.height = unit(0.15, "npc")
    ) +
    coord_fixed(ratio=1)
}

# ---- 2.6.12 Split CIBERSORT cells ----
all_cells <- sort(unique(as.character(df_ciber$Variable)))
cat("All CIBERSORT cells (", length(all_cells), "):\n", paste(all_cells, collapse=", "), "\n\n")

# Lymphoid group (B, T, NK, Plasma)
lymphoid <- c("B naive","B memory","Plasma","CD8","CD4(R)","CD4(A)","Tfh","Treg","NK(R)")
lymphoid <- intersect(lymphoid, all_cells)

# Myeloid + others (Mono, M0-2, DC, Mast, Eosino, Neutro, NK(A), T(gd))
myeloid <- setdiff(all_cells, lymphoid)

cat("Lymphoid:", paste(lymphoid, collapse=", "), "\n")
cat("Myeloid:", paste(myeloid, collapse=", "), "\n")

# ---- 2.6.13 Generate figures ----
# Figure 1: ESTIMATE
p_est <- make_heatmap(df_est, "4-Gene Signature vs TME Scores", 3.5, 45, 9)
pdf(file.path(out_dir, "Fig_ESTIMATE_4metrics.pdf"), width=12, height=8)
print(p_est)
dev.off()
png(file.path(out_dir, "Fig_ESTIMATE_4metrics.png"), width=12*150, height=8*150, res=150)
print(p_est)
dev.off()
cat("ESTIMATE saved.\n")

# Figure 2: CIBERSORT Lymphoid
df_ly <- subset(df_ciber, Variable %in% lymphoid)
df_ly$Variable <- factor(df_ly$Variable, levels=lymphoid)
cat("Lymphoid df rows:", nrow(df_ly), "\n")

p_ly <- make_heatmap(df_ly, "4-Gene vs Lymphoid Cells (CIBERSORT)", 2.6, 45, 7)
pdf(file.path(out_dir, "Fig_CIBERSORT_Lymphoid.pdf"), width=12, height=8)
print(p_ly)
dev.off()
png(file.path(out_dir, "Fig_CIBERSORT_Lymphoid.png"), width=12*150, height=8*150, res=150)
print(p_ly)
dev.off()
cat("CIBERSORT Lymphoid saved.\n")

# Figure 3: CIBERSORT Myeloid
df_my <- subset(df_ciber, Variable %in% myeloid)
df_my$Variable <- factor(df_my$Variable, levels=myeloid)
cat("Myeloid df rows:", nrow(df_my), "\n")

p_my <- make_heatmap(df_my, "4-Gene vs Myeloid & Other Cells (CIBERSORT)", 2.8, 45, 7)
pdf(file.path(out_dir, "Fig_CIBERSORT_Myeloid.pdf"), width=12, height=8)
print(p_my)
dev.off()
png(file.path(out_dir, "Fig_CIBERSORT_Myeloid.png"), width=12*150, height=8*150, res=150)
print(p_my)
dev.off()
cat("CIBERSORT Myeloid saved.\n")

cat("\nAll done!\n")
