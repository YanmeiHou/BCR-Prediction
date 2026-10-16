# ==== 2.4 单基因GSEA ====
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
# Single-Gene GSEA Analysis — Hallmark Pathways (Optimized)
#
# For each of the 4 core genes (DNAH8, TOP2A, CHRDL1, PGM5):
#   - Rank all genes by Spearman correlation with the gene of interest
#   - Run GSEA against HALLMARK gene sets
#   - Generate publication-quality plots
#
# Key optimizations:
#   - Ensembl → Symbol mapping via org.Hs.eg.db
#   - cor() for vectorized Spearman correlation (much faster than apply)
#
# Output: SingleGene_GSEA/
################################################################################

library(Biobase)
library(fgsea)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd(legacy_output())
dir.create("SingleGene_GSEA", showWarnings = FALSE)

# ===========================================================================
# ---- 2.4.1 Parameters ----
# ===========================================================================
core_genes <- c("DNAH8", "TOP2A", "CHRDL1", "PGM5")
core_ensembl <- c(
  "ENSG00000124721" = "DNAH8",
  "ENSG00000131747" = "TOP2A",
  "ENSG00000101938" = "CHRDL1",
  "ENSG00000154330" = "PGM5"
)

pal_gene <- c("DNAH8" = "#E41A1C", "TOP2A" = "#377EB8",
              "CHRDL1" = "#4DAF4A", "PGM5" = "#984EA3")

# ===========================================================================
# ---- 2.4.2 Load & Map Expression Data ----
# ===========================================================================
cat("═══ Loading TCGA-PRAD ═══\n")
eset_tcga <- readRDS(project_path("data", "external", "modelling", "TCGA-PRAD_eSet.RDS"))
expr_tcga <- exprs(eset_tcga)
cat(sprintf("%d genes x %d samples\n", nrow(expr_tcga), ncol(expr_tcga)))

cat("\n═══ Mapping Ensembl → Symbol ═══\n")
ensembl_ids <- rownames(expr_tcga)
symbol_map <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                     column = "SYMBOL", keytype = "ENSEMBL",
                     multiVals = "first")

keep <- which(!is.na(symbol_map))
nc <- ncol(expr_tcga)
expr_mapped <- expr_tcga[keep, TRUE]
symbols_mapped <- symbol_map[keep]
cat(sprintf("Mapped: %d / %d genes\n", length(symbols_mapped), length(ensembl_ids)))

# Deduplicate symbols
dup_symbols <- names(which(table(symbols_mapped) > 1))
if (length(dup_symbols) > 0) {
  row_vars <- apply(expr_mapped, 1, var)
  rows_to_drop <- integer(0)
  for (sym in dup_symbols) {
    idx <- which(symbols_mapped == sym)
    best_idx <- idx[which.max(row_vars[idx])]
    rows_to_drop <- c(rows_to_drop, setdiff(idx, best_idx))
  }
  if (length(rows_to_drop) > 0) {
    keep_rows <- setdiff(seq_len(nrow(expr_mapped)), rows_to_drop)
    expr_mapped <- expr_mapped[keep_rows, TRUE]
    symbols_mapped <- symbols_mapped[keep_rows]
  }
}
rownames(expr_mapped) <- symbols_mapped
cat(sprintf("After dedup: %d unique symbols\n", nrow(expr_mapped)))

# ===========================================================================
# ---- 2.4.3 Load Hallmark Gene Sets ----
# ===========================================================================
cat("\n═══ Loading Hallmark GMT ═══\n")
gmt_path <- project_path("data", "external", "h.all.v2026.1.Hs.symbols.gmt")
hallmark_pathways <- gmtPathways(gmt_path)
cat(sprintf("%d pathways loaded\n", length(hallmark_pathways)))

hallmark_genes <- unique(unlist(hallmark_pathways))
overlap <- intersect(hallmark_genes, rownames(expr_mapped))
cat(sprintf("Overlap: %d / %d Hallmark genes in expression\n",
            length(overlap), length(hallmark_genes)))

# ===========================================================================
# ---- 2.4.4 Compute Spearman Correlations (vectorized) ----
# ===========================================================================
cat("\n═══ Computing Spearman Correlations ═══\n")

# Transpose to samples × genes for cor()
expr_t <- t(expr_mapped)
gene_expr_subset <- expr_t[, core_genes, drop = FALSE]

system.time({
  all_cors <- cor(expr_t, gene_expr_subset, method = "spearman")
})
cat(sprintf("Correlation matrix: %d genes x %d query genes\n",
            nrow(all_cors), ncol(all_cors)))

# ===========================================================================
# ---- 2.4.5 Run fgsea for Each Gene ----
# ===========================================================================
cat("\n═══ Running GSEA ═══\n")

gsea_results <- list()
gsea_tables  <- list()
ranked_lists <- list()

for (g in core_genes) {
  cat(sprintf("\n─── %s ───\n", g))

  # Extract ranked list
  cors <- all_cors[, g, drop = TRUE]
  names(cors) <- rownames(all_cors)
  cors <- sort(cors[!is.na(cors)], decreasing = TRUE)
  cat(sprintf("  Ranked genes: %d (corr range: %.3f ~ %.3f)\n",
              length(cors), min(cors), max(cors)))

  ranked_lists[[g]] <- cors

  # Run fgsea
  set.seed(777)
  fgsea_res <- fgsea(
    pathways = hallmark_pathways,
    stats    = cors,
    minSize  = 10,
    maxSize  = 500,
    nperm = 1000
  )

  fgsea_res <- fgsea_res[order(-NES)]

  n_up <- sum(fgsea_res$padj < 0.05 & fgsea_res$NES > 0, na.rm = TRUE)
  n_dn <- sum(fgsea_res$padj < 0.05 & fgsea_res$NES < 0, na.rm = TRUE)
  cat(sprintf("  Significant: %d UP, %d DOWN (FDR < 0.05)\n", n_up, n_dn))

  gsea_results[[g]] <- fgsea_res
  gsea_tables[[g]] <- as.data.frame(fgsea_res)[, c("pathway", "pval", "padj", "NES", "ES", "size")]

  write.table(gsea_tables[[g]],
              file = sprintf("SingleGene_GSEA/GSEA_table_%s.txt", g),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# ===========================================================================
# ---- 2.4.6 Generate Plots ----
# ===========================================================================
cat("\n═══ Generating Plots ═══\n")

# ─────────────────────────────────────────────────────────────────
# 5A. Running Enrichment Plots
# ─────────────────────────────────────────────────────────────────
for (g in core_genes) {
  dir.create(sprintf("SingleGene_GSEA/%s", g),
             showWarnings = FALSE, recursive = TRUE)

  res <- gsea_results[[g]]
  cors <- ranked_lists[[g]]

  res_up <- res[NES > 0][order(-NES)]
  res_dn <- res[NES < 0][order(NES)]

  plot_paths <- unique(c(head(res_up$pathway, 4), head(res_dn$pathway, 4)))
  plot_paths <- plot_paths[!is.na(plot_paths)]

  for (pw in plot_paths) {
    nes_val <- round(res[pathway == pw, NES], 2)
    pval_val <- res[pathway == pw, padj]
    pw_short <- gsub("HALLMARK_", "", pw)
    pw_short <- gsub("_", " ", pw_short)

    p <- plotEnrichment(hallmark_pathways[[pw]], cors) +
      ggtitle(paste(g, "—", pw_short)) +
      labs(subtitle = bquote(NES == .(nes_val) ~ "  FDR =" ~
                               .(format(pval_val, digits = 2, scientific = TRUE)))) +
      xlab("Gene Rank (by Spearman correlation)") +
      ylab("Running Enrichment Score") +
      theme_minimal(base_family = "serif") +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9)
      )

    safe_name <- gsub("[^A-Za-z0-9_]", "_", pw)
    ggsave(sprintf("SingleGene_GSEA/%s/GSEA_Running_%s.pdf", g, safe_name),
           p, width = 8, height = 5, device = "pdf")
    ggsave(sprintf("SingleGene_GSEA/%s/GSEA_Running_%s.png", g, safe_name),
           p, width = 8, height = 5, dpi = 150, device = "png")
  }
  cat(sprintf("  %s: %d running score plots\n", g, length(plot_paths)))
}

# ─────────────────────────────────────────────────────────────────
# 5B. Bar Plots (top 30 by |NES|)
# ─────────────────────────────────────────────────────────────────
for (g in core_genes) {
  res <- as.data.frame(gsea_results[[g]])
  res_plot <- res[order(-abs(res$NES)), TRUE]
  res_plot <- head(res_plot, 30)
  res_plot <- res_plot[order(res_plot$NES), TRUE]

  res_plot$pathway_short <- gsub("HALLMARK_", "", res_plot$pathway)
  res_plot$pathway_short <- gsub("_", " ", res_plot$pathway_short)
  res_plot$direction <- ifelse(res_plot$NES > 0,
                               "Positively correlated",
                               "Negatively correlated")
  res_plot$pathway_short <- factor(res_plot$pathway_short,
                                    levels = res_plot$pathway_short)
  res_plot$significant <- ifelse(res_plot$padj < 0.05, "FDR < 0.05", "ns")

  p_bar <- ggplot(res_plot, aes(x = NES, y = pathway_short,
                                 fill = direction, alpha = significant)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("Positively correlated" = "#B2182B",
                                  "Negatively correlated" = "#2166AC")) +
    scale_alpha_manual(values = c("FDR < 0.05" = 1, "ns" = 0.5)) +
    geom_text(aes(label = sprintf("NES=%.2f", NES),
                  x = ifelse(NES > 0, NES + 0.1, NES - 0.1),
                  hjust = ifelse(NES > 0, 0, 1)),
              size = 2.8, color = "grey30") +
    labs(x = "Normalized Enrichment Score (NES)", y = "",
         title = paste(g, "— Hallmark GSEA (Top 30 by |NES|)"),
         subtitle = "TCGA-PRAD (n=547) | Spearman correlation ranking",
         fill = "Direction", alpha = "Significance") +
    theme_minimal(base_family = "serif") +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )

  ggsave(sprintf("SingleGene_GSEA/%s/GSEA_BarPlot_%s.pdf", g, g),
         p_bar, width = 11, height = 10, device = "pdf")
  ggsave(sprintf("SingleGene_GSEA/%s/GSEA_BarPlot_%s.png", g, g),
         p_bar, width = 11, height = 10, dpi = 150, device = "png")

  cat(sprintf("  %s: bar plot (%d pathways)\n", g, nrow(res_plot)))
}

# ─────────────────────────────────────────────────────────────────
# 5C. NES Heatmap
# ─────────────────────────────────────────────────────────────────
cat("  Generating NES heatmap...\n")

all_pathways <- unique(unlist(lapply(gsea_results, function(x) x$pathway)))
nes_mat <- matrix(NA, length(all_pathways), length(core_genes),
                  dimnames = list(all_pathways, core_genes))
for (g in core_genes) {
  res_dt <- as.data.frame(gsea_results[[g]])
  rownames(res_dt) <- res_dt$pathway
  nes_mat[res_dt$pathway, g] <- res_dt$NES
}

nes_hm <- nes_mat
rn <- gsub("HALLMARK_", "", rownames(nes_hm))
rn <- gsub("_", " ", rn)
rownames(nes_hm) <- rn

max_abs <- max(abs(nes_hm), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs),
                       c("#2166AC", "#F7F7F7", "#B2182B"))

fdr_mat <- matrix(NA, nrow = nrow(nes_hm), ncol = ncol(nes_hm),
                  dimnames = dimnames(nes_hm))
for (g in core_genes) {
  res <- as.data.frame(gsea_results[[g]])
  rownames(res) <- res$pathway
  for (pw in all_pathways) {
    if (pw %in% res$pathway) {
      fdr_mat[gsub("HALLMARK_", "", pw), g] <- res[pw, "padj"]
    }
  }
}

ht <- Heatmap(
  nes_hm,
  name = "NES",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_names_max_width = unit(14, "cm"),
  width = unit(4, "cm"),
  height = unit(max(nrow(nes_hm) * 0.6, 8), "cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    fv <- fdr_mat[i, j]
    if (!is.na(fv) && fv < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 14, fontface = "bold", col = "white"))
    }
  },
  heatmap_legend_param = list(title = "NES", direction = "vertical"),
  column_title = "Single-Gene GSEA — Hallmark Pathways (50 gene sets)",
  column_title_gp = gpar(fontsize = 13, fontface = "bold")
)

pdf("SingleGene_GSEA/Figure_GSEA_NES_Heatmap.pdf", width = 8, height = 14)
draw(ht, heatmap_legend_side = "right")
dev.off()

png("SingleGene_GSEA/Figure_GSEA_NES_Heatmap.png",
    width = 8 * 150, height = 14 * 150, res = 150)
draw(ht, heatmap_legend_side = "right")
dev.off()
cat("  NES heatmap saved\n")

# ─────────────────────────────────────────────────────────────────
# 5D. Dot Plot
# ─────────────────────────────────────────────────────────────────
cat("  Generating dot plot...\n")

top_n_pw <- 10
dot_list <- lapply(core_genes, function(g) {
  res <- as.data.frame(gsea_results[[g]])
  res <- res[order(-abs(res$NES)), TRUE]
  top_up <- head(res[res$NES > 0, TRUE], top_n_pw)
  top_dn <- head(res[res$NES < 0, TRUE], top_n_pw)
  df <- rbind(top_up, top_dn)
  if (nrow(df) > 0) df$Gene <- g
  df
})

dot_df <- do.call(rbind, dot_list)

if (!is.null(dot_df) && nrow(dot_df) > 0) {
  dot_df <- dot_df[!is.na(dot_df$pathway), TRUE]
  if (nrow(dot_df) > 0) {
    dot_df$pathway_short <- gsub("HALLMARK_", "", dot_df$pathway)
    dot_df$pathway_short <- gsub("_", " ", dot_df$pathway_short)
    dot_df$neg_log10_padj <- -log10(dot_df$padj)
    max_finite <- suppressWarnings(
      max(dot_df$neg_log10_padj[is.finite(dot_df$neg_log10_padj)], na.rm = TRUE))
    if (!is.finite(max_finite)) max_finite <- 3
    dot_df$neg_log10_padj[is.infinite(dot_df$neg_log10_padj)] <- max_finite

    pw_order <- dot_df %>%
      group_by(pathway_short) %>%
      summarise(mean_abs_NES = mean(abs(NES), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_abs_NES))
    dot_df$pathway_short <- factor(dot_df$pathway_short,
                                    levels = rev(pw_order$pathway_short))

    p_dot <- ggplot(dot_df, aes(x = NES, y = pathway_short,
                                 color = Gene, size = neg_log10_padj)) +
      geom_point(alpha = 0.85) +
      scale_color_manual(values = pal_gene) +
      scale_size_continuous(name = "-log10(FDR)", range = c(2, 7)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
      labs(x = "Normalized Enrichment Score (NES)", y = "",
           title = "Single-Gene GSEA — Hallmark Pathways",
           subtitle = "TCGA-PRAD (n=547) | Top 10 by |NES| per direction") +
      theme_minimal(base_family = "serif") +
      theme(
        plot.title = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 11, face = "bold"),
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3)
      )

    ggsave("SingleGene_GSEA/Figure_GSEA_DotPlot.pdf", p_dot,
           width = 14, height = 16, device = "pdf")
    ggsave("SingleGene_GSEA/Figure_GSEA_DotPlot.png", p_dot,
           width = 14, height = 16, dpi = 150, device = "png")
    cat("  Dot plot saved\n")
  }
} else {
  cat("  Skipped dot plot (no data)\n")
}

# ─────────────────────────────────────────────────────────────────
# 5E. Multi-Gene Comparison Bar Chart
# ─────────────────────────────────────────────────────────────────
cat("  Generating multi-gene comparison...\n")

nes_all <- as.data.frame(nes_mat)
nes_all$pathway <- rownames(nes_all)
nes_long <- melt(nes_all, id.vars = "pathway",
                 variable.name = "Gene", value.name = "NES")
nes_long <- nes_long[!is.na(nes_long$NES), TRUE]
nes_long$pathway_short <- gsub("HALLMARK_", "", nes_long$pathway)
nes_long$pathway_short <- gsub("_", " ", nes_long$pathway_short)

pw_rank <- nes_long %>%
  group_by(pathway_short) %>%
  summarise(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_NES))

top_n_pw2 <- min(30, nrow(pw_rank))
top_pws <- pw_rank$pathway_short[1:top_n_pw2]
nes_top <- nes_long[nes_long$pathway_short %in% top_pws, TRUE]
nes_top$pathway_short <- factor(nes_top$pathway_short, levels = rev(top_pws))

p_multi <- ggplot(nes_top, aes(x = NES, y = pathway_short, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = pal_gene) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  labs(x = "Normalized Enrichment Score (NES)", y = "",
       title = "Single-Gene GSEA: Hallmark Pathway Enrichment Across Four Genes",
       subtitle = "TCGA-PRAD (n=547) | DNAH8 · TOP2A · CHRDL1 · PGM5") +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    axis.text.y = element_text(size = 8.5),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ggsave("SingleGene_GSEA/Figure_GSEA_MultiGene_Comparison.pdf",
       p_multi, width = 13, height = max(top_n_pw2 * 0.35, 6), device = "pdf")
ggsave("SingleGene_GSEA/Figure_GSEA_MultiGene_Comparison.png",
       p_multi, width = 13, height = max(top_n_pw2 * 0.35, 6), dpi = 150, device = "png")

cat("  Multi-gene comparison saved\n")

# ===========================================================================
# ---- 2.4.7 Summary ----
# ===========================================================================
cat("\n")
cat("═══════════════════════════════════════════════════\n")
cat("     Single-Gene GSEA — COMPLETE\n")
cat("═══════════════════════════════════════════════════\n\n")

for (g in core_genes) {
  res <- as.data.frame(gsea_results[[g]])
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)

  cat(sprintf("─── %s ───\n", g))
  cat(sprintf("  Significant (FDR<0.05): %d\n", n_sig))

  top3_up <- head(res[res$NES > 0, TRUE], 3)
  top3_dn <- head(res[res$NES < 0, TRUE], 3)
  if (nrow(top3_up) > 0) {
    cat(sprintf("  Top UP:   %s (NES: %s)\n",
                paste(gsub("HALLMARK_", "", top3_up$pathway), collapse = ", "),
                paste(round(top3_up$NES, 2), collapse = ", ")))
  }
  if (nrow(top3_dn) > 0) {
    cat(sprintf("  Top DOWN: %s (NES: %s)\n",
                paste(gsub("HALLMARK_", "", top3_dn$pathway), collapse = ", "),
                paste(round(top3_dn$NES, 2), collapse = ", ")))
  }
  cat("\n")
}

cat(sprintf("Results saved to: %s/SingleGene_GSEA/\n", getwd()))
cat("═══════════════════════════════════════════════════\n")
