# ==== 2.10 合并队列列线图与DCA ====
# Execution scope: 历史备选版本
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
# Combined + Individual Nomograms — 4-Gene BCR Prediction
# Output: Nomogram/
################################################################################

library(rms)
library(pROC)

setwd(legacy_output())
prep <- readRDS("prepared_data.RDS")
expr   <- prep$expr
label  <- prep$label
genes  <- prep$genes

cohorts <- c("TCGA", "DKFZ", "Taylor", "CancerMap")
out_dir <- "Nomogram"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# PART 1: COMBINED NOMOGRAM (all 4 datasets)
# ============================================================================
cat("========== Combined Nomogram ==========\n")

all_scaled <- as.data.frame(scale(expr))
all_scaled$BCR <- as.numeric(as.character(label$outcome))
all_scaled$Cohort <- label$Cohort

comb_model <- glm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                  data = all_scaled, family = binomial)
cat("Combined model:\n")
print(summary(comb_model)$coefficients)

dd_comb <- datadist(all_scaled)
options(datadist = "dd_comb")

lrm_comb <- lrm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                data = all_scaled, x = TRUE, y = TRUE)

nom_comb <- nomogram(lrm_comb,
                     fun = function(x) 1/(1+exp(-x)),
                     fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                     funlabel = "Predicted BCR Risk", lp = FALSE)

pdf(file.path(out_dir, "Fig_Combined_Nomogram.pdf"), width = 10, height = 7)
par(mar = c(4, 4, 3, 2))
plot(nom_comb, lplabel = "Linear Predictor", xfrac = 0.35,
     cex.axis = 0.9, cex.var = 1.1, col.grid = c("grey80", "grey90"))
title(main = "BCR Prediction Nomogram — 4 Cohorts Combined", cex.main = 1.2, font = 2)
dev.off()

png(file.path(out_dir, "Fig_Combined_Nomogram.png"), width = 10*150, height = 7*150, res = 150)
par(mar = c(4, 4, 3, 2))
plot(nom_comb, lplabel = "Linear Predictor", xfrac = 0.35,
     cex.axis = 0.9, cex.var = 1.1, col.grid = c("grey80", "grey90"))
title(main = "BCR Prediction Nomogram — 4 Cohorts Combined", cex.main = 1.2, font = 2)
dev.off()
cat("Combined nomogram saved.\n")

# Predictions
pred_comb <- predict(comb_model, newdata = all_scaled, type = "response")
true_y    <- as.numeric(as.character(label$outcome))

# ---- 2.10.1 Combined Calibration ----
pdf(file.path(out_dir, "Fig_Combined_Calibration.pdf"), width = 11, height = 9)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  pp <- pred_comb[idx]
  ty <- true_y[idx]

  n_grp <- 10
  brks <- quantile(pp, probs = seq(0, 1, length.out = n_grp + 1), na.rm = TRUE)
  grp <- cut(pp, breaks = brks, include.lowest = TRUE)

  obs_mean  <- tapply(ty, grp, mean)
  pred_mean <- tapply(pp, grp, mean)
  n_grp_obs <- tapply(ty, grp, length)

  plot(pred_mean, obs_mean,
       xlim = c(0, max(c(pred_mean, obs_mean), na.rm = TRUE) * 1.2),
       ylim = c(0, max(c(pred_mean, obs_mean), na.rm = TRUE) * 1.2),
       xlab = "Predicted Probability", ylab = "Observed Proportion",
       main = paste0("Calibration (Combined) — ", cn),
       pch = 19, col = "#377EB8", cex = n_grp_obs / max(n_grp_obs) * 3 + 1.5,
       cex.main = 1.2, cex.lab = 1.1, las = 1)
  abline(0, 1, lty = 2, col = "grey50", lwd = 2)

  ok <- !is.na(pred_mean) & !is.na(obs_mean)
  if (sum(ok) >= 3) {
    lo <- loess(obs_mean[ok] ~ pred_mean[ok], span = 1.5)
    px <- seq(min(pred_mean, na.rm = TRUE), max(pred_mean, na.rm = TRUE), length = 50)
    py <- predict(lo, newdata = data.frame(pred_mean = px))
    lines(px, py, col = "#D64045", lwd = 3)
  }

  brier <- mean((pp - ty)^2)
  roc_obj <- roc(ty, pp, quiet = TRUE)
  legend("bottomright",
         c(paste0("n = ", length(idx)),
           paste0("AUC = ", round(auc(roc_obj), 3)),
           paste0("Brier = ", round(brier, 4)),
           "Ideal", "Loess"),
         col = c("#377EB8", NA, NA, "grey50", "#D64045"),
         lty = c(NA, NA, NA, 2, 1), lwd = c(NA, NA, NA, 2, 3),
         pch = c(19, NA, NA, NA, NA), bty = "n", cex = 0.85)
}
dev.off()

png(file.path(out_dir, "Fig_Combined_Calibration.png"), width = 11*150, height = 9*150, res = 150)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))
for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  pp <- pred_comb[idx]; ty <- true_y[idx]
  n_grp <- 10
  brks <- quantile(pp, probs = seq(0, 1, length.out = n_grp + 1), na.rm = TRUE)
  grp <- cut(pp, breaks = brks, include.lowest = TRUE)
  obs_mean  <- tapply(ty, grp, mean)
  pred_mean <- tapply(pp, grp, mean)
  n_grp_obs <- tapply(ty, grp, length)
  plot(pred_mean, obs_mean,
       xlim = c(0, max(c(pred_mean, obs_mean), na.rm = TRUE) * 1.2),
       ylim = c(0, max(c(pred_mean, obs_mean), na.rm = TRUE) * 1.2),
       xlab = "Predicted Probability", ylab = "Observed Proportion",
       main = paste0("Calibration (Combined) — ", cn),
       pch = 19, col = "#377EB8", cex = n_grp_obs / max(n_grp_obs) * 3 + 1.5,
       cex.main = 1.2, cex.lab = 1.1, las = 1)
  abline(0, 1, lty = 2, col = "grey50", lwd = 2)
  ok <- !is.na(pred_mean) & !is.na(obs_mean)
  if (sum(ok) >= 3) {
    lo <- loess(obs_mean[ok] ~ pred_mean[ok], span = 1.5)
    px <- seq(min(pred_mean, na.rm = TRUE), max(pred_mean, na.rm = TRUE), length = 50)
    py <- predict(lo, newdata = data.frame(pred_mean = px))
    lines(px, py, col = "#D64045", lwd = 3)
  }
  brier <- mean((pp - ty)^2)
  roc_obj <- roc(ty, pp, quiet = TRUE)
  legend("bottomright",
         c(paste0("n = ", length(idx)),
           paste0("AUC = ", round(auc(roc_obj), 3)),
           paste0("Brier = ", round(brier, 4)),
           "Ideal", "Loess"),
         col = c("#377EB8", NA, NA, "grey50", "#D64045"),
         lty = c(NA, NA, NA, 2, 1), lwd = c(NA, NA, NA, 2, 3),
         pch = c(19, NA, NA, NA, NA), bty = "n", cex = 0.85)
}
dev.off()
cat("Combined calibration saved.\n")

# ---- 2.10.2 Combined DCA ----
# ---- 2.10.3 compute_nb ----
compute_nb <- function(pp, ty, th) {
  sapply(th, function(t) {
    tp <- sum(pp >= t & ty == 1); fp <- sum(pp >= t & ty == 0); n <- length(ty)
    if (n == 0) return(NA)
    tp / n - (fp / n) * (t / (1 - t))
  })
}

thresholds <- seq(0.01, 0.50, by = 0.01)

pdf(file.path(out_dir, "Fig_Combined_DCA.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  pp <- pred_comb[idx]; ty <- true_y[idx]

  nb_model <- compute_nb(pp, ty, thresholds)
  nb_all   <- sapply(thresholds, function(t) {
    prev <- sum(ty == 1) / length(ty); prev - (1 - prev) * (t / (1 - t))
  })

  ymax <- max(c(nb_model, nb_all, 0), na.rm = TRUE)
  ymin <- min(c(nb_model, nb_all, 0), na.rm = TRUE)

  plot(thresholds, nb_model, type = "l", lwd = 3, col = "#D64045",
       xlab = "Threshold Probability", ylab = "Net Benefit",
       main = paste0("DCA (Combined) — ", cn),
       ylim = c(ymin, ymax * 1.1), cex.main = 1.2, cex.lab = 1.1, las = 1)
  lines(thresholds, nb_all, lwd = 2.5, col = "#377EB8", lty = 2)
  abline(h = 0, lty = 3, col = "grey50", lwd = 1.5)
  legend("topright", legend = c("4-Gene Model", "Treat All", "Treat None"),
         col = c("#D64045", "#377EB8", "grey50"),
         lty = c(1, 2, 3), lwd = c(3, 2.5, 1.5), bty = "n", cex = 0.9)
}
dev.off()

png(file.path(out_dir, "Fig_Combined_DCA.png"), width = 12*150, height = 10*150, res = 150)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))
for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  pp <- pred_comb[idx]; ty <- true_y[idx]
  nb_model <- compute_nb(pp, ty, thresholds)
  nb_all   <- sapply(thresholds, function(t) {
    prev <- sum(ty == 1) / length(ty); prev - (1 - prev) * (t / (1 - t))
  })
  ymax <- max(c(nb_model, nb_all, 0), na.rm = TRUE)
  ymin <- min(c(nb_model, nb_all, 0), na.rm = TRUE)
  plot(thresholds, nb_model, type = "l", lwd = 3, col = "#D64045",
       xlab = "Threshold Probability", ylab = "Net Benefit",
       main = paste0("DCA (Combined) — ", cn),
       ylim = c(ymin, ymax * 1.1), cex.main = 1.2, cex.lab = 1.1, las = 1)
  lines(thresholds, nb_all, lwd = 2.5, col = "#377EB8", lty = 2)
  abline(h = 0, lty = 3, col = "grey50", lwd = 1.5)
  legend("topright", legend = c("4-Gene Model", "Treat All", "Treat None"),
         col = c("#D64045", "#377EB8", "grey50"),
         lty = c(1, 2, 3), lwd = c(3, 2.5, 1.5), bty = "n", cex = 0.9)
}
dev.off()
cat("Combined DCA saved.\n\n")

# ============================================================================
# PART 2: INDIVIDUAL NOMOGRAMS (one per dataset)
# ============================================================================
cat("========== Individual Nomograms ==========\n")

individual_auc <- c()

# 4-panel nomogram
pdf(file.path(out_dir, "Fig_Individual_Nomograms.pdf"), width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  n <- length(idx)

  ds_scaled <- as.data.frame(scale(expr[idx, , drop = FALSE]))
  ds_scaled$BCR <- as.numeric(as.character(label$outcome[idx]))

  ds_model <- glm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                  data = ds_scaled, family = binomial)

  dd_ds <- datadist(ds_scaled)
  options(datadist = "dd_ds")

  lrm_ds <- lrm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                data = ds_scaled, x = TRUE, y = TRUE)

  nom_ds <- nomogram(lrm_ds,
                     fun = function(x) 1/(1+exp(-x)),
                     fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                     funlabel = "Predicted BCR Risk", lp = FALSE)

  plot(nom_ds, lplabel = "Linear Predictor", xfrac = 0.35,
       cex.axis = 0.75, cex.var = 0.9, col.grid = c("grey80", "grey90"))
  title(main = paste0("Nomogram — ", cn, " (n=", n, ")"), cex.main = 1.2, font = 2)

  pred_ds <- predict(ds_model, newdata = ds_scaled, type = "response")
  roc_obj <- roc(ds_scaled$BCR, pred_ds, quiet = TRUE)
  individual_auc[cn] <- auc(roc_obj)
}
dev.off()

png(file.path(out_dir, "Fig_Individual_Nomograms.png"), width = 16*150, height = 12*150, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  n <- length(idx)
  ds_scaled <- as.data.frame(scale(expr[idx, , drop = FALSE]))
  ds_scaled$BCR <- as.numeric(as.character(label$outcome[idx]))
  ds_model <- glm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                  data = ds_scaled, family = binomial)
  dd_ds <- datadist(ds_scaled)
  options(datadist = "dd_ds")
  lrm_ds <- lrm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                data = ds_scaled, x = TRUE, y = TRUE)
  nom_ds <- nomogram(lrm_ds,
                     fun = function(x) 1/(1+exp(-x)),
                     fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                     funlabel = "Predicted BCR Risk", lp = FALSE)
  plot(nom_ds, lplabel = "Linear Predictor", xfrac = 0.35,
       cex.axis = 0.75, cex.var = 0.9, col.grid = c("grey80", "grey90"))
  title(main = paste0("Nomogram — ", cn, " (n=", n, ")"), cex.main = 1.2, font = 2)
}
dev.off()
cat("Individual nomograms saved.\n")

# ---- 2.10.4 Individual Calibration ----
pdf(file.path(out_dir, "Fig_Individual_Calibration.pdf"), width = 11, height = 9)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  ds_scaled <- as.data.frame(scale(expr[idx, , drop = FALSE]))
  ds_scaled$BCR <- as.numeric(as.character(label$outcome[idx]))

  ds_model <- glm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                  data = ds_scaled, family = binomial)
  pp <- predict(ds_model, newdata = ds_scaled, type = "response")
  ty <- ds_scaled$BCR

  n_grp <- 10
  brks <- quantile(pp, probs = seq(0, 1, length.out = n_grp + 1), na.rm = TRUE)
  if (length(unique(brks)) < 3) next
  grp <- cut(pp, breaks = brks, include.lowest = TRUE)

  obs_mean  <- tapply(ty, grp, mean)
  pred_mean <- tapply(pp, grp, mean)
  n_grp_obs <- tapply(ty, grp, length)

  plot(pred_mean, obs_mean,
       xlim = c(0, max(c(pred_mean, obs_mean), na.rm = TRUE) * 1.2),
       ylim = c(0, max(c(pred_mean, obs_mean), na.rm = TRUE) * 1.2),
       xlab = "Predicted Probability", ylab = "Observed Proportion",
       main = paste0("Calibration (Own Model) — ", cn),
       pch = 19, col = "#4DAF4A", cex = n_grp_obs / max(n_grp_obs) * 3 + 1.5,
       cex.main = 1.2, cex.lab = 1.1, las = 1)
  abline(0, 1, lty = 2, col = "grey50", lwd = 2)

  ok <- !is.na(pred_mean) & !is.na(obs_mean)
  if (sum(ok) >= 3) {
    lo <- loess(obs_mean[ok] ~ pred_mean[ok], span = 1.5)
    px <- seq(min(pred_mean, na.rm = TRUE), max(pred_mean, na.rm = TRUE), length = 50)
    py <- predict(lo, newdata = data.frame(pred_mean = px))
    lines(px, py, col = "#984EA3", lwd = 3)
  }

  brier <- mean((pp - ty)^2)
  legend("bottomright",
         c(paste0("n = ", length(idx)),
           paste0("AUC = ", round(individual_auc[cn], 3)),
           paste0("Brier = ", round(brier, 4))),
         bty = "n", cex = 0.85)
}
dev.off()

png(file.path(out_dir, "Fig_Individual_Calibration.png"), width = 11*150, height = 9*150, res = 150)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))
for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  ds_scaled <- as.data.frame(scale(expr[idx, , drop = FALSE]))
  ds_scaled$BCR <- as.numeric(as.character(label$outcome[idx]))
  ds_model <- glm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                  data = ds_scaled, family = binomial)
  pp <- predict(ds_model, newdata = ds_scaled, type = "response")
  ty <- ds_scaled$BCR
  n_grp <- 10
  brks <- quantile(pp, probs = seq(0, 1, length.out = n_grp + 1), na.rm = TRUE)
  if (length(unique(brks)) < 3) next
  grp <- cut(pp, breaks = brks, include.lowest = TRUE)
  obs_mean  <- tapply(ty, grp, mean)
  pred_mean <- tapply(pp, grp, mean)
  n_grp_obs <- tapply(ty, grp, length)
  plot(pred_mean, obs_mean,
       xlim = c(0, max(c(pred_mean, obs_mean), na.rm = TRUE) * 1.2),
       ylim = c(0, max(c(pred_mean, obs_mean), na.rm = TRUE) * 1.2),
       xlab = "Predicted Probability", ylab = "Observed Proportion",
       main = paste0("Calibration (Own Model) — ", cn),
       pch = 19, col = "#4DAF4A", cex = n_grp_obs / max(n_grp_obs) * 3 + 1.5,
       cex.main = 1.2, cex.lab = 1.1, las = 1)
  abline(0, 1, lty = 2, col = "grey50", lwd = 2)
  ok <- !is.na(pred_mean) & !is.na(obs_mean)
  if (sum(ok) >= 3) {
    lo <- loess(obs_mean[ok] ~ pred_mean[ok], span = 1.5)
    px <- seq(min(pred_mean, na.rm = TRUE), max(pred_mean, na.rm = TRUE), length = 50)
    py <- predict(lo, newdata = data.frame(pred_mean = px))
    lines(px, py, col = "#984EA3", lwd = 3)
  }
  brier <- mean((pp - ty)^2)
  legend("bottomright",
         c(paste0("n = ", length(idx)),
           paste0("AUC = ", round(individual_auc[cn], 3)),
           paste0("Brier = ", round(brier, 4))),
         bty = "n", cex = 0.85)
}
dev.off()
cat("Individual calibration saved.\n")

# ---- 2.10.5 Individual DCA ----
pdf(file.path(out_dir, "Fig_Individual_DCA.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  ds_scaled <- as.data.frame(scale(expr[idx, , drop = FALSE]))
  ds_scaled$BCR <- as.numeric(as.character(label$outcome[idx]))
  ds_model <- glm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                  data = ds_scaled, family = binomial)
  pp <- predict(ds_model, newdata = ds_scaled, type = "response")
  ty <- ds_scaled$BCR

  nb_model <- compute_nb(pp, ty, thresholds)
  nb_all   <- sapply(thresholds, function(t) {
    prev <- sum(ty == 1) / length(ty); prev - (1 - prev) * (t / (1 - t))
  })

  ymax <- max(c(nb_model, nb_all, 0), na.rm = TRUE)
  ymin <- min(c(nb_model, nb_all, 0), na.rm = TRUE)

  plot(thresholds, nb_model, type = "l", lwd = 3, col = "#4DAF4A",
       xlab = "Threshold Probability", ylab = "Net Benefit",
       main = paste0("DCA (Own Model) — ", cn),
       ylim = c(ymin, ymax * 1.1), cex.main = 1.2, cex.lab = 1.1, las = 1)
  lines(thresholds, nb_all, lwd = 2.5, col = "#377EB8", lty = 2)
  abline(h = 0, lty = 3, col = "grey50", lwd = 1.5)
  legend("topright", legend = c("4-Gene Model", "Treat All", "Treat None"),
         col = c("#4DAF4A", "#377EB8", "grey50"),
         lty = c(1, 2, 3), lwd = c(3, 2.5, 1.5), bty = "n", cex = 0.9)
}
dev.off()

png(file.path(out_dir, "Fig_Individual_DCA.png"), width = 12*150, height = 10*150, res = 150)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))
for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  ds_scaled <- as.data.frame(scale(expr[idx, , drop = FALSE]))
  ds_scaled$BCR <- as.numeric(as.character(label$outcome[idx]))
  ds_model <- glm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
                  data = ds_scaled, family = binomial)
  pp <- predict(ds_model, newdata = ds_scaled, type = "response")
  ty <- ds_scaled$BCR
  nb_model <- compute_nb(pp, ty, thresholds)
  nb_all   <- sapply(thresholds, function(t) {
    prev <- sum(ty == 1) / length(ty); prev - (1 - prev) * (t / (1 - t))
  })
  ymax <- max(c(nb_model, nb_all, 0), na.rm = TRUE)
  ymin <- min(c(nb_model, nb_all, 0), na.rm = TRUE)
  plot(thresholds, nb_model, type = "l", lwd = 3, col = "#4DAF4A",
       xlab = "Threshold Probability", ylab = "Net Benefit",
       main = paste0("DCA (Own Model) — ", cn),
       ylim = c(ymin, ymax * 1.1), cex.main = 1.2, cex.lab = 1.1, las = 1)
  lines(thresholds, nb_all, lwd = 2.5, col = "#377EB8", lty = 2)
  abline(h = 0, lty = 3, col = "grey50", lwd = 1.5)
  legend("topright", legend = c("4-Gene Model", "Treat All", "Treat None"),
         col = c("#4DAF4A", "#377EB8", "grey50"),
         lty = c(1, 2, 3), lwd = c(3, 2.5, 1.5), bty = "n", cex = 0.9)
}
dev.off()
cat("Individual DCA saved.\n\n")

# ============================================================================
# Summary
# ============================================================================
cat("========== Summary ==========\n")
cat("\nCombined Model AUC:\n")
for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  roc_obj <- roc(true_y[idx], pred_comb[idx], quiet = TRUE)
  cat(sprintf("  %-12s AUC = %.3f\n", cn, auc(roc_obj)))
}

cat("\nIndividual Model AUC:\n")
for (cn in cohorts) {
  cat(sprintf("  %-12s AUC = %.3f\n", cn, individual_auc[cn]))
}

cat("\nDone! All files in Nomogram/\n")
