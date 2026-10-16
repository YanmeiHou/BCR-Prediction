out <- project_path('outputs','original_pooled_dca')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source_data <- data_path('historical_dca','prepared_data.RDS')
prep <- readRDS(source_data)
expr <- prep$expr
label <- prep$label
stopifnot(nrow(expr) == nrow(label), identical(rownames(expr), rownames(label)))
all_scaled <- as.data.frame(scale(expr))
all_scaled$BCR <- as.numeric(as.character(label$outcome))
all_scaled$Cohort <- label$Cohort
stopifnot(!anyNA(all_scaled), all(all_scaled$BCR %in% c(0, 1)))
# Reproduce the original Combined model; all four cohorts participated in fitting.
fit <- glm(BCR ~ DNAH8 + TOP2A + CHRDL1 + PGM5,
           data = all_scaled, family = binomial)
pp <- predict(fit, newdata = all_scaled, type = 'response')
write.csv(data.frame(Sample_ID = rownames(expr), Cohort = label$Cohort,
                    BCR = all_scaled$BCR, Predicted_probability = pp),
          file.path(out, 'Original_DCA_predictions.csv'), row.names = FALSE)
thresholds <- seq(0.01, 0.50, by = 0.01)
cohorts <- c('TCGA', 'DKFZ', 'Taylor', 'CancerMap')
details <- list(); summaries <- list()
interval_text <- function(ok) {
  inds <- which(ok)
  if (!length(inds)) return('None')
  blocks <- split(inds, cumsum(c(TRUE, diff(inds) != 1)))
  paste(vapply(blocks, function(b) {
    v <- round(thresholds[b] * 100)
    if (length(v) == 1) paste0(v, '%') else paste0(min(v), '%-', max(v), '%')
  }, character(1)), collapse = '; ')
}
for (cn in cohorts) {
  idx <- which(label$Cohort == cn)
  pred <- pp[idx]; y <- all_scaled$BCR[idx]; n <- length(y)
  tp <- vapply(thresholds, function(t) sum(pred >= t & y == 1), integer(1))
  fp <- vapply(thresholds, function(t) sum(pred >= t & y == 0), integer(1))
  model <- tp/n - (fp/n) * thresholds/(1-thresholds)
  prev <- mean(y)
  treat_all <- prev - (1-prev) * thresholds/(1-thresholds)
  # Exclude floating-point noise at theoretical equalities.
  good <- model > 1e-12 & model - treat_all > 1e-12
  details[[cn]] <- data.frame(Cohort = cn, N = n, Events = sum(y),
    Threshold = thresholds, TP = tp, FP = fp, Net_benefit_model = model,
    Net_benefit_treat_all = treat_all, Net_benefit_treat_none = 0,
    Model_minus_all = model-treat_all, Benefit_over_both = good,
    Strict_floating_comparison = model > 0 & model > treat_all)
  summaries[[cn]] <- data.frame(Cohort = cn, N = n, Events = sum(y),
    Benefit_intervals = interval_text(good), Beneficial_grid_points = sum(good),
    Evaluated_range = '1%-50%', Step = '1%',
    Evaluation = 'Apparent performance: four cohorts pooled for model fitting')
}
detail <- do.call(rbind, details); summary_table <- do.call(rbind, summaries)
write.csv(detail, file.path(out, 'Original_DCA_threshold_data.csv'), row.names = FALSE)
write.csv(summary_table, file.path(out, 'Original_DCA_benefit_intervals.csv'), row.names = FALSE)
write.csv(data.frame(Term = names(coef(fit)), Coefficient = unname(coef(fit))),
          file.path(out, 'Original_DCA_model_coefficients.csv'), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
print(summary_table, row.names = FALSE)
cat('Floating-point equality exclusions:', sum(detail$Strict_floating_comparison != detail$Benefit_over_both), '\n')
