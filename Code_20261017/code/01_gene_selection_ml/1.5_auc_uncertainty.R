# ==== 1.5 AUC患者簇Bootstrap置信区间 ====
# Execution scope: 工具已检查；需最终预测表
# See CHAPTERS.md for order, dependencies and entry points.

# Optional future calculation. Not run by run_reproduction.R.
# Resamples patient clusters within each cohort, holding supplied scores fixed.
# Does NOT repeat discovery/model selection and does NOT resolve tissue inclusion.
# ---- 1.5.1 cluster_bootstrap_auc ----
cluster_bootstrap_auc <- function(y,score,patient_id,B=2000L,seed=777L) {
  if(anyNA(patient_id)||any(!nzchar(patient_id)))stop("Patient IDs are required")
  if(length(y)!=length(patient_id)||length(score)!=length(y))stop("Mismatched input lengths")
  if(B<200L)stop("Use at least 200 resamples for diagnostic use; at least 2000 for a finalized analysis")
  groups<-split(seq_along(y),patient_id);ids<-names(groups)
  set.seed(seed)
  values<-rep(NA_real_,B)
  for(b in seq_len(B)) {
    selected<-sample(ids,length(ids),replace=TRUE)
    idx<-unlist(groups[selected],use.names=FALSE)
    values[b]<-auc_fixed(y[idx],score[idx])
  }
  valid<-is.finite(values)
  if(sum(valid)<.95*B)stop("Too many bootstrap samples have an undefined AUC")
  interval<-quantile(values[valid],c(.025,.975),names=FALSE)
  data.frame(AUC=auc_fixed(y,score),CI95_low=interval[1],CI95_high=interval[2],valid_resamples=sum(valid),requested_resamples=B,
    method="patient-cluster percentile bootstrap; fixed supplied scores",seed=seed)
}
