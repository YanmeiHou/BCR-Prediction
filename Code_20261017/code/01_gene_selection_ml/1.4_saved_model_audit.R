# ==== 1.4 现有预测AUC、模型记录与筛选交集核对 ====
# Execution scope: 已复算
# See CHAPTERS.md for order, dependencies and entry points.

# AUC is calculated with a fixed direction: higher score means BCR=1.
# This reproduces supplied probabilities; it does not train or select a model.
# ---- 1.4.1 auc_fixed ----
auc_fixed <- function(y, score) {
  if(length(y)!=length(score) || anyNA(y) || any(!is.finite(score)) || !all(y %in% c(0,1))) stop("Invalid AUC input")
  n1 <- sum(y==1); n0 <- sum(y==0)
  if(!n1 || !n0) return(NA_real_)
  (sum(rank(score,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}
# ---- 1.4.2 evaluate_predictions ----
evaluate_predictions <- function(predictions,samples) {
  assert_unique(predictions,c("sample_id","prediction_type"))
  idx <- match(predictions$sample_id,samples$model_sample_id)
  stopifnot(!anyNA(idx),all(predictions$cohort==samples$cohort[idx]),all(predictions$BCR==samples$BCR[idx]))
  if(any(!is.finite(predictions$prediction)) || any(predictions$prediction<0 | predictions$prediction>1)) stop("Invalid probability")
  keys <- unique(predictions[c("cohort","prediction_type","model_id")]); result <- list()
  for(i in seq_len(nrow(keys))) {
    k <- keys[i,]; d <- predictions[predictions$cohort==k$cohort & predictions$prediction_type==k$prediction_type & predictions$model_id==k$model_id,]
    result[[i]] <- data.frame(k,n_samples=nrow(d),BCR=sum(d$BCR==1),nonBCR=sum(d$BCR==0),AUC=auc_fixed(d$BCR,d$prediction))
  }
  do.call(rbind,result)
}
# ---- 1.4.3 audit_model_records ----
audit_model_records <- function() {
  config <- read_data("algorithm_configuration.csv"); models <- read_data("saved_models.csv")
  auc <- read_data("model_auc.csv"); ranking <- read_data("saved_ranking.csv")
  assert_unique(config,"model_id");assert_unique(models,"model_id")
  auc$mean_including_training <- rowMeans(auc[c("TCGA","DKFZ","Taylor","CancerMap")])
  auc$external_descriptive_mean <- rowMeans(auc[c("DKFZ","Taylor","CancerMap")])
  idx <- match(auc$model_id,ranking$model_id)
  stopifnot(!anyNA(idx),max(abs(auc$mean_including_training-ranking$saved_mean_AUC[idx]))<1e-10)
  write_result(auc,"model_auc_audit.csv")
  write_result(config[!config$model_id %in% models$model_id,],"unsaved_configurations.csv")
  sets <- read_data("feature_sets.csv")
  # Preserve original feature identifiers. Author-reported relabeling is archived,
  # not used as evidence that the four selected genes were recovered independently.
  original <- read.csv(project_path("data","reference","feature_label_originals","feature_sets.csv"),stringsAsFactors=FALSE)
  stopifnot(identical(sets,original))
  getset <- function(label) sets$gene[sets$set==label]
  intersection <- intersect(getset("rf_top20"),getset("lasso_min"))
  write_result(data.frame(gene=intersection),"recomputed_RF_LASSO_intersection.csv")
  stopifnot(setequal(intersection,getset("intersection_min")))
}
