# ==== 0.3 样本数量与患者重复核对 ====
# Execution scope: 已复算
# See CHAPTERS.md for order, dependencies and entry points.

# Describe the supplied sample set; do not silently redefine inclusion or deduplicate.
# ---- 0.3.1 audit_samples ----
audit_samples <- function(samples) {
  assert_columns(samples,c("cohort","sample_id","model_sample_id","patient_id","BCR","pcadb_group","primary_tumor_sensitivity"))
  assert_unique(samples,c("cohort","sample_id")); assert_unique(samples,"model_sample_id")
  stopifnot(all(samples$BCR %in% c(0,1)), !anyNA(samples$primary_tumor_sensitivity))
  result <- list()
  for (population in c("original_1025","primary_tumor_only")) {
    x <- if(population=="original_1025") samples else samples[samples$primary_tumor_sensitivity,]
    for(cohort in unique(x$cohort)) {
      a <- x[x$cohort==cohort,]; n <- table(a$patient_id)
      result[[length(result)+1]] <- data.frame(population,cohort,n_samples=nrow(a),BCR=sum(a$BCR==1),nonBCR=sum(a$BCR==0),
        n_patients=length(n),repeated_patients=sum(n>1),extra_samples=sum(pmax(n-1,0)),
        tissue_types=paste(names(table(a$pcadb_group)),as.integer(table(a$pcadb_group)),sep=":",collapse=";"))
    }
  }
  do.call(rbind,result)
}
