if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# ==== 0.4 计算最终预测的AUC区间 ====
# Execution scope: Entry point; prerequisites remain unchanged
# See CHAPTERS.md for order, dependencies and entry points.

# Run ONLY after the final cohort/model and cluster unit have been resolved.
# Create data/external/final_predictions.csv with the listed columns.
# This calculates fixed-score uncertainty, not full-pipeline validation.
args<-commandArgs(trailingOnly=TRUE)
if(!"--confirmed-final-input" %in% args)stop("Provide reviewed final_predictions.csv and use --confirmed-final-input; do not treat this as full-pipeline validation.")
cmd<-commandArgs(trailingOnly=FALSE);script<-sub("^--file=","",cmd[grepl("^--file=",cmd)])
CODE_ROOT<-dirname(normalizePath(script,winslash="/",mustWork=TRUE))
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"));for(m in c("io","model_audit","auc_uncertainty"))source(module_path(m))
p<-project_path("data","external","final_predictions.csv")
if(!file.exists(p))stop("Missing reviewed input: ",p)
d<-read.csv(p,stringsAsFactors=FALSE)
assert_columns(d,c("cohort","sample_id","patient_id","BCR","prediction","prediction_type","model_id"))
assert_unique(d,c("cohort","sample_id","prediction_type","model_id"))
keys<-unique(d[c("cohort","prediction_type","model_id")]);result<-list()
for(i in seq_len(nrow(keys))) {
 k<-keys[i,];x<-d[d$cohort==k$cohort & d$prediction_type==k$prediction_type & d$model_id==k$model_id,]
 result[[i]]<-cbind(k,cluster_bootstrap_auc(x$BCR,x$prediction,x$patient_id))
}
dir.create(project_path("outputs","revision"),recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,result),project_path("outputs","revision","auc_cluster_bootstrap.csv"),row.names=FALSE)
