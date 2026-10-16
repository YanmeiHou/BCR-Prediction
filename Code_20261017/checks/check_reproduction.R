if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# ==== 8.1 现有数据复算结果检查 ====
# Execution scope: 检查脚本
# See CHAPTERS.md for order, dependencies and entry points.

# After Rscript run_reproduction.R, run: Rscript checks/check_reproduction.R
CODE_ROOT<-normalizePath(".",winslash="/",mustWork=TRUE)
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"));source(project_path("code/00_setup/0.2_data_io.R"));source(project_path("code/01_gene_selection_ml/1.4_saved_model_audit.R"))
# ---- 8.1.1 read_output ----
read_output<-function(name)read.csv(output_path(name),check.names=FALSE)
# ---- 8.1.2 read_ref ----
read_ref<-function(name)read.csv(project_path("data","reference",name),check.names=FALSE)
# ---- 8.1.3 same_numeric ----
same_numeric<-function(actual,expected,tolerance=1e-9){stopifnot(identical(is.na(actual),is.na(expected)));ok<-is.finite(expected);stopifnot(max(abs(actual[ok]-expected[ok]),0)<tolerance)}
ref<-read_ref("auc.csv");x<-read_output("auc.csv")
k<-match(paste(ref$cohort,ref$prediction_type),paste(x$cohort,x$prediction_type));stopifnot(!anyNA(k))
same_numeric(x$AUC[k],ref$recomputed_AUC)
refs<-rbind(read_ref("immune_pooled.csv"),read_ref("immune_by_cohort.csv"));x<-read_output("immune_tests.csv")
# ---- 8.1.4 key ----
key<-function(z)paste(z$population,z$cohort,z$cell_type,sep="|")
k<-match(key(refs),key(x));stopifnot(!anyNA(k),nrow(x)==220)
for(col in c("welch_p_raw","wilcoxon_p_raw","welch_p_BH22","wilcoxon_p_BH22"))same_numeric(x[[col]][k],refs[[col]])
refs<-rbind(read_ref("immune_meta_cohen.csv"),read_ref("immune_meta_hedges.csv"));x<-read_output("immune_meta.csv")
# ---- 8.1.5 key ----
key<-function(z)paste(z$population,z$effect_type,z$cell_type,sep="|")
k<-match(key(refs),key(x));stopifnot(!anyNA(k),nrow(x)==88)
for(col in c("pooled_effect","CI95_low","CI95_high","p_raw","p_BH22","tau2","I2_percent"))same_numeric(x[[col]][k],refs[[col]])
# AUC direction, ties and missing-class behavior have deliberately explicit tests.
stopifnot(auc_fixed(c(0,0,1,1),c(.1,.2,.8,.9))==1,auc_fixed(c(0,1),c(.5,.5))==.5,
          auc_fixed(c(0,1),c(.9,.1))==0,is.na(auc_fixed(c(1,1),c(.1,.2))))
# Repeated patient observations stay together when bootstrap clusters are sampled.
source(project_path("code/01_gene_selection_ml/1.5_auc_uncertainty.R"))
bootstrap<-cluster_bootstrap_auc(rep(c(0,1),each=20),c(seq(0,.4,length.out=20),seq(.6,1,length.out=20)),rep(paste0("p",1:20),each=2),B=200L)
stopifnot(bootstrap$AUC==1,bootstrap$CI95_low==1,bootstrap$CI95_high==1)
cat("PASS: five AUC rows; 220 immune comparisons; 88 meta rows; AUC edge cases; cluster bootstrap diagnostic.\n")
