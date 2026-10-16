if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# ==== 8.2 单细胞合成数据流程检查 ====
# Execution scope: 检查脚本
# See CHAPTERS.md for order, dependencies and entry points.

# Synthetic implementation test only. No synthetic data are saved in study outputs.
CODE_ROOT<-normalizePath(".",winslash="/",mustWork=TRUE)
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"));source(project_path("code/03_single_cell_knockout/3.1_single_cell_config.R"))
for(m in c("io","single_cell"))source(module_path(m))
require_packages(c("Seurat","harmony","Matrix"))
set.seed(1)
counts<-matrix(rpois(600*100,lambda=1),600,100)
rownames(counts)<-c(paste0("MT-",1:10),paste0("RPL",1:10),paste0("gene",1:580));colnames(counts)<-paste0("cell",1:100)
cfg<-single_cell_config;cfg$variable_features<-200L;cfg$resolution<-.2
o<-Seurat::CreateSeuratObject(counts=Matrix::Matrix(counts,sparse=TRUE));o$sample_id<-rep(c("batch1","batch2"),each=50)
q<-qc_single_cell(o,cfg)
before<-sc_counts(q$object)
r<-reduce_single_cell(q$object,cfg)
after<-sc_counts(r)
stopifnot(identical(dim(before),dim(after)),all(before==after),"harmony" %in% names(r@reductions),"umap" %in% names(r@reductions))
# Invalid normalized values must be rejected as counts.
bad<-counts;bad[1,1]<-.5;stopifnot(inherits(try(validate_counts(bad),silent=TRUE),"try-error"))
cat("PASS: synthetic QC, Harmony-based neighbors/UMAP, unchanged counts, invalid-count rejection. scTenifoldKnk is not tested here.\n")
