if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# ==== 0.2 运行单细胞流程 ====
# Execution scope: Entry point; prerequisites remain unchanged
# See CHAPTERS.md for order, dependencies and entry points.

# Reference-based reconstruction. Requires real input and checked parameters.
cmd<-commandArgs(trailingOnly=FALSE);script<-sub("^--file=","",cmd[grepl("^--file=",cmd)])
CODE_ROOT<-dirname(normalizePath(script,winslash="/",mustWork=TRUE))
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"));source(project_path("code/03_single_cell_knockout/3.1_single_cell_config.R"))
for(m in c("io","single_cell"))source(module_path(m))
cfg<-single_cell_config
if(!isTRUE(cfg$parameters_confirmed))stop("Confirm the study's actual input and QC/integration parameters in code/03_single_cell_knockout/3.1_single_cell_config.R; defaults only reflect the supplied references.")
object<-load_single_cell(cfg);qc<-qc_single_cell(object,cfg)
object<-reduce_single_cell(qc$object,cfg)
save_single_cell_outputs(object,qc$qc,cfg,project_path("outputs","single_cell"))
