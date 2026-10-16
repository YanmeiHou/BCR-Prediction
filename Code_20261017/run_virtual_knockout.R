if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# ==== 0.3 运行虚拟敲除 ====
# Execution scope: Entry point; prerequisites remain unchanged
# See CHAPTERS.md for order, dependencies and entry points.

# Requires the checked single-cell output and an installed scTenifoldKnk package.
cmd<-commandArgs(trailingOnly=FALSE);script<-sub("^--file=","",cmd[grepl("^--file=",cmd)])
CODE_ROOT<-dirname(normalizePath(script,winslash="/",mustWork=TRUE))
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"));source(project_path("code/03_single_cell_knockout/3.1_single_cell_config.R"))
for(m in c("io","single_cell","virtual_knockout"))source(module_path(m))
cfg<-single_cell_config
if(!isTRUE(cfg$parameters_confirmed)||!isTRUE(cfg$ko_parameters_confirmed))stop("Confirm actual single-cell and knockout parameters first")
p<-project_path("outputs","single_cell","single_cell_processed.rds")
if(!file.exists(p))stop("Run the checked single-cell workflow first")
virtual_knockout(readRDS(p),cfg,project_path("outputs","virtual_knockout"))
