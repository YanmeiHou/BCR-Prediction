if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# 1.6 Optional original-model OOB plots.
cmd <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
CODE_ROOT <- if(length(script)) dirname(normalizePath(script,winslash="/",mustWork=TRUE)) else normalizePath(".",winslash="/",mustWork=TRUE)
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"))
source(module_path("io"))
source(project_path("code/01_gene_selection_ml/1.6_original_OOB_figures.R"))
