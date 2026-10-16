if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# 2. Immune statistics from supplied fractions; no CIBERSORT rerun.
cmd <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
CODE_ROOT <- if(length(script)) dirname(normalizePath(script,winslash="/",mustWork=TRUE)) else normalizePath(".",winslash="/",mustWork=TRUE)
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"))
source(module_path("io"))

source(module_path("immune_statistics"))
immune_statistics(read_data("immune_fractions.csv"),read_data("samples.csv"))
x <- read.csv(output_path("immune_meta.csv"))
focus <- c("T cells follicular helper","T cells CD8","Plasma cells","Macrophages M2","Dendritic cells resting","T cells regulatory (Tregs)")
x <- x[x$population=="original_1025" & x$effect_type=="Cohen_d",]
x <- x[match(focus,x$cell_type),c("cell_type","k_cohorts","pooled_effect","p_raw")]
stopifnot(!anyNA(x))
write_result(x,"Figure3B_raw_p.csv")
capture.output(sessionInfo(),file=output_path("sessionInfo_immune.txt"))
