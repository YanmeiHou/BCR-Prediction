if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# Check outputs added for the approved revision; run from repository root.
CODE_ROOT <- normalizePath(".",winslash="/",mustWork=TRUE)
source(file.path(CODE_ROOT,"code/00_setup/0.1_project_paths.R"))
p <- read.csv(output_path("Figure3B_raw_p.csv"))
expected <- c(0.00548630101307578,0.0235979259108711,0.114094637767485,0.0513117101154839,0.513812257611678,0.107608273309937)
stopifnot(nrow(p)==6L,max(abs(p$p_raw-expected))<1e-9)
cat("PASS: six Figure 3B raw p values.\n")
