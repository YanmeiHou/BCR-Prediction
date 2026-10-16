# ==== 0.1 项目路径设置 ====
# Execution scope: 基础函数
# See CHAPTERS.md for order, dependencies and entry points.

# Paths are relative to the release directory; no user-specific paths are used.
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  locale_result <- suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
  if (!nzchar(locale_result)) stop("A UTF-8 locale is required to read the supplied Unicode filenames.")
}
options(encoding="UTF-8") # Source chapter files consistently on Windows and Unix.
# ---- 0.1.1 project_path ----
project_path <- function(...) file.path(CODE_ROOT, ...)
# Restore bundled inputs on demand. Existing local files are never overwritten.
# The compact upload contains data/analysis_inputs.zip instead of 80 loose files.
input_archive <- project_path("data", "analysis_inputs.zip")
if(file.exists(input_archive)) {
  entries <- utils::unzip(input_archive,list=TRUE)$Name
  unsafe <- grepl("^/|^[A-Za-z]:|(^|/)\\.\\.(/|$)|\\\\",entries)
  if(any(unsafe)) stop("Unsafe path in bundled data archive")
  entries <- entries[!grepl("/$",entries)]
  missing <- entries[!file.exists(project_path("data",entries))]
  if(length(missing)) {
    utils::unzip(input_archive,files=missing,exdir=project_path("data"),overwrite=FALSE)
    if(!all(file.exists(project_path("data",entries))))stop("Input extraction was incomplete")
    message("Restored ",length(missing)," bundled input files.")
  }
}

# ---- 0.1.2 data_path ----
data_path <- function(...) project_path("data", "derived", ...)
# ---- 0.1.3 output_path ----
output_path <- function(...) project_path("outputs", "reproduction", ...)
# ---- 0.1.4 legacy_output ----
legacy_output <- function() {
  p <- project_path("outputs", "historical_replay")
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
  p
}

# ---- 0.1.5 Numbered module locations ----
module_path <- function(name) {
  routes <- c(
    io="code/00_setup/0.2_data_io.R",
    sample_audit="code/00_setup/0.3_sample_audit.R",
    model_audit="code/01_gene_selection_ml/1.4_saved_model_audit.R",
    auc_uncertainty="code/01_gene_selection_ml/1.5_auc_uncertainty.R",
    immune_statistics="code/02_clinical_enrichment_immune/2.7_immune_statistics.R",
    single_cell="code/03_single_cell_knockout/3.2_single_cell.R",
    virtual_knockout="code/03_single_cell_knockout/3.3_virtual_knockout.R"
  )
  if(length(name)!=1L || !name %in% names(routes)) stop("Unknown module: ",name)
  project_path(unname(routes[[name]]))
}
