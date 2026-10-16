# ==== 0.2 数据读取与输入检查 ====
# Execution scope: 基础函数
# See CHAPTERS.md for order, dependencies and entry points.

# Shared input validation and CSV output. No packages are installed automatically.
# ---- 0.2.1 require_packages ----
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Install required packages first: ", paste(missing, collapse=", "))
}
# ---- 0.2.2 read_data ----
read_data <- function(filename) {
  p <- data_path(filename)
  if (!file.exists(p)) stop("Missing input: ", p)
  data<-read.csv(p, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8")
  for(column in intersect(c("primary_tumor_sensitivity","restored_zero_column","saved_in_Publication"),names(data))) {
    values<-toupper(as.character(data[[column]]))
    if(anyNA(values)||!all(values %in% c("TRUE","FALSE")))stop("Invalid logical flag in ",column)
    data[[column]]<-values=="TRUE"
  }
  data
}
# ---- 0.2.3 write_result ----
write_result <- function(x, filename) {
  dir.create(output_path(), recursive=TRUE, showWarnings=FALSE)
  write.csv(x, output_path(filename), row.names=FALSE, na="", fileEncoding="UTF-8")
}
# ---- 0.2.4 assert_unique ----
assert_unique <- function(data, columns) {
  if (anyDuplicated(data[columns])) stop("Duplicate key: ", paste(columns, collapse=" + "))
}
# ---- 0.2.5 assert_columns ----
assert_columns <- function(data, columns) {
  if (!all(columns %in% names(data))) stop("Missing columns: ", paste(setdiff(columns,names(data)),collapse=", "))
}
