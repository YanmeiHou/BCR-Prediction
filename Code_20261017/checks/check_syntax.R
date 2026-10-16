if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
# Parse scripts only; no analysis runs or automatic package installation.
files <- c(list.files("code",pattern="[.]R$",recursive=TRUE,full.names=TRUE),list.files(".",pattern="[.]R$",full.names=TRUE),list.files("checks",pattern="[.]R$",full.names=TRUE))
for(f in files) tryCatch(parse(f,encoding="UTF-8"),error=function(e)stop(f,": ",conditionMessage(e)))
cat("PASS: ",length(files)," R scripts parsed.\n",sep="")
