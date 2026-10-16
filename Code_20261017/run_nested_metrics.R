if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
options(encoding="UTF-8")
cmd <- commandArgs(FALSE)
script <- sub('^--file=', '', cmd[grepl('^--file=', cmd)])
CODE_ROOT <- if(length(script)) dirname(normalizePath(script,winslash='/')) else normalizePath('.',winslash='/')
source(file.path(CODE_ROOT,'code/00_setup/0.1_project_paths.R'),encoding='UTF-8')
out <- project_path('outputs','nested_tuning')
for(s in c('tables','figures')) dir.create(file.path(out,s),recursive=TRUE,showWarnings=FALSE)
for(n in c('outer_predictions.csv','final_model_predictions.csv')) {
 p<-file.path(out,'tables',n)
 if(!file.exists(p)) file.copy(project_path('data','reference','nested_tuning',n),p)
}
source(project_path('code/01_gene_selection_ml/1.9_nested_metrics_figures.R'),encoding='UTF-8')
