# ==== 3.3 四个基因的虚拟敲除 ====
# Execution scope: 需真实对象及scTenifoldKnk
# See CHAPTERS.md for order, dependencies and entry points.

# Separate knockout runs for each target; no unrelated result is loaded afterward.
# FC is a perturbation-distance statistic, NOT signed expression fold change.
# ---- 3.3.1 virtual_knockout ----
virtual_knockout <- function(object,cfg,destination) {
  require_packages(c("scTenifoldKnk","Seurat","SeuratObject","Matrix"))
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  if(!is.null(cfg$ko_cell_type)) {
    if(!"cell_type" %in% names(object[[]]))stop("Cell-type annotations are required for a restricted knockout analysis")
    cells<-colnames(object)[object$cell_type %in% cfg$ko_cell_type]
    if(!length(cells))stop("No cells match ko_cell_type")
    object<-object[,cells]
  }
  counts<-sc_counts(object);validate_counts(counts)
  # Explicit QC before selecting a feasible network subsample size.
  mt<-grepl(cfg$mt_pattern,rownames(counts));library_size<-Matrix::colSums(counts)
  ratio<-Matrix::colSums(counts[mt,,drop=FALSE])/pmax(library_size,1)
  keep_cells<-library_size>=cfg$ko_qc_min_library_size & ratio<=cfg$ko_qc_max_mt_ratio
  if(sum(keep_cells)<50)stop("Fewer than 50 cells pass the configured knockout QC")
  counts<-counts[,keep_cells,drop=FALSE];object<-object[,colnames(counts)]
  if(cfg$ko_use_hvg)object<-Seurat::FindVariableFeatures(object,selection.method="vst",nfeatures=cfg$ko_n_hvg,verbose=FALSE)
  expressed<-rownames(counts)[Matrix::rowMeans(counts>0)>=cfg$ko_min_fraction]
  engine<-scTenifoldKnk::scTenifoldKnk
  formal_names<-names(formals(engine))
  if(!all(c("countMatrix","gKO","qc","nc_nNet","nc_nCells") %in% formal_names))stop("Unsupported scTenifoldKnk API; review installed version")
  audit<-list()
  for(target in cfg$target_genes) {
    if(!target %in% rownames(counts))stop("Target is absent from raw counts: ",target)
    if(!any(counts[target,]>0))stop("Target has zero observed expression: ",target)
    keep<-unique(c(expressed,target))
    if(cfg$ko_use_hvg)keep<-intersect(keep,unique(c(Seurat::VariableFeatures(object),target)))
    matrix<-counts[keep,,drop=FALSE]
    # Counts are unchanged; only rows/cells are selected. No artificial expression.
    stopifnot(identical(as.matrix(matrix[target,,drop=FALSE]),as.matrix(counts[target,,drop=FALSE])))
    set.seed(cfg$ko_seed)
    args<-list(countMatrix=matrix,gKO=target,qc=FALSE,nc_nNet=cfg$ko_n_networks,nc_nCells=min(cfg$ko_cells_per_network,ncol(matrix)))
    if("nCores" %in% formal_names)args$nCores<-cfg$ko_cores
    if("seed" %in% formal_names)args$seed<-cfg$ko_seed
    result<-do.call(engine,args)
    saveRDS(result,file.path(destination,paste0(target,"_knockout.rds")))
    tab<-as.data.frame(result$diffRegulation)
    assert_columns(tab,c("gene","distance","FC","p.value"))
    if(!"p.adj" %in% names(tab))tab$p.adj<-p.adjust(tab$p.value,method="BH")
    tab<-tab[tab$gene!=target,,drop=FALSE]
    tab<-tab[order(tab$distance,decreasing=TRUE),,drop=FALSE]
    write.csv(tab,file.path(destination,paste0(target,"_all_regulation.csv")),row.names=FALSE)
    write.csv(tab[!is.na(tab$p.adj)&tab$p.adj<.05,,drop=FALSE],file.path(destination,paste0(target,"_FDR05.csv")),row.names=FALSE)
    write.csv(head(tab,30),file.path(destination,paste0(target,"_top30_by_distance.csv")),row.names=FALSE)
    writeLines(rownames(matrix),file.path(destination,paste0(target,"_tested_gene_universe.txt")))
    audit[[target]]<-data.frame(target=target,cells=ncol(matrix),genes=nrow(matrix),expressed_fraction=mean(matrix[target,]>0),seed=cfg$ko_seed,n_networks=cfg$ko_n_networks,n_cells_per_network=args$nc_nCells)
  }
  write.csv(do.call(rbind,audit),file.path(destination,"knockout_run_summary.csv"),row.names=FALSE)
  write.csv(data.frame(cell_id=colnames(counts)),file.path(destination,"included_cells.csv"),row.names=FALSE)
  dput(cfg,file.path(destination,"used_config.R"));capture.output(sessionInfo(),file=file.path(destination,"sessionInfo.txt"))
}
