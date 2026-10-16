# ==== 3.2 单细胞读取、QC、Harmony、聚类与绘图 ====
# Execution scope: 合成数据检查通过；需真实输入
# See CHAPTERS.md for order, dependencies and entry points.

# Refactored from the four user-supplied references.
# No target-specific expression values are increased, decreased, or injected.
# ---- 3.2.1 sc_counts ----
sc_counts <- function(object) {
  if(utils::packageVersion("SeuratObject")>="5.0.0") {
    layers<-SeuratObject::Layers(object[["RNA"]],search="^counts")
    if(length(layers)>1)object<-SeuratObject::JoinLayers(object,assay="RNA",layers=layers,new="counts")
    return(SeuratObject::GetAssayData(object,assay="RNA",layer="counts"))
  }
  SeuratObject::GetAssayData(object,assay="RNA",slot="counts")
}
# ---- 3.2.2 validate_counts ----
validate_counts <- function(counts) {
  if(is.null(rownames(counts))||is.null(colnames(counts))||anyDuplicated(rownames(counts))||anyDuplicated(colnames(counts)))stop("Unique gene and cell names are required")
  values<-if(inherits(counts,"sparseMatrix"))counts@x else as.vector(counts)
  if(any(!is.finite(values))||any(values<0)||any(abs(values-round(values))>1e-8))stop("Input must be non-negative raw integer counts, not log-normalized expression")
}
# ---- 3.2.3 load_single_cell ----
load_single_cell <- function(cfg) {
  require_packages(c("Seurat","SeuratObject","Matrix"))
  if(cfg$input_format=="10x") {
    manifest<-read.csv(project_path(cfg$manifest),stringsAsFactors=FALSE)
    assert_columns(manifest,c("sample_id","path"));assert_unique(manifest,"sample_id")
    if(!nrow(manifest))stop("Fill the sample manifest with the actual input paths")
    objects<-list()
    for(i in seq_len(nrow(manifest))) {
      sample_id<-as.character(manifest$sample_id[i])
      if(!grepl("^[A-Za-z0-9.-]+$",sample_id))stop("Use stable alphanumeric sample IDs")
      p<-project_path(manifest$path[i]);if(!dir.exists(p))stop("Missing 10X directory: ",p)
      counts<-Seurat::Read10X(p)
      if(is.list(counts)) {
        if(!"Gene Expression" %in% names(counts))stop("10X multimodal input has no Gene Expression matrix")
        counts<-counts[["Gene Expression"]]
      }
      validate_counts(counts)
      object<-Seurat::CreateSeuratObject(counts=counts,project=sample_id,min.cells=cfg$create_min_cells,min.features=cfg$create_min_features)
      object$sample_id<-sample_id
      objects[[sample_id]]<-Seurat::RenameCells(object,add.cell.id=sample_id)
    }
    object<-if(length(objects)==1)objects[[1]]else merge(objects[[1]],y=objects[-1],merge.data=FALSE)
  } else if(cfg$input_format=="csv") {
    counts<-as.matrix(read.csv(project_path(cfg$csv_counts),row.names=1,check.names=FALSE))
    storage.mode(counts)<-"numeric";validate_counts(counts)
    metadata<-read.csv(project_path(cfg$csv_metadata),stringsAsFactors=FALSE)
    assert_columns(metadata,c("cell_id","sample_id"));assert_unique(metadata,"cell_id")
    idx<-match(colnames(counts),metadata$cell_id)
    if(anyNA(idx)||anyNA(metadata$sample_id[idx])||any(!nzchar(metadata$sample_id[idx])))stop("Every count-matrix cell needs a sample_id")
    object<-Seurat::CreateSeuratObject(counts=Matrix::Matrix(counts,sparse=TRUE),min.cells=cfg$create_min_cells,min.features=cfg$create_min_features)
    object$sample_id<-metadata$sample_id[match(colnames(object),metadata$cell_id)]
  } else stop("input_format must be 10x or csv")
  if(utils::packageVersion("SeuratObject")>="5.0.0")object<-SeuratObject::JoinLayers(object,assay="RNA")
  object
}
# ---- 3.2.4 qc_single_cell ----
qc_single_cell <- function(object,cfg) {
  if(!any(grepl(cfg$mt_pattern,rownames(object))))stop("No mitochondrial symbols matched; confirm gene ID mapping before QC")
  if(!any(grepl(cfg$ribo_pattern,rownames(object))))stop("No ribosomal symbols matched; confirm gene ID mapping before QC")
  object[["percent.mt"]]<-Seurat::PercentageFeatureSet(object,pattern=cfg$mt_pattern)
  object[["percent.rb"]]<-Seurat::PercentageFeatureSet(object,pattern=cfg$ribo_pattern)
  md<-object[[]]
  keep<-md$nFeature_RNA>cfg$min_features_exclusive & md$nFeature_RNA<cfg$max_features_exclusive & md$percent.mt<cfg$max_mt_percent_exclusive & md$percent.rb<=cfg$max_ribo_percent_inclusive
  if(!is.null(cfg$min_counts))keep<-keep & md$nCount_RNA>=cfg$min_counts
  if(!any(keep))stop("No cells pass the configured QC thresholds")
  qc<-data.frame(cell_id=rownames(md),md[c("sample_id","nCount_RNA","nFeature_RNA","percent.mt","percent.rb")],pass_QC=keep,row.names=NULL)
  list(object=object[,keep],qc=qc)
}
# ---- 3.2.5 reduce_single_cell ----
reduce_single_cell <- function(object,cfg) {
  require_packages(c("Seurat","harmony"));set.seed(cfg$seed)
  object<-Seurat::NormalizeData(object,normalization.method="LogNormalize",scale.factor=cfg$scale_factor,verbose=FALSE)
  object<-Seurat::FindVariableFeatures(object,selection.method="vst",nfeatures=cfg$variable_features,verbose=FALSE)
  object<-Seurat::ScaleData(object,features=rownames(object),verbose=FALSE)
  npcs<-max(20L,max(cfg$dimensions))
  if(ncol(object)<=npcs+1L)stop("Too few cells for the configured PCA dimensions")
  object<-Seurat::RunPCA(object,features=Seurat::VariableFeatures(object),npcs=npcs,seed.use=cfg$seed,verbose=FALSE)
  reduction<-"pca"
  if(cfg$use_harmony) {
    if(length(unique(object$sample_id))<2L)stop("Harmony requested but fewer than two batches are present")
    object<-harmony::RunHarmony(object,group.by.vars="sample_id",reduction.use="pca",dims.use=cfg$dimensions,reduction.save="harmony",verbose=FALSE)
    reduction<-"harmony"
  }
  object<-Seurat::FindNeighbors(object,reduction=reduction,dims=cfg$dimensions,verbose=FALSE)
  object<-Seurat::FindClusters(object,resolution=cfg$resolution,random.seed=cfg$seed,verbose=FALSE)
  object<-Seurat::RunUMAP(object,reduction=reduction,dims=cfg$dimensions,n.neighbors=cfg$umap_neighbors,seed.use=cfg$seed,verbose=FALSE)
  object
}
# ---- 3.2.6 save_single_cell_outputs ----
save_single_cell_outputs <- function(object,qc,cfg,destination) {
  require_packages("ggplot2");dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  write.csv(qc,file.path(destination,"cell_QC.csv"),row.names=FALSE)
  markers<-Seurat::FindAllMarkers(object,only.pos=TRUE,min.pct=cfg$marker_min_pct,logfc.threshold=cfg$marker_logfc)
  write.csv(markers,file.path(destination,"cluster_markers.csv"),row.names=FALSE)
  # Cell types require an explicit author-supplied mapping, not automatic relabeling.
  annotation_path<-project_path(cfg$annotation_file)
  if(file.exists(annotation_path)) {
    annotation<-read.csv(annotation_path,stringsAsFactors=FALSE)
    assert_columns(annotation,c("cluster","cell_type"));assert_unique(annotation,"cluster")
    idx<-match(as.character(object$seurat_clusters),as.character(annotation$cluster))
    if(anyNA(idx)||anyNA(annotation$cell_type[idx])||any(!nzchar(annotation$cell_type[idx])))stop("Annotation must cover every cluster")
    object$cell_type<-annotation$cell_type[idx]
  }
  ggplot2::ggsave(file.path(destination,"UMAP_clusters.pdf"),Seurat::DimPlot(object,reduction="umap",label=TRUE),width=8,height=6)
  genes<-intersect(cfg$target_genes,rownames(object))
  if(length(genes))ggplot2::ggsave(file.path(destination,"target_gene_expression.pdf"),Seurat::FeaturePlot(object,features=genes,ncol=2),width=10,height=8)
  saveRDS(object,file.path(destination,"single_cell_processed.rds"))
  dput(cfg,file=file.path(destination,"used_config.R"))
  capture.output(sessionInfo(),file=file.path(destination,"sessionInfo.txt"))
  invisible(object)
}
