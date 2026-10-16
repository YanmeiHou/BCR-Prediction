# ==== 3.1 单细胞及虚拟敲除参数 ====
# Execution scope: 需核对实际参数
# See CHAPTERS.md for order, dependencies and entry points.

# Reference-derived settings, not confirmed historical settings for the manuscript.
# Set parameters_confirmed=TRUE only after checking the actual analysis records.
single_cell_config <- list(
  parameters_confirmed=FALSE,
  input_format="10x", # "10x" or "csv"
  manifest="data/external/single_cell/sample_manifest.csv",
  csv_counts="data/external/single_cell/counts.csv", # genes x cells, symbols as rows
  csv_metadata="data/external/single_cell/cell_metadata.csv", # cell_id,sample_id
  create_min_cells=3L, create_min_features=200L,
  # Effective intersection of BOTH QC filters in the reference script:
  min_features_exclusive=200L, max_features_exclusive=2500L,
  max_mt_percent_exclusive=5, max_ribo_percent_inclusive=20,
  # nCount>=1000 was defined but NOT included in the reference filter.
  min_counts=NULL, mt_pattern="^MT-", ribo_pattern="^RP",
  variable_features=2000L, scale_factor=10000, dimensions=1:10,
  use_harmony=TRUE, resolution=1.2, umap_neighbors=30L, seed=777L,
  marker_min_pct=.25, marker_logfc=.25,
  target_genes=c("TOP2A","DNAH8","CHRDL1","PGM5"),
  annotation_file="data/external/single_cell/cluster_annotation.csv",
  ko_parameters_confirmed=FALSE,
  ko_cell_type=NULL, # Explicitly choose a cell population or document all cells.
  ko_min_fraction=.01, # The reference CODE used 1%; its 5% comment was inconsistent.
  ko_use_hvg=TRUE, ko_n_hvg=5000L,
  ko_qc_max_mt_ratio=.1, ko_qc_min_library_size=500L,
  ko_n_networks=10L, ko_cells_per_network=2000L,
  ko_cores=2L, ko_seed=777L
)
