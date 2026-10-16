# Numbered R chapters

Use the root entry points listed in docs/RUNNING.md; chapter order is organizational.

| Chapter | R file | Description |
|---|---|---|
| 0.1 | `code/00_setup/0.1_project_paths.R` | 0.1 项目路径设置 |
| 0.2 | `code/00_setup/0.2_data_io.R` | 0.2 数据读取与输入检查 |
| 0.3 | `code/00_setup/0.3_sample_audit.R` | 0.3 样本数量与患者重复核对 |
| 1.1 | `code/01_gene_selection_ml/1.1_feature_discovery.R` | 1.1 差异表达、LASSO与随机森林筛选 |
| 1.2 | `code/01_gene_selection_ml/1.2_ML_functions.R` | 1.2 机器学习公共函数 |
| 1.3 | `code/01_gene_selection_ml/1.3_model_benchmark.R` | 1.3 模型训练、比较与验证 |
| 1.4 | `code/01_gene_selection_ml/1.4_saved_model_audit.R` | 1.4 现有预测AUC、模型记录与筛选交集核对 |
| 1.5 | `code/01_gene_selection_ml/1.5_auc_uncertainty.R` | 1.5 AUC患者簇Bootstrap置信区间 |
| 1.6 | `code/01_gene_selection_ml/1.6_original_OOB_figures.R` | 1.6 Original saved-model OOB and AUC figures |
| 1.7 | `code/01_gene_selection_ml/1.7_TCGA495_fixed_RF_OOB.R` | 1.7 TCGA495 fixed RF OOB |
| 1.8 | `code/01_gene_selection_ml/1.8_nested_model_tuning.R` | 1.8 nested model tuning |
| 1.9 | `code/01_gene_selection_ml/1.9_nested_metrics_figures.R` | 1.9 nested metrics figures |
| 2.1 | `code/02_clinical_enrichment_immune/2.1_survival.R` | 2.1 生存分析 |
| 2.2 | `code/02_clinical_enrichment_immune/2.2_cox_forest.R` | 2.2 Cox回归与森林图 |
| 2.3 | `code/02_clinical_enrichment_immune/2.3_single_gene_ROC.R` | 2.3 单基因ROC |
| 2.4 | `code/02_clinical_enrichment_immune/2.4_single_gene_GSEA.R` | 2.4 单基因GSEA |
| 2.5 | `code/02_clinical_enrichment_immune/2.5_immune_deconvolution.R` | 2.5 免疫细胞比例推断 |
| 2.6 | `code/02_clinical_enrichment_immune/2.6_immune_correlation.R` | 2.6 基因与免疫细胞相关性 |
| 2.7 | `code/02_clinical_enrichment_immune/2.7_immune_statistics.R` | 2.7 免疫组间比较、Meta分析与BH校正 |
| 2.8 | `code/02_clinical_enrichment_immune/2.8_nomogram_DCA.R` | 2.8 Nomogram/DCA audit and reconstruction from original expression/clinical data. |
| 2.9 | `code/02_clinical_enrichment_immune/2.9_historical_nomogram_TCGA.R` | 2.9 TCGA列线图、校准与DCA |
| 2.10 | `code/02_clinical_enrichment_immune/2.10_historical_nomogram_combined.R` | 2.10 合并队列列线图与DCA |
| 2.11 | `code/02_clinical_enrichment_immune/2.11_Figure2B_ROC.R` | Figure 2B: cohort-specific ROC curves, preserving the original five-model comparison. |
| 2.12 | `code/02_clinical_enrichment_immune/2.12_original_pooled_DCA_intervals.R` | out <- project_path('outputs','original_pooled_dca') |
| 2.13 | `code/02_clinical_enrichment_immune/2.13_historical_drug_sensitivity.R` | 2.13 Historical oncoPredict implementation; manuscript linkage unconfirmed. |
| 3.1 | `code/03_single_cell_knockout/3.1_single_cell_config.R` | 3.1 单细胞及虚拟敲除参数 |
| 3.2 | `code/03_single_cell_knockout/3.2_single_cell.R` | 3.2 单细胞读取、QC、Harmony、聚类与绘图 |
| 3.3 | `code/03_single_cell_knockout/3.3_virtual_knockout.R` | 3.3 四个基因的虚拟敲除 |
