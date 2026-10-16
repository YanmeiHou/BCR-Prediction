# Execution guide

Run from the package root. Inputs are read-only except optional missing-input extraction branches; outputs go to `outputs/`.

| Command after `Rscript --vanilla` | Purpose | Inputs / limits |
|---|---|---|
| `run_reproduction.R` | Saved-score AUC, original feature audit, immune statistics | Bundled CSV; metafor |
| `run_gene_ml.R` | Saved-score/model record audit only | Bundled CSV |
| `run_immune.R` | Immune comparison/meta-analysis; raw p values exported separately | Bundled CSV; metafor |
| `run_oob.R` | Historical 547-record RF-stage OOB diagnostics | Bundled saved predictions; not 495-patient analysis |
| `run_tcga495.R` | Retrain fixed historical 3-predictor RF on 495 unique primary tumors | Bundled 495-patient input; seed 777, 1000 trees |
| `run_nested_metrics.R` | Recalculate nested-run metrics and figures | Bundled saved outer/final predictions; no new tuning |
| `run_nested_tuning.R --retrain` | Full nested tuning: 5 outer, 4 inner folds | Bundled 495-patient input; expensive; fixed historical features |
| `run_clinical.R` | Seven-variable nomogram, bootstrap, separate four-gene DCA | Bundled clinical data; 500 bootstrap refits per configured model |
| `run_figure2b.R` | 20 cohort/model ROC comparisons from saved scores | Bundled scores; historical training/cross-cohort results |
| `run_original_dca.R` | Original four-cohort pooled-fit DCA intervals | Bundled prepared data; APPARENT performance |
| `run_historical_drug.R` | Historical oncoPredict source, not confirmed final manuscript workflow | Opt-in + external cohort and GDSC2 data |
| `run_single_cell.R` | Real-data single-cell reconstruction | External data, annotation and confirmed config required |
| `run_virtual_knockout.R` | Virtual knockout | Processed object + confirmed cell population/parameters required |
| `run_auc_ci.R --confirmed-final-input` | Patient-cluster bootstrap AUC interval | Separately reviewed prediction table required |

Historical chapters 1.1, 1.3, 2.1–2.6, 2.9–2.10 and 2.13 require review before `PCA_ENABLE_HISTORICAL=1`. Do not set that flag globally for routine reproduction.
