# Historical implementation disclosure

These scripts preserve selected supplied implementation logic while standardizing UTF-8 encoding, filenames, input paths and output locations. They are not part of `run_reproduction.R` and have not been rerun on the missing original objects. They must not be cited as a corrected, unbiased modelling pipeline.

## Why execution is opt-in

Historical decisions remain unresolved: predefined candidate genes differ from the saved RF/LASSO intersection; model ranking uses external AUCs; the best saved forest has three inputs; LODO uses only the classifier half of a combination; and tissue/patient inclusion has not been finalized. We have not silently replaced those choices with invented history.

The pooled ROC plotting branch uses pre-reordering `test_idx` positions on reordered data. The main tested reproduction route instead matches actual sample IDs and groups by cohort, avoiding that plotting branch. This does not repair historical figures automatically.

The historical nominal/final AUC tables, full-cohort scaling, selector fallbacks and swallowed fitting errors are retained for disclosure, not endorsed as final validation. The historical nomogram sources fit four-gene logistic models; one combines all cohorts. They must not be represented as a verified seven-variable clinical nomogram.

## Input layout

- `data/external/discovery/`: the root-level original TCGA ExpressionSet and PCaDB gene annotation object used by feature-discovery code.
- `data/external/modelling/`: the four ExpressionSets used by the modelling code. Keep these separate from discovery objects unless checksums/content establish equivalence.
- `data/external/CIBERSORT/`: appropriately obtained `CIBERSORT.R` and `LM22.txt`, required only for historical deconvolution.
- `data/external/h.all.v2026.1.Hs.symbols.gmt`: the matching Hallmark resource, required only for GSEA.
- Outputs/intermediate objects are under `outputs/historical_replay/`; downstream scripts expect their corresponding prepared/model objects there.

No bundled third-party signature or GMT license is newly granted. Obtain these resources under their applicable terms. Existing CIBERSORT fractions are sufficient for the tested reproduction route.

For an intentional historical replay after reviewing the missing inputs and limitations:

```r
Sys.setenv(PCA_ENABLE_HISTORICAL="1", PCA_CODE_ROOT=normalizePath("."))
source("code/01_gene_selection_ml/1.3_model_benchmark.R", encoding="UTF-8")
```

This is historical replay, not the prescribed future reanalysis. Do not execute all scripts indiscriminately. Full model binaries and historical package versions are not included. `docs/source_manifest.csv` maps chapter sources; `docs/original_archive_inventory.csv` preserves the original archive inventory and its earlier inclusion decisions. Current scope and updates are described in the README and execution guide.
