# Validation — 17 October 2026

## Executed in this release

- All R files parsed successfully (final count in `validation/release_summary.json`).
- The dry-lab reproduction entry point completed. Tests passed for five saved AUC
  rows, 220 immune comparisons, 88 meta-analysis rows, AUC direction/ties and a
  clustered-bootstrap diagnostic; the six Figure 3B raw p values matched references.
- The current clinical reconstruction ran with 500 bootstrap refits per configured
  model. Seven-variable TCGA apparent AUC 0.7403324 and corrected AUC 0.7158090,
  DKFZ AUC 0.7803819; clinical metrics and exact DCA intervals match references.
- One fixed RF was retrained on all 495 unique TCGA primary-tumor patients.
  Training AUC 0.9940057; RF-stage OOB AUC 0.6426640. All per-patient predictions,
  metrics and the threshold-0.50 confusion matrix match saved references.
- Historical 547-record saved-model OOB figures regenerated (OOB AUC 0.6747610).
  This analysis is explicitly separate from the 495-patient reconstruction.
- Figure 2B reproduced all 20 cohort/model AUCs from saved scores.
- Original pooled-fit DCA threshold tables and intervals reproduced exactly.
- Nested-run metrics/figures and paired fixed-prediction bootstrap comparisons
  were regenerated from stored predictions. The full nested search was NOT rerun.
- Fourteen output tables matched saved references within numerical tolerance;
  see `validation/output_checks.csv`.

## Limits

The historical discovery/ML pipeline, raw GSEA/deconvolution, drug-sensitivity
pipeline, real single-cell workflow and scTenifoldKnk were not rerun here.
scTenifoldKnk is absent in this host environment; install a documented compatible
version before executing that optional workflow. Actual single-cell input/settings
remain unconfirmed. A syntax check is not a scientific validation of these modules.
Nested tuning uses fixed historical features; it is not nested feature discovery.
Reference agreement does not resolve feature identity or dataset-selection issues.

The host emits unsupported inherited-locale warnings on R startup. Entrypoints now
establish Windows UTF-8 before sourcing chapters; the listed runs completed.
`validation/package_versions.csv` and `sessionInfo.txt` describe this validation
environment, not an invented historical software environment.

Generated plots/models are retained locally outside the upload package. The ZIP
contains code, small inputs, reference tables, documentation and validation records.
