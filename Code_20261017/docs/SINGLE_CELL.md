# Single-cell reconstruction from author-supplied references

The author confirmed on 2026-10-14 that the earlier target-expression modification scripts were not used in the study. Those scripts and workflows reading their `up3` objects are excluded. The four newly supplied scripts are treated as workflow references, not proof of exact historical parameters or figure provenance.

## Run order

1. Provide genuine 10X matrices and a filled sample manifest, OR a genes-by-cells integer count CSV with a separate cell-to-sample metadata table. The example files are schemas, not study data.
2. Check `code/03_single_cell_knockout/3.1_single_cell_config.R` against the actual study. Set `parameters_confirmed=TRUE` only after checking those values.
3. Run `Rscript run_single_cell.R`.
4. Supply the actual cluster-to-cell-type mapping. No automatic mapping is invented. Run again to save annotated output if needed.
5. Confirm the knockout population, gene-filter and network parameters; set `ko_parameters_confirmed=TRUE` and run `Rscript run_virtual_knockout.R`.

## Explicit corrections to the references

- Input paths and sample IDs are explicit. `GSE158631_count.csv` in the CSV example is not assumed to be the manuscript dataset.
- The original 10X script applied two QC filters sequentially. The reference-derived defaults retain their combined effect: nFeature >200 and <2500; mitochondrial fraction <5%; ribosomal fraction <=20%. Its nCount>=1000 condition was defined but never used; `min_counts=NULL` records that behavior. These are reference-derived values, not confirmed manuscript settings.
- PCA/normalization are performed once. When Harmony is requested, neighbors and UMAP explicitly use the Harmony reduction instead of silently reverting to PCA. This is a reconstruction correction, not a claim that previous figures used Harmony correctly.
- Sample IDs are preserved by manifest/metadata matching, not by guessing identity levels. The pipeline never invents biological sample labels.
- Each knockout uses its own target-specific saved result. The reference's FOS target / FOXO1 output / TP53 reload mismatch is removed; no unrelated result is loaded.
- The reference's gene-filter code used 1%, although its comment said 5%. This is now an explicit configuration value.
- Knockout QC is applied explicitly before choosing the network subsample size. The wrapper then calls `qc=FALSE` to avoid undocumented double filtering and incompatible QC argument names. This is a prospective implementation change; its configuration and included cells are saved.
- Both nominal p values and BH-adjusted p values are retained. Top 30 means the largest perturbation distances, with no artificial padding if fewer genes exist. The tested gene universe is saved for subsequent enrichment.
- The `FC` statistic is not interpreted as signed expression up/down regulation. It describes a distance-based perturbation measure; see the [official differential-regulation documentation](https://github.com/cailab-tamu/scTenifoldKnk/blob/master/man/dRegulation.Rd).

## Current validation boundary

QC, dimensional reduction and preservation of raw count values are tested with explicitly synthetic data. The real count matrix, final cell annotations and historical analysis settings are not included. scTenifoldKnk is not installed in the packaging environment, so its full engine has not been executed; its wrapper is syntax-checked against the documented API. No original UMAP, knockout result or Fig. 4 has been represented as reproduced.

GO/KEGG enrichment of knockout results still requires a finalized tested-gene universe, identifier mapping and actual enrichment output. No enrichment statistics are fabricated.
