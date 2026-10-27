# Analysis scripts

Run from the repository root; see the top-level README for the environment variables and
`inputs/README.md` for the input files. `run_all.R` runs scripts 01-25 in order
(`Rscript analysis/run_all.R 19 25` runs a range). Every script sources
`00_load_features.R` (paths, key, loaders, display labels). Output file names keep the table
(`T..`) and figure (`F..`) codes used internally; the mapping to manuscript items is below.

| Script | Content | Manuscript items |
|---|---|---|
| `prepare_household_tables.R` | household tables from the long questionnaire exports | input to 05, 15-17 |
| 01 | single-copy marker KOs (KEGG BRITE ko03011) | Table S4 |
| 02 | cohort inventory, clinical ARG abundance and tests, filter sensitivity, alpha diversity | Fig 2, Fig S2, Tables S2, S8, S11, S12 |
| 03 | PERMANOVA, dispersion, within-soil household distances | Fig S1, Tables S9, S10 |
| 04 | equal-effort containment and feature counts | Fig 3, Table S13 |
| 05 | exposure and soil-property associations | Figs S4, S5, Data S2, S6 |
| 06 | MicrobeCensus against single-copy normalization | Table S5 |
| 07 | rare-feature household sharing and influence checks | Fig 5, Fig S3, Data S1 |
| 08 | QC and batch audit; plate-restricted permutation | Supp. Methods S2.2 |
| 09 | read count vs covered-fraction presence rules | Supp. Methods S6.2 |
| 10 | clinical subset definitions | Table S1 |
| 11 | presence-only and abundance-weighted distances | Fig 1, Tables S9, S10, Supp. Methods S6.3 |
| 12 | ICRAF prediction reliability | Table S6 |
| 13 | catalog-source sensitivity of ARG-variant containment | Supp. Methods S6.1 |
| 14 | containment at 95, 90, 85, 80% catalog identity | Supp. Methods S6.4 |
| 15 | exposure associations with composition | Data S5 |
| 16 | child covariates | Data S3 |
| 17 | diet | Data S4 |
| 18-21 | matched non-ARG background: target profile, matched sets, containment, figure | Fig 4, Supp. Methods S5.2 |
| 22 | pair classifier and SHAP | Figs S7, S8, Table S15 |
| 23 | FEAST source tracking | Fig S6, Table S14 |
| 24 | contamination-filter, ratio and taxonomy source tables | Tables S2, S3 |
| 25 | renders Tables S1-S15 (.docx) and Data S1-S6 (.csv) | all supplementary items |

Table S7 (software versions) is read from `data/metadata/software_versions.tsv`.

## Run notes

- 19 takes about 12 minutes and 20 about 14; 22 takes about 100 minutes at `R8_NCORE=6` and is
  memory-bound (run it with nothing else in R); 23 takes about 3 minutes at `R5_NCORE=6`.
- Scripts that call a figure builder (`util_fig_*.R`) draw from tables already on disk and
  compute nothing, so a label fix never changes a result.
- Gene-level analyses use the gene-level exclusion list and family-level analyses the symbol-level
  list; never apply one to the other.
