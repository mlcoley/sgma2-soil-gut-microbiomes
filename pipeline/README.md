# Cluster pipeline

SLURM scripts, run from the repository root after editing `config/paths.sh` and the
`#SBATCH` headers (see the top-level README). Stages run in numeric order; within a stage,
scripts are numbered in run order. Array scripts skip samples with a done marker, so a
failed array can be resubmitted. `$SGMA2_BASE` is the scratch root.

| Stage | Scripts | What it does | Main output |
|---|---|---|---|
| 00_setup | 01-04 | one-time conda environments (`biobakery`, `assembly`, `genecat`, `amr`) and the AMRFinderPlus database; `environment_metagenomic.yml` (fastp 1.3.3, FastQC 0.12.1, MultiQC 1.35, minimap2 2.30, samtools 1.6, seqtk, seqkit 2.13.0, bowtie2 2.5.5, pigz, Woltka 0.1.7) and `environment_microbecensus2.yml` (Python 2.7) | environments |
| 01_host-filter-qc | `01_host_filter_array.sh`, `pairfilt.awk` | fastp (adapters, minimum length 100), then minimap2 against cattle, GRCh38 and T2T-Y, keeping unmapped pairs; reads `fastq/sample_list.tsv`. With `SRA_HUMAN_ONLY=true` there is no trimming and no cow-filter step, which gives the data prepared for SRA upload | `host_filter/filtered_fastq/<sample>_{1,2}.filt.fq.gz`; SRA version in `host_filter_human_only/filtered_fastq/` |
| 02_post-host-qc-accounting | 01-03, `build_accounting_table.py` | sample sheet, seqkit and FastQC per library, MultiQC, read-accounting table | `posthost_qc/sgma2_step1_read_accounting.tsv` |
| 03_taxonomy-woltka | 01-02, `merge_woltka_tables.py` | bowtie2 to Web of Life 2, Woltka classification, per-genome coverage, cohort tables | `taxonomy/merged/cohort_*.tsv` |
| 04_taxonomy-metaphlan | 01-02, `metaphlan_classified_summary.py` | MetaPhlAn 4.1.2 profiles and merged tables | `biobakery/merged_metaphlan/` |
| 06_assembly-megahit | `01_assembly_array.sh` | per-sample MEGAHIT assembly (contigs >= 1000 bp) and QUAST; needs `assembly/assembly_manifest.tsv` | `assembly/assemblies/<sample>/` |
| 07_genecatalog-prodigal-mmseqs | 01-02 | Prodigal genes with sample-prefixed ids; MMseqs2 `easy-linclust` catalog at 95% identity and 90% coverage; gene length table | `assembly/genecat/catalog_*` |
| 08_abundance-coverm | 01-02, `merge_coverm.py` | CoverM read count and covered fraction for every catalog gene per sample; merge to one sparse long table (`MANIFEST` and `DONE` can be overridden to map the control libraries) | `genecat/gene_abundance_long.tsv.gz` |
| 09_function-eggnog | 01-04, `build_ko_table.py` | catalog proteins, protein chunks, eggNOG-mapper 2.1.13, merged annotations, KO x sample table (`MIN_COVERED_FRACTION=0.5` for the breadth-filtered table) | `genecat/eggnog/ko_abundance_wide.tsv.gz` |
| 10_amr-amrfinder | 01-03, `merge_amr.py`, `build_arg_table.py` | AMRFinderPlus 4.2.7 on the protein chunks, merged hits, ARG x sample table | `genecat/amr/arg_abundance_wide.tsv.gz`, `arg_hits_long.tsv.gz` |
| 11_feature-tables | six scripts | breadth-filtered KO and ARG tables; ARG gene abundance slice; control (negative blank) feature profiles and candidate exclusion lists; ARG re-clustering at 90, 85 and 80% identity | `*_cf50.tsv.gz`, `arg_gene_abundance.tsv.gz`, `controls/`, `arg_reclust/` |
| 12_background-genes | two scripts | per-gene summary of the catalog (prevalence, abundance, length, eggNOG OG and KO) and abundance rows for the selected background genes; see the background steps in `analysis/run_all.R` | `background/` |
| 13_microbecensus | array, merge | MicrobeCensus average genome size and genome equivalents (environment `microbecensus2`, Python 2.7) | `microbecensus_ge.tsv` |

Stage 05 is not part of this release. HUMAnN, the assembly pilot and the residual non-bovine
host screen were not part of the reported analysis and are not included.

## Notes

- Control libraries (extraction blanks, PCR blanks, positive aliquots, mocks) are mapped to the
  catalog with a separate marker directory (`DONE`):
  `MANIFEST=$SGMA2_BASE/assembly/control_manifest.tsv DONE=$SGMA2_BASE/assembly/done_coverm_controls sbatch --array=1-N%8 pipeline/08_abundance-coverm/01_coverm_array.sh`.
  Control outputs stay in `assembly/coverm/`. `gene_abundance_long.tsv.gz` holds
  `assembly_manifest.tsv` samples only.
- `11_feature-tables/build_control_features.sbatch` and `build_control_gene_features.sbatch`
  write `*.CANDIDATE.tsv` exclusion lists to `assembly/genecat/controls/`. Each lists features
  detected at >= 5 reads in any extraction or PCR blank (`DEFAULT_MIN_READS`); the
  `*_impact.tsv` files give counts at 1 to 50 reads. The lists applied in the analysis hold only
  the ARG rows with `max_count_in_negative_control` >= 5 (KO features are not filtered); save them in
  `feature_tables/` as `control_features_blocklist.tsv` or `control_arg_genes_blocklist.tsv`.
- `03_taxonomy-woltka/calculate-coverage.py` computes per-genome coverage from Woltka `--outcov`
  output. It is copied from [qiita-spots/qp-woltka](https://github.com/qiita-spots/qp-woltka/blob/222d68d84e0000fce63c5b4ead567e04ae4dbaf2/qp_woltka/util.py) (BSD 3-clause).
- `arg_gene_abundance.tsv.gz` is not breadth-filtered. It carries `covered_fraction` for
  filtering in analysis scripts.
- The T2T-Y step in `01_host_filter_array.sh` writes final pair files without `pairfilt.awk`;
  cattle and GRCh38 steps use it.
- CoverM 0.7.0 in `genecat` bundles minimap2 2.31; `metagenomic_env` has minimap2 2.30.
- `genecat/gene_abundance_long.tsv.gz` has no header; the scripts that read it detect a header
  rather than assume one.
- Shared databases under `$SGMA2_DB_DIR`: `human/{grch38-db-may2024,t2ty-db}.mmi`,
  `cow/{taurus-indicus,taurus,indicus}-db.mmi`, `wol2/`, `metaphlan/` (index
  mpa_vJun23_CHOCOPhlAnSGB_202307), `eggnog/` (eggNOG 5.0; eggnog.db and eggnog_proteins.dmnd).
  eggNOG-mapper needs `--dmnd_iterate no` with the format-v2 diamond index.
- Do not submit an array twice while it is running: the MetaPhlAn array has no in-progress lock.
