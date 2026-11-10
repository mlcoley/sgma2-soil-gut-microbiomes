# SGMA2: code for analysis of child stool, cattle dung and household soil resistome

Code for the analyses in the manuscript on 194 shotgun metagenomes from 39 households in
Kisumu county, western Kenya. Code has two parts:

- `pipeline/`: SLURM scripts that take host-filtered reads to generate gene-by-sample feature tables
  (host filtering, taxonomy, assembly, gene catalog, abundance, annotation, resistance genes).
- `analysis/`: R scripts that implement the analyses, using as inputs the feature tables and household metadata.

Reads with human sequences removed are deposited in SRA [PRJNA1551843]. Feature tables and
de-identified household, child and soil metadata are archived at OSF [10.17605/OSF.IO/3G2Z5](https://doi.org/10.17605/OSF.IO/3G2Z5);
`inputs/README.md` lists each input file and its columns.

## Layout

```
config/paths.sh        scratch root and database directory for the cluster scripts
pipeline/              cluster stages 00-13 (see pipeline/README.md)
analysis/              R scripts 01-27, run_all.R, utilities (see analysis/README.md)
inputs/README.md       every input file, its columns and where it comes from
```

## Running the pipeline

1. Edit `config/paths.sh` (`SGMA2_BASE`, `SGMA2_DB_DIR`).
2. In every `#SBATCH` header replace `YOUR_ACCOUNT` and `YOUR_EMAIL` and adjust the partition
   and QOS (`condo` on the cluster used here).
3. Submit each script from the repository root so `logs/` and `config/paths.sh` resolve:
   `mkdir -p logs && sbatch pipeline/<stage>/<script>`. Array scripts take `--array=1-N`
   where N is the number of samples (or assembly manifest rows, or protein chunks).

Stages and their order are in `pipeline/README.md`.

## Running the analysis

```
export SGMA2_DATA=/path/to/data          # feature tables, QC tables, metadata (see inputs/README.md)
export SGMA2_HOUSEHOLD=/path/to/household  # household metadata
export SGMA2_OUT=/path/to/output
Rscript analysis/run_all.R
```

Run from the repository root. Defaults are `./data`, `./data/household` and `./output`.
Outputs are written to `$SGMA2_OUT/tables`, `figures`, `logs`, `tables_docx` and
`supplementary_data`. The background comparison in scripts 18-23 requires HPC;
`analysis/run_all.R` documents the steps. Download the OSF project files and point `SGMA2_DATA` and `SGMA2_HOUSEHOLD` at its folders.

## Software

- HPC: fastp 1.3.3, minimap2 2.30, samtools 1.6, seqkit 2.13.0, FastQC 0.12.1, MultiQC 1.35,
  bowtie2 2.5.5, Woltka 0.1.7 (conda env `metagenomic_env`); MetaPhlAn 4.1.2
  (`biobakery`); MEGAHIT 1.2.9, QUAST 5.3.0 (`assembly`); Prodigal, MMseqs2 18.8cc5c, CoverM 0.7.0,
  eggNOG-mapper 2.1.13 (`genecat`); AMRFinderPlus 4.2.7 (`amr`); MicrobeCensus (`microbecensus2`). Environment files: `pipeline/00_setup/`.
  Full list: Supplementary Table S7.
- Analysis: R 4.5.2 with vegan 2.7.5, ggplot2 4.0.3, patchwork 1.3.2, dplyr 1.2.1, tidyr 1.3.2,
  scales 1.4.0, lme4 2.0.1, data.table 1.18.4, xgboost 3.2.1.1, FEAST 0.1.0 (GitHub, commit
  2f8f3df8051e0e08341f597a9f4693bfb76b3bf6), officer 0.7.4, flextable 0.9.11, readr 2.2.0,
  purrr 1.2.2, stringr 1.6.0. `run_all.R` writes `sessionInfo.txt` to the log directory.
  openxlsx is loaded only by an unused workbook writer in `analysis/util_manuscript_tables.R`.

## Reproducibility notes

- Every script that draws random numbers sets its own seed.
- `run_all.R` was run end to end on the deposited inputs. Every manuscript table and figure matches
  the version used in the manuscript, with three exceptions:
  - The deposited household files give child age in whole months, not dates of birth, so the
    child-age rows of Data S3 (script 16) differ in the third decimal place. No conclusion changes.
  - The pair classifier (24) forks workers for its permutation null, which is not
    seed-reproducible. Its AUCs and its calibration table (`T40c`) match exactly; its permutation
    p-values differ slightly between runs (Table S15).
  - The jittered points in Figure 2 fall in different positions; the data, boxes and statistics
    are identical.
- SRA reads are the raw reads filtered against GRCh38 and T2T-Y only (`SRA_HUMAN_ONLY=true` in `pipeline/01_host-filter-qc/01_host_filter_array.sh`). All analyses use the fastp-trimmed, cattle- and human-filtered reads.
- One library was re-run through host filtering with `pairfilt.awk`; the other 255 ran before it
  was added. The filter passes clean read pairs through unchanged, so the two versions give the
  same output on libraries without zero-length reads.
