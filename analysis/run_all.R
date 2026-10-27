#!/usr/bin/env Rscript
# Runs the analysis scripts in dependency order from the repository root.
#   Rscript analysis/run_all.R            # scripts 01-25
#   Rscript analysis/run_all.R 19 25      # scripts 19 through 25
# 00_load_features.R and the util_*.R files are sourced by the scripts that need them.
# The background comparison (18-21) needs a cluster round trip, so run it in stages:
#   1. 18_background_target_profile.R writes T23c and T23d; copy them to
#      $SGMA2_BASE/assembly and run pipeline/12_background-genes/background_gene_summary.sbatch
#      (and copy background_gene_summary.tsv.gz and arg_gene_lengths.tsv to data/feature_tables).
#   2. 19_background_match.R writes background_gene_ids.txt; copy it to
#      $SGMA2_BASE/assembly, run pipeline/12_background-genes/extract_background_abundance.sbatch
#      and copy background_abundance*.tsv.gz to data/feature_tables.
#   3. 20_background_containment.R, 21_fig4.R.
# 22_pair_classifier.R takes about 100 minutes at R8_NCORE=6 and is memory-bound; run it
# with nothing else in R alongside it.
# Scripts 15-17 read household metadata (see inputs/README.md); 05 must run before 15.
set.seed(1)

scripts <- c(
  "01_build_single_copy_markers.R",
  "02_descriptive.R",
  "03_diversity.R",
  "04_nestedness.R",
  "05_exposure.R",
  "06_normalization.R",
  "07_rare_feature_sharing.R",
  "08_qc_batch_audit.R",
  "09_breadth_sensitivity.R",
  "10_clinical_subset.R",
  "11_unweighted_distance.R",
  "12_soil_analyte_reliability.R",
  "13_allele_source_sensitivity.R",
  "14_allele_threshold_sweep.R",
  "15_exposure_multivariate.R",
  "16_child_covariates.R",
  "17_diet.R",
  "18_background_target_profile.R",
  "19_background_match.R",
  "20_background_containment.R",
  "21_fig4.R",
  "22_pair_classifier.R",
  "23_feast_source_tracking.R",
  "24_supplement_source_tables.R",
  "25_manuscript_tables.R"
)

args  <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
first <- if (length(args) >= 1 && !is.na(args[1])) args[1] else 1L
last  <- if (length(args) >= 2 && !is.na(args[2])) args[2] else length(scripts)
stopifnot(first >= 1, last <= length(scripts), first <= last)

for (s in scripts[first:last]) {
  cat("\n================  running", s, " ================\n")
  source(file.path("analysis", s), local = new.env())
}
dir.create(file.path(Sys.getenv("SGMA2_OUT", "output"), "logs"), recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(Sys.getenv("SGMA2_OUT", "output"), "logs", "sessionInfo.txt"))
cat("\nrun_all complete.\n")
