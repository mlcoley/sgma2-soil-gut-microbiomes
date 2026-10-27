#!/usr/bin/env Rscript
# Target profile for the matched non-ARG background comparison (step 1 of 4).
# Gut and soil harbour different phyla, so any 95%-identity gene cluster tends to be
# compartment-specific; the containment null holds the feature set fixed and cannot
# say whether ARGs behave differently from other genes. The comparison needs non-ARG
# catalog genes matched to the clinical ARG variants. This script characterises the
# ARG side only (the distribution background genes must be matched to) and prepares
# the two files the cluster summary job needs.
# Matching axes: per-compartment prevalence (the triple, not the pooled scalar),
# abundance, gene length, and family granularity (clusters per ARG family; gene length
# is only available on the cluster).
# Steps: this script -> pipeline/12_background-genes (summary job) -> 19_background_match.R
# -> pipeline/12_background-genes/extract_background_abundance.sbatch ->
# 20_background_containment.R -> 21_fig4.R.
# Outputs: T23a (ARG target profile), T23b (family granularity), T23c (cohort map) and
# T23d (AMRFinderPlus gene ids), the last two uploaded to $SGMA2_BASE/assembly on the
# cluster.
# Usage: Rscript analysis/18_background_target_profile.R   (from the repository root)

source("analysis/00_load_features.R")

OUT_PREFIX <- "T23"
GATE_GENES <- TRUE   # apply the gene-level negative-control exclusion list

cat("\n=== ARG target profile ===\n\n")

## ---- cohort --------------------------------------------------------------

key  <- load_key()
# load_key() returns 195: it still carries P181 (dung, 25,052 reads), which the
# analysis drops by intersecting with the feature tables. Do that explicitly here,
# since the cohort list this script emits drives the cluster job.
ab_ids <- unique(read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")),
                            stringsAsFactors = FALSE)$sample_id)
dropped <- setdiff(key$sample_id, ab_ids)
if (length(dropped)) {
  cat(sprintf("Dropping %d key sample(s) absent from the feature tables: %s\n",
              length(dropped), paste(dropped, collapse = ", ")))
  key <- key[key$sample_id %in% ab_ids, ]
}
stopifnot(nrow(key) == 194)
comp <- setNames(as.character(key$compartment), key$sample_id)
n_by_comp <- table(factor(comp, levels = COMP_LEVELS))
cat(sprintf("Cohort: %d samples  (gut %d, dung %d, soil %d)\n",
            length(comp), n_by_comp["gut"], n_by_comp["dung"], n_by_comp["soil"]))

## ---- ARG gene -> symbol map ---------------------------------------------

hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")),
                   stringsAsFactors = FALSE)
# One gene can carry multiple hit rows (same symbol, different method/accession).
map <- unique(hits[, c("gene_id", "arg_symbol")])
dup <- map$gene_id[duplicated(map$gene_id)]
if (length(dup)) {
  cat(sprintf("%d catalog genes carry >1 arg_symbol; keeping the first per gene.\n",
              length(unique(dup))))
  map <- map[!duplicated(map$gene_id), ]
}
cat(sprintf("ARG-bearing catalog genes (all symbols): %d across %d symbols\n",
            nrow(map), length(unique(map$arg_symbol))))

## ---- clinical restriction -----------------------------------------------

clin <- clinical_args()
map$is_clinical <- map$arg_symbol %in% clin
cat(sprintf("  clinical subset (AMR|core): %d genes across %d symbols\n",
            sum(map$is_clinical), length(unique(map$arg_symbol[map$is_clinical]))))

## ---- abundance, gated ----------------------------------------------------

ab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")),
                 stringsAsFactors = FALSE)
cat(sprintf("Abundance rows (sparse, count>0): %d\n", nrow(ab)))

if (GATE_GENES) {
  bl <- read.delim(file.path(FT_DIR, "control_arg_genes_blocklist.tsv"),
                   stringsAsFactors = FALSE)
  n0 <- nrow(ab)
  ab <- ab[!(ab$gene_id %in% bl$gene_id), ]
  cat(sprintf("Gene-level control blocklist: %d genes listed, %d abundance rows dropped\n",
              nrow(bl), n0 - nrow(ab)))
}

ab <- ab[ab$sample_id %in% names(comp), ]
ab$compartment <- comp[ab$sample_id]
cat(sprintf("Rows on the analysed cohort: %d\n\n", nrow(ab)))

## ---- per-gene profile ----------------------------------------------------
# Presence is count >= 1; the breadth rule (covered_fraction >= 0.5) is carried as a
# parallel column.

ab$present_cf50 <- ab$covered_fraction >= 0.5

prev_tab <- function(rows, tag) {
  t <- table(factor(rows$gene_id, levels = unique(ab$gene_id)),
             factor(rows$compartment, levels = COMP_LEVELS))
  m <- as.data.frame.matrix(t)
  names(m) <- paste0(tag, "_", names(m))
  m
}

genes <- sort(unique(ab$gene_id))
P1  <- prev_tab(ab, "prev")
P50 <- prev_tab(ab[ab$present_cf50, ], "prevcf50")

agg <- function(f) tapply(ab$count, ab$gene_id, f)
prof <- data.frame(
  gene_id           = rownames(P1),
  P1[, paste0("prev_", COMP_LEVELS)],
  P50[rownames(P1), paste0("prevcf50_", COMP_LEVELS)],
  total_count       = as.numeric(agg(sum))[match(rownames(P1), names(agg(sum)))],
  mean_nonzero      = as.numeric(agg(mean))[match(rownames(P1), names(agg(mean)))],
  stringsAsFactors  = FALSE
)
prof$prev_total <- rowSums(prof[, paste0("prev_", COMP_LEVELS)])
prof$arg_symbol  <- map$arg_symbol[match(prof$gene_id, map$gene_id)]
prof$is_clinical <- prof$arg_symbol %in% clin
prof <- prof[!is.na(prof$arg_symbol), ]

cat(sprintf("Detected ARG genes (>=1 read, gated, on cohort): %d  (clinical %d)\n",
            nrow(prof), sum(prof$is_clinical)))

cl <- prof[prof$is_clinical, ]

## ---- family granularity ---------------------------------------------------
# Clusters per family: the catalog-wide count (arg_annotation$n_genes) and the detected
# count (genes seen in >= 1 cohort sample), which is what the containment test sees.

ann <- load_arg_annotation()
det <- as.data.frame(table(cl$arg_symbol), stringsAsFactors = FALSE)
names(det) <- c("arg_symbol", "n_genes_detected")
det$n_genes_catalog <- ann$n_genes[match(det$arg_symbol, ann$arg_symbol)]
det <- det[order(-det$n_genes_detected), ]

cat(sprintf("\nClinical symbols with >=1 detected gene: %d (median %d detected clusters per family)\n",
            nrow(det), median(det$n_genes_detected)))

## ---- write ---------------------------------------------------------------

save_tab(prof, sprintf("%sa_arg_target_profile.tsv", OUT_PREFIX))
save_tab(det,  sprintf("%sb_arg_family_granularity.tsv", OUT_PREFIX))

# Two small files the cluster job needs.
cohort <- data.frame(sample_id = sort(names(comp)),
                     compartment = comp[sort(names(comp))],
                     stringsAsFactors = FALSE)
write.table(cohort, file.path(TAB_DIR, "T23c_cohort_compartments.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nWrote cohort map (%d samples) -> T23c_cohort_compartments.tsv\n",
            nrow(cohort)))

# Every AMRFinderPlus-called catalog gene at any element type and scope, so the
# background pool cannot contain a gene recognised as a resistance, stress or
# virulence determinant.
arg_ids <- sort(unique(hits$gene_id))
writeLines(arg_ids, file.path(TAB_DIR, "T23d_amrfinder_gene_ids.txt"))
cat(sprintf("Wrote AMRFinderPlus gene-id exclusion list (%d) -> T23d_amrfinder_gene_ids.txt\n",
            length(arg_ids)))
cat("\nDone.\n")
