#!/usr/bin/env Rscript
# Source tables for three supplementary tables.
#   T44a  what the contamination exclusion lists removed, at each resolution
#   T44b  clinical ARG abundance ratios between compartments with and without the
#         exclusion list (needs T_arg_abundance_filter_sensitivity.tsv from 02_descriptive.R)
#   T44c  reference-based taxonomic classification rates by compartment and in the mock
#         communities (from feature_tables/metaphlan_classified_summary.tsv)
# Usage: Rscript analysis/26_supplement_source_tables.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/00_load_features.R")

set.seed(1)

log_lines <- c(sprintf("Supplement source tables - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

ann <- load_arg_annotation()
key <- load_key()

## ---------------------------------------------------------------------------
## A. What contamination filtering removed at each resolution
## ---------------------------------------------------------------------------
GAT  <- load_gating()                                     # family level
GGAT <- load_gating("control_arg_genes_exclusion_list.tsv", "gene-level")

M   <- load_matrix(file.path(FT_DIR, LAYERS$arg$file))
M   <- M[rownames(M) %in% key$sample_id, , drop = FALSE]
det <- colnames(M)[colSums(M) > 0]                        # detected in the analysed cohort
cl  <- unique(clinical_args(ann))

gab  <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")), stringsAsFactors = FALSE)
gab  <- gab[gab$sample_id %in% key$sample_id, ]
hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")), stringsAsFactors = FALSE)
g2s  <- setNames(hits$arg_symbol, hits$gene_id)

genes_all <- unique(gab$gene_id)
genes_cl  <- unique(gab$gene_id[gab$gene_id %in% names(g2s)[g2s %in% cl]])

t44a <- data.frame(
  resolution = c("ARG family", "ARG family", "ARG variant", "ARG variant"),
  feature_set = c("all detected", "clinical subset", "all detected", "clinical subset"),
  exclusion_list = c(rep("family level", 2), rep("gene level", 2)),
  n_before = c(length(det), length(intersect(det, cl)),
               length(genes_all), length(genes_cl)),
  n_excluded = c(length(intersect(det, GAT$features)),
                 length(intersect(intersect(det, cl), GAT$features)),
                 length(intersect(genes_all, GGAT$features)),
                 length(intersect(genes_cl, GGAT$features))),
  stringsAsFactors = FALSE)
t44a$n_after <- t44a$n_before - t44a$n_excluded
t44a$pct_removed <- round(100 * t44a$n_excluded / t44a$n_before, 1)

stopifnot(t44a$n_before[1] == 303, t44a$n_after[1] == 288,
          t44a$n_before[2] == 181, t44a$n_after[2] == 170,
          t44a$n_after[4] == 1681)
save_tab(t44a, "T44a_contamination_filter_features.tsv")

## ---------------------------------------------------------------------------
## B. What contamination filtering did to the between-compartment contrast
## ---------------------------------------------------------------------------
fs <- read.delim(file.path(TAB_DIR, "T_arg_abundance_filter_sensitivity.tsv"),
                 stringsAsFactors = FALSE)
stopifnot(setequal(fs$compartment, COMP_LEVELS))
g <- function(cp, col) fs[[col]][match(cp, fs$compartment)]

t44b <- data.frame(
  contrast = c("child stool / soil", "child stool / cattle dung"),
  ratio_filtered = c(g("gut", "clinical_filtered")   / g("soil", "clinical_filtered"),
                     g("gut", "clinical_filtered")   / g("dung", "clinical_filtered")),
  ratio_unfiltered = c(g("gut", "clinical_unfiltered") / g("soil", "clinical_unfiltered"),
                       g("gut", "clinical_unfiltered") / g("dung", "clinical_unfiltered")),
  stringsAsFactors = FALSE)
t44b$direction <- ifelse(t44b$ratio_filtered > t44b$ratio_unfiltered,
                         "filtering widens the contrast", "filtering narrows the contrast")
stopifnot(length(unique(t44b$direction)) == 2)      #  they differ
save_tab(t44b, "T44b_contamination_filter_contrast.tsv")

## ---------------------------------------------------------------------------
## C. Reference-based taxonomic classification rates
## ---------------------------------------------------------------------------
# mock libraries are the comparator: classification is near-complete on a known community
mp <- read.delim(file.path(FT_DIR, "metaphlan_classified_summary.tsv"),
                 stringsAsFactors = FALSE)
mock <- mp[mp$compartment == "control_mock", ]

# the analysed cohort is the key intersected with the feature tables
analysed <- intersect(key$sample_id, rownames(M))
stopifnot(length(analysed) == 194)
bio <- mp[mp$sample_id %in% analysed & mp$compartment %in% COMP_LEVELS, ]

missing_mp <- setdiff(analysed, mp$sample_id)
stopifnot(nrow(mock) == 4, nrow(bio) == 193, identical(missing_mp, "P151"))
logmsg("taxonomic rates cover %d of the %d analysed samples; %s has no MetaPhlAn profile",
       nrow(bio), length(analysed), paste(missing_mp, collapse = ", "))

t44c <- bio %>%
  mutate(compartment = factor(compartment, levels = COMP_LEVELS)) %>%
  group_by(compartment) %>%
  summarise(n = n(),
            pct_classified_median = median(pct_classified),
            pct_classified_min    = min(pct_classified),
            pct_classified_max    = max(pct_classified),
            n_species_median      = median(n_species),
            .groups = "drop") %>%
  mutate(sample_type = "study sample") %>%
  as.data.frame()

t44c <- rbind(t44c, data.frame(
  compartment = "mock community",
  n = nrow(mock),
  pct_classified_median = median(mock$pct_classified),
  pct_classified_min    = min(mock$pct_classified),
  pct_classified_max    = max(mock$pct_classified),
  n_species_median      = median(mock$n_species),
  sample_type = "ZymoBIOMICS mock",
  stringsAsFactors = FALSE))
t44c <- t44c[, c("sample_type", "compartment", "n", "pct_classified_median",
                 "pct_classified_min", "pct_classified_max", "n_species_median")]


stopifnot(abs(round(t44c$pct_classified_median[1:3]) - c(83, 10, 3)) < 1e-9,
          floor(t44c$pct_classified_min[4]) == 84,
          floor(t44c$pct_classified_max[4]) == 93)
save_tab(t44c, "T44c_taxonomic_classification.tsv")

logmsg("")
logmsg("== contamination filtering ==")
for (i in seq_len(nrow(t44a)))
  logmsg("  %-12s %-16s %5d -> %5d  (%d removed, %.1f%%)", t44a$resolution[i],
         t44a$feature_set[i], t44a$n_before[i], t44a$n_after[i],
         t44a$n_excluded[i], t44a$pct_removed[i])
for (i in seq_len(nrow(t44b)))
  logmsg("  %-26s %.2fx filtered vs %.2fx unfiltered; %s", t44b$contrast[i],
         t44b$ratio_filtered[i], t44b$ratio_unfiltered[i], t44b$direction[i])
logmsg("")
logmsg("== taxonomic classification ==")
for (i in seq_len(nrow(t44c)))
  logmsg("  %-16s n=%3d  median %5.1f%% classified (range %.1f-%.1f), median %.0f species",
         t44c$compartment[i], t44c$n[i], t44c$pct_classified_median[i],
         t44c$pct_classified_min[i], t44c$pct_classified_max[i], t44c$n_species_median[i])

writeLines(unlist(log_lines), file.path(LOG_DIR, "26_supplement_source_tables.log"))
cat("\nSupplement source tables complete.\n")
