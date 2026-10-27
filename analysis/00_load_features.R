#!/usr/bin/env Rscript
# Shared paths, constants and loaders. Sourced by every analysis script.
# Run all scripts from the repository root.
#
# Environment variables (all optional):
#   SGMA2_DATA       input directory (default: ./data)
#   SGMA2_HOUSEHOLD  household metadata directory (default: $SGMA2_DATA/household)
#   SGMA2_OUT        output directory (default: ./output)
# Inputs and layout are described in inputs/README.md.

suppressPackageStartupMessages({
  library(vegan); library(ggplot2)
})

## ---- paths ---------------------------------------------------------------
SGMA2_ROOT <- getwd()
DATA_DIR   <- Sys.getenv("SGMA2_DATA", file.path(SGMA2_ROOT, "data"))
OUT_DIR    <- Sys.getenv("SGMA2_OUT",  file.path(SGMA2_ROOT, "output"))
FT_DIR     <- file.path(DATA_DIR, "feature_tables")
QC_DIR     <- file.path(DATA_DIR, "posthost_qc")
META_DIR   <- Sys.getenv("SGMA2_HOUSEHOLD", file.path(DATA_DIR, "household"))
KEY_FILE   <- file.path(DATA_DIR, "metadata", "sample_household_key.tsv")
FIG_DIR    <- file.path(OUT_DIR, "figures")
TAB_DIR    <- file.path(OUT_DIR, "tables")
LOG_DIR    <- file.path(OUT_DIR, "logs")
for (.sgma2_dir in c(FIG_DIR, TAB_DIR, LOG_DIR))
  dir.create(.sgma2_dir, recursive = TRUE, showWarnings = FALSE)
rm(.sgma2_dir)

# Feature layers. `scope` is where the layer is usable.
LAYERS <- list(
  tax_phylum  = list(file = "cohort_phylum.tsv",         type = "woltka", scope = "gut"),
  tax_genus   = list(file = "cohort_genus.tsv",          type = "woltka", scope = "gut"),
  tax_species = list(file = "cohort_species.tsv",        type = "woltka", scope = "gut"),
  ko_all      = list(file = "ko_abundance_wide.tsv.gz",  type = "woltka", scope = "all"),
  arg         = list(file = "arg_abundance_wide.tsv.gz", type = "woltka", scope = "all")
)

## ---- constants -----------------------------------------------------------
COMP_LEVELS <- c("gut", "dung", "soil")
COMP_COLS   <- c(gut = "#B5462F", dung = "#C08A2A", soil = "#3F6B57")
# Presence is count > MIN_COUNT. Breadth-filtered (covered fraction >= 0.5) tables
# are the parallel *_cf50 files and are reported alongside, not instead.
MIN_COUNT   <- 0

## ---- sample key ----------------------------------------------------------
# Biological samples only. soil_primary_only keeps soil microsites S01-S03; the
# excluded ids are attached as attr(key, "soil_excluded").
load_key <- function(biological_only = TRUE, soil_primary_only = TRUE) {
  stopifnot(file.exists(KEY_FILE))
  key <- read.delim(KEY_FILE, stringsAsFactors = FALSE)
  if (biological_only) key <- key[key$sample_type == "true_sample", ]
  key$soil_site <- suppressWarnings(as.integer(sub("^.*-S0*([0-9]+)$", "\\1", key$field_name)))
  key$soil_site[key$compartment != "soil"] <- NA
  excluded <- character(0)
  if (soil_primary_only) {
    drop <- key$compartment == "soil" & !is.na(key$soil_site) & key$soil_site > 3
    excluded <- key$sample_id[drop]
    key <- key[!drop, ]
  }
  key$compartment <- factor(key$compartment, levels = COMP_LEVELS)
  attr(key, "soil_excluded") <- excluded
  key
}

## ---- feature-table loader ------------------------------------------------
# Returns a samples x features matrix. type "metaphlan" reads the headerless merged
# MetaPhlAn tables, taking sample ids from metaphlan_classified_summary.tsv.
load_matrix <- function(path, type = "woltka") {
  stopifnot(file.exists(path))
  if (type == "metaphlan") {
    df <- read.delim(path, header = FALSE, check.names = FALSE,
                     comment.char = "#", quote = "")
    rownames(df) <- df[[1]]; df[[1]] <- NULL
    sm <- file.path(dirname(path), "metaphlan_classified_summary.tsv")
    if (!file.exists(sm))
      stop("metaphlan table has no header; need metaphlan_classified_summary.tsv alongside it")
    summ <- read.delim(sm, stringsAsFactors = FALSE)
    if (nrow(summ) != ncol(df))
      stop(sprintf("column count (%d) != classified_summary rows (%d)", ncol(df), nrow(summ)))
    if (!all(colSums(df > 0) == summ$n_species))
      stop("column order does not match classified_summary (n_species mismatch)")
    colnames(df) <- summ$sample_id
    m <- t(as.matrix(df)); storage.mode(m) <- "double"; m[is.na(m)] <- 0
    return(m)
  }
  # generic (Woltka counts / KO / ARG / HUMAnN CPM); read.delim reads .gz directly
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  raw   <- readLines(con, n = 50)
  hdr_i <- which(grepl("\t", raw))[1]
  df <- read.delim(path, header = TRUE, check.names = FALSE,
                   skip = hdr_i - 1, comment.char = "", quote = "")
  rownames(df) <- df[[1]]; df[[1]] <- NULL
  colnames(df) <- sub("^([Pp][0-9]+).*$", "\\1", colnames(df))   # normalize to bare P###
  num <- vapply(df, is.numeric, logical(1))
  m <- t(as.matrix(df[, num, drop = FALSE]))                     # samples x features
  storage.mode(m) <- "double"; m[is.na(m)] <- 0
  m
}

# Loads a named layer aligned to the key, dropping all-zero features. Returns
# list(M, key, name, scope). restrict_scope keeps only the layer's own compartment.
get_layer <- function(name, key = load_key(), restrict_scope = FALSE) {
  spec <- LAYERS[[name]]
  if (is.null(spec)) stop("unknown layer: ", name)
  M <- load_matrix(file.path(FT_DIR, spec$file), spec$type)
  common <- intersect(rownames(M), key$sample_id)
  if (length(common) < 5)
    stop(sprintf("only %d samples matched key for layer %s — check naming", length(common), name))
  M   <- M[common, , drop = FALSE]
  kk  <- key[match(common, key$sample_id), ]
  if (restrict_scope && spec$scope != "all") {
    keep <- kk$compartment == spec$scope
    M <- M[keep, , drop = FALSE]; kk <- kk[keep, ]
  }
  M <- M[, colSums(M) > 0, drop = FALSE]
  list(M = M, key = kk, name = name, scope = spec$scope)
}

## ---- ARG annotation ------------------------------------------------------
load_arg_annotation <- function() {
  f <- file.path(FT_DIR, "arg_annotation.tsv")
  stopifnot(file.exists(f))
  read.delim(f, stringsAsFactors = FALSE)
}

# Clinical subset: acquired AMR determinants at core scope.
clinical_args <- function(ann = load_arg_annotation()) {
  ann$arg_symbol[ann$element_type == "AMR" & ann$scope == "core"]
}

# Resistance categories: clinical AMR, other AMR, metal, biocide/other stress.
METAL_RE <- "MERCURY|COPPER|SILVER|ARSENIC|TELLURIUM|NICKEL|ZINC|CADMIUM|COBALT|GOLD|IRON|ALUMINUM|LEAD"
arg_category <- function(ann = load_arg_annotation()) {
  dc <- toupper(ifelse(is.na(ann$drug_class), "", ann$drug_class))
  cat <- ifelse(ann$element_type == "VIRULENCE", "Virulence",
         ifelse(ann$element_type == "AMR" & ann$scope == "core", "Clinical AMR",
         ifelse(ann$element_type == "AMR", "Other AMR (--plus)",
         ifelse(grepl(METAL_RE, dc), "Metal",
                "Biocide/other stress"))))
  setNames(cat, ann$arg_symbol)
}

## ---- negative-control filtering ------------------------------------------
# The symbol-level exclusion list holds ARG features detected at >= 5 reads in any
# negative control. KO is not filtered. Gene-level analyses use
# control_arg_genes_blocklist.tsv instead; never apply one list to the other level.
# The list's key column is `feature` (symbol level) or `gene_id` (gene level).
load_gating <- function(file = "control_features_blocklist.tsv", label = "symbol-level") {
  f <- file.path(FT_DIR, file)
  if (!file.exists(f))
    return(list(applied = FALSE, features = character(0),
                note = sprintf("NOT FILTERED (%s) — %s not present", label, file)))
  bl <- read.delim(f, stringsAsFactors = FALSE)
  ids <- if (!is.null(bl$feature)) bl$feature else bl[[1]]
  list(applied = TRUE, features = ids,
       note = sprintf("contamination-filtered against negative controls (%s, %d features)",
                      label, nrow(bl)))
}

# Drop blocked columns from a samples x features matrix.
drop_gated <- function(M, gating) {
  if (!gating$applied) return(M)
  M[, !colnames(M) %in% gating$features, drop = FALSE]
}

# Apply gating to a get_layer() result in place, re-dropping all-zero features
# and warning if any sample is emptied (which would break distance metrics).
gate_layer <- function(layer, gating = load_gating()) {
  if (!gating$applied) { layer$gating <- gating; return(layer) }
  n0 <- ncol(layer$M)
  M <- drop_gated(layer$M, gating)
  M <- M[, colSums(M) > 0, drop = FALSE]
  empty <- rownames(M)[rowSums(M) == 0]
  if (length(empty))
    warning(sprintf("gating emptied %d sample(s): %s", length(empty),
                    paste(empty, collapse = ", ")), call. = FALSE)
  layer$M <- M
  layer$gating <- gating
  layer$n_blocked <- n0 - ncol(M)
  layer
}

# Reader-facing label for a table's filter column.
gating_col <- function(gating) if (gating$applied) "contamination-filtered" else "NOT FILTERED"

## ---- reader-facing layer labels ------------------------------------------
# `layer` values written into tables are join keys and keep their names; only the
# displayed label changes. Unmapped values pass through unchanged.
LAYER_DISPLAY <- c(
  "ARG (clinical)"        = "ARG family (clinical)",
  "ARG (all)"             = "ARG family (all)",
  "ARG allele (clinical)" = "ARG variant (clinical)",
  "ARG allele (all)"      = "ARG variant (all)",
  "ARG clinical"          = "ARG family (clinical)",
  "ARG (resistome)"       = "ARG family",
  "ARG (filtered)"        = "ARG family (filtered)",
  "ARG"                   = "ARG family",
  "resistome (ARG)"       = "resistome (ARG)",
  "gut resistome (ARG)"   = "gut resistome (ARG)",
  "KO (function)"         = "KO (function)",
  "function (KO)"         = "KO (function)",
  "gut function (KO)"     = "gut KO (function)",
  "KO"                    = "KO (function)"
)

# Layer keys to display labels. wrap = TRUE breaks long labels over two lines.
layer_display <- function(x, wrap = FALSE) {
  out <- unname(LAYER_DISPLAY[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  # "KO (function)" is short enough to sit on one line; breaking it just makes
  # the strip taller for no gain, so only the two long resistance labels wrap.
  if (wrap) {
    long <- grepl("family|variant", out)
    out[long] <- sub(" \\(", "\n(", out[long])
  }
  out
}

# Ordered factor carrying display labels (family -> variant -> KO by default).
layer_factor <- function(x, levels = NULL, wrap = FALSE) {
  if (is.null(levels)) {
    known <- names(LAYER_DISPLAY)[names(LAYER_DISPLAY) %in% as.character(x)]
    levels <- c(known, setdiff(unique(as.character(x)), known))
  }
  stopifnot(all(as.character(x) %in% levels))
  factor(layer_display(x, wrap = wrap),
         levels = unique(layer_display(levels, wrap = wrap)))
}

# Facet labeller for a layer column that is already ordered: facet_wrap(~layer, labeller = layer_labeller())
layer_labeller <- function(wrap = TRUE) {
  ggplot2::as_labeller(function(x) layer_display(x, wrap = wrap))
}

## ---- reader-facing display names ------------------------------------------
# Join keys (gut, dung, soil; variable names) are never renamed in the data; figures
# and tables draw these display labels. Unmapped values pass through unchanged.
.map_display <- function(x, map) {
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

COMPARTMENT_DISPLAY <- c(gut = "child stool", dung = "cattle dung", soil = "soil")
compartment_display <- function(x) .map_display(x, COMPARTMENT_DISPLAY)

EXPOSURE_DISPLAY <- c(
  cows_head              = "Cattle (head)",
  total_livestock_head   = "Cattle and chickens (head)",
  cows_in_compound       = "Cattle kept in a pen or stable in the compound",
  child_touch_feces      = "Child touches livestock faeces",
  child_antibiotics_2mo  = "Child given antibiotics in the past two months",
  vet_meds_cows          = "Veterinary medicines given to cattle",
  chickens_head          = "Chickens (head)",
  child_eats_soil        = "Child eats soil",
  soap_available         = "Soap available at the handwashing place",
  water_untreated        = "Drinking water untreated",
  hh_leftover_abx        = "Household uses leftover antibiotics",
  ses_income_crops       = "Household income from crops",
  ses_income_ownbusiness = "Household income from own business",
  ses_land_over_1acre    = "Household owns 1 acre or more",
  ses_kitchen_garden     = "Kitchen garden in the compound")
exposure_display <- function(x) .map_display(x, EXPOSURE_DISPLAY)

RESPONSE_DISPLAY <- c(
  arg_abundance = "Clinical ARG abundance per genome",
  arg_richness  = "Clinical ARG richness",
  ko_shannon    = "KO Shannon diversity")
response_display <- function(x) .map_display(x, RESPONSE_DISPLAY)

## ---- per-genome normalization --------------------------------------------
# Single-copy factor per sample = median abundance of the 57 bacterial ribosomal-
# protein marker KOs in the KO table; dividing by it gives copies per average genome.
load_single_copy_markers <- function() {
  f <- file.path(FT_DIR, "single_copy_marker_kos.tsv")
  if (!file.exists(f)) stop("run: Rscript analysis/01_build_single_copy_markers.R")
  read.delim(f, stringsAsFactors = FALSE)$ko
}

# Denominator for the per-genome views:
#   "single_copy"   median of the marker KOs (cancels catalog capture rate)
#   "microbecensus" genome equivalents from microbecensus_ge.tsv
# The two are compared in 06_normalization.R.
NORM_METHOD <- "single_copy"

# Returns list(factor = named per-sample vector, method = label), following
# NORM_METHOD. Falls back to the single-copy marker median with a warning if
# the MicrobeCensus table is requested but absent.
genome_factor <- function(key = load_key(), method = NORM_METHOD) {
  ge_f <- file.path(FT_DIR, "microbecensus_ge.tsv")
  if (method == "microbecensus") {
    if (!file.exists(ge_f)) {
      warning("NORM_METHOD='microbecensus' but microbecensus_ge.tsv is absent; ",
              "falling back to the single-copy marker median", call. = FALSE)
    } else {
      ge  <- read.delim(ge_f, stringsAsFactors = FALSE)
      fac <- setNames(ge$genome_equivalents, ge$sample_id)
      miss <- setdiff(key$sample_id, names(fac))
      if (length(miss))
        warning(sprintf("%d key sample(s) absent from microbecensus_ge.tsv: %s",
                        length(miss), paste(miss, collapse = ", ")), call. = FALSE)
      return(list(factor = fac[intersect(names(fac), key$sample_id)],
                  method = "MicrobeCensus genome equivalents"))
    }
  }
  ko <- get_layer("ko_all", key)
  mk <- intersect(load_single_copy_markers(), colnames(ko$M))
  fac <- apply(ko$M[, mk, drop = FALSE], 1, median)   # per-sample median marker count
  list(factor = fac, method = "single-copy ribosomal-protein median")
}

# Normalize a samples x features matrix to copies-per-genome using a per-sample
# factor (rows divided by the sample's factor). Samples aligned by name.
per_genome <- function(M, gf = genome_factor()) {
  fac <- gf$factor[rownames(M)]
  sweep(M, 1, pmax(fac, .Machine$double.eps), "/")
}

## ---- small helpers ---------------------------------------------------------
# p-value for a figure annotation: decimals up to `dp`, "p < 10^-dp" below that.
fmt_p <- function(p, dp = 4) {
  floor_p <- 10^-dp
  if (!is.finite(p)) return("p = NA")
  if (p < floor_p) return(sprintf("p < %.*f", dp, floor_p))
  s <- sub("\\.$", "", sub("0+$", "", sprintf("%.*f", dp, p)))
  sprintf("p = %s", s)
}

tss <- function(M) sweep(M, 1, pmax(rowSums(M), 1e-9), "/")   # within-sample relative abundance
presence <- function(M, min_count = MIN_COUNT) M > min_count  # boolean presence matrix

# Save a ggplot with the project defaults, returning the path (for logging).
save_fig <- function(plot, file, width = 7, height = 5, dpi = 150) {
  p <- file.path(FIG_DIR, file)
  ggsave(p, plot, width = width, height = height, dpi = dpi)
  p
}
save_tab <- function(df, file) {
  p <- file.path(TAB_DIR, file)
  write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE)
  p
}

