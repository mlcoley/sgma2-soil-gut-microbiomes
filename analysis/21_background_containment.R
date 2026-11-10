#!/usr/bin/env Rscript
# ARG-variant containment against the matched background genes.
# Inputs: T25 sets (19_background_match.R) and feature_tables/background_abundance*.tsv.gz
# (pipeline/12_background-genes/extract_background_abundance.sbatch).
# Outputs: T26a (summary), T26b (per draw).
# Usage: Rscript analysis/21_background_containment.R   (from the repository root)

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
})
source("analysis/00_load_features.R")

set.seed(95)
B_SUB  <- 100                       # equal-effort subsamples
THR    <- 1                         # read threshold
ARMS   <- c("arm1_magnitude", "arm1_shape", "arm2_family")
PRIMARY_ARM <- "arm1_magnitude"
DIRS   <- c("C(dung -> soil)", "C(gut -> soil)", "C(gut -> dung)")

# Matching files read and de-duplicated.
BG_ABUN <- sort(Sys.glob(file.path(FT_DIR, "background_abundance*.tsv.gz")))
if (!length(BG_ABUN))
  stop("no background_abundance*.tsv.gz in ", FT_DIR,
       "\n  built by pipeline/12_background-genes/extract_background_abundance.sbatch")

cat("\n=== ARG vs matched background containment ===\n\n")

## ---- cohort --------------------------------------------------------------

key <- load_key()

## ---- ARG reference layer ---------------------------

GENE_GATING <- load_gating("control_arg_genes_exclusion_list.tsv", "gene-level")
gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")),
                  stringsAsFactors = FALSE)
gab <- gab[gab$sample_id %in% key$sample_id, ]
gab <- gab[!gab$gene_id %in% GENE_GATING$features, ]

hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")),
                   stringsAsFactors = FALSE)
hits <- hits[!duplicated(hits$gene_id), ]
clin_genes <- hits$gene_id[hits$arg_symbol %in% clinical_args()]

samples <- sort(unique(gab$sample_id))
kk <- key[match(samples, key$sample_id), ]
n_eq <- min(table(kk$compartment))
cat(sprintf("cohort %d samples (gut %d, dung %d, soil %d); equal effort n = %d\n",
            length(samples), sum(kk$compartment == "gut"), sum(kk$compartment == "dung"),
            sum(kk$compartment == "soil"), n_eq))

# presence matrices are samples x features logical
build_presence <- function(dt, gene_col, samp_col, cnt_col, genes) {
  M <- matrix(FALSE, nrow = length(samples), ncol = length(genes),
              dimnames = list(samples, genes))
  ok <- dt[[cnt_col]] >= THR & dt[[samp_col]] %in% samples & dt[[gene_col]] %in% genes
  M[cbind(as.character(dt[[samp_col]][ok]), as.character(dt[[gene_col]][ok]))] <- TRUE
  M
}

arg_genes <- sort(intersect(unique(gab$gene_id), clin_genes))
ARG <- build_presence(gab, "gene_id", "sample_id", "count", arg_genes)
cat(sprintf("ARG reference layer: %d clinical alleles\n", ncol(ARG)))

## ---- background presence -------------------------------------------------

bga <- unique(rbindlist(lapply(BG_ABUN, function(f) {
  d <- fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t",
             header = TRUE, showProgress = FALSE)
  cat(sprintf("  %-42s %s rows\n", basename(f), format(nrow(d), big.mark = ",")))
  d
})))
cat(sprintf("background abundance rows (deduplicated): %s\n",
            format(nrow(bga), big.mark = ",")))


need <- unique(unlist(lapply(ARMS, function(a)
  fread(cmd = paste("gzip -dc", shQuote(file.path(TAB_DIR,
        sprintf("T25_%s_sets.tsv.gz", a)))), select = "gene_id",
        showProgress = FALSE)$gene_id)))
absent <- setdiff(need, unique(bga$gene_id))
if (length(absent))
  stop(sprintf("%s of %s matched genes have no abundance rows; extract them before running",
               format(length(absent), big.mark = ","), format(length(need), big.mark = ",")))
cat(sprintf("all %s matched genes have abundance rows\n", format(length(need), big.mark = ",")))
bg_genes <- sort(unique(bga$gene_id))
BG <- build_presence(bga, "gene_id", "sample_id", "count", bg_genes)
cat(sprintf("background presence matrix: %d samples x %s genes\n",
            nrow(BG), format(ncol(BG), big.mark = ",")))
rm(bga); invisible(gc())

## ---- containment machinery -------------------------

containment <- function(A, B) if (length(A)) length(intersect(A, B)) / length(A) else NA_real_

# exact hypergeometric null: E[C] = K/N, sd in closed form
contain_with_null <- function(S, A, B, label) {
  N <- length(unique(unlist(S))); K <- length(B); n <- length(A)
  obs <- containment(A, B)
  mu  <- K / N
  sdC <- if (n > 0 && N > 1) sqrt(n * mu * (1 - mu) * (N - n) / (N - 1)) / n else NA_real_
  data.frame(direction = label, observed = obs, null_mean = mu,
             SES = if (is.na(sdC) || sdC == 0) NA_real_ else (obs - mu) / sdC)
}

comp_idx <- lapply(setNames(COMP_LEVELS, COMP_LEVELS),
                   function(cc) which(kk$compartment == cc))

run_layer <- function(M) {
  reps <- lapply(seq_len(B_SUB), function(b) {
    si <- lapply(comp_idx, function(i) sample(i, n_eq))
    S <- lapply(si, function(i) which(colSums(M[i, , drop = FALSE]) > 0))
    rbind(
      contain_with_null(S, S$dung, S$soil, DIRS[1]),
      contain_with_null(S, S$gut,  S$soil, DIRS[2]),
      contain_with_null(S, S$gut,  S$dung, DIRS[3]))
  })
  do.call(rbind, reps) %>% group_by(direction) %>%
    summarise(observed = mean(observed), null_mean = mean(null_mean),
              SES = mean(SES), .groups = "drop")
}

## ---- ARG reference value -------------------------------------------------

cat("\ncomputing ARG reference ...\n")
arg_res <- run_layer(ARG) %>% mutate(n_features = ncol(ARG))
print(as.data.frame(arg_res), row.names = FALSE)

## ---- background draws ----------------------------------------------------

per_draw <- list()
for (arm in ARMS) {
  f <- file.path(TAB_DIR, sprintf("T25_%s_sets.tsv.gz", arm))
  sets <- fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t", showProgress = FALSE)
  draws <- split(sets$gene_id, sets$draw)
  cat(sprintf("\n%s: %d draws x %d genes\n", arm, length(draws), length(draws[[1]])))
  res <- vector("list", length(draws))
  for (b in seq_along(draws)) {
    g <- intersect(draws[[b]], colnames(BG))
    res[[b]] <- run_layer(BG[, g, drop = FALSE]) %>%
      mutate(draw = b, arm = arm, n_features = length(g))
    if (b %% 25 == 0) cat(sprintf("  draw %d/%d\n", b, length(draws)))
  }
  per_draw[[arm]] <- bind_rows(res)
}
per_draw <- bind_rows(per_draw)

## ---- locate the ARG value in each background distribution ---------------
# empirical two-sided p

summ <- per_draw %>%
  left_join(arg_res %>% select(direction, arg_SES = SES, arg_obs = observed,
                               arg_null = null_mean, arg_n = n_features),
            by = "direction") %>%
  group_by(arm, direction) %>%
  summarise(
    n_draw      = n(),
    bg_SES_mean = mean(SES), bg_SES_sd = sd(SES),
    bg_obs_mean = mean(observed), bg_obs_sd = sd(observed),
    bg_n_feat   = mean(n_features),
    arg_SES = first(arg_SES), arg_obs = first(arg_obs),
    arg_null = first(arg_null), arg_n = first(arg_n),
    # standardised position of the ARG value within the background envelope
    z_vs_bg = (first(arg_SES) - mean(SES)) / sd(SES),
    p_emp   = 2 * min(mean(SES <= first(arg_SES)), mean(SES >= first(arg_SES))),
    .groups = "drop") %>%
  mutate(p_emp = pmin(1, pmax(p_emp, 1 / n_draw)))

# BH across the three compartment pairs of the primary arm only
summ$q_BH <- NA_real_
i <- summ$arm == PRIMARY_ARM
summ$q_BH[i] <- p.adjust(summ$p_emp[i], method = "BH")

## ---- feature-count guard -------------------------------------------------
bad <- summ %>% filter(abs(bg_n_feat - arg_n) > 0.5)
if (nrow(bad)) {
  cat("\nfeature counts differ; these rows are not SES-comparable (read observed containment):\n")
  print(as.data.frame(bad %>% select(arm, direction, bg_n_feat, arg_n)), row.names = FALSE)
}

## ---- report --------------------------------------------------------------

fmt <- function(d) {
  d %>% transmute(arm, direction,
                  ARG_SES = round(arg_SES, 2),
                  BG_SES = sprintf("%.2f ± %.2f", bg_SES_mean, bg_SES_sd),
                  z_vs_bg = round(z_vs_bg, 2),
                  p_emp = signif(p_emp, 3), q_BH = ifelse(is.na(q_BH), "", signif(q_BH, 3)),
                  ARG_C = round(arg_obs, 3),
                  BG_C = sprintf("%.3f ± %.3f", bg_obs_mean, bg_obs_sd)) %>% as.data.frame()
}
for (arm in ARMS) {
  cat(sprintf("\n--- %s%s ---\n", arm,
              if (arm == PRIMARY_ARM) "  (primary, BH across 3 directions)" else "  (sensitivity, raw values only)"))
  print(fmt(summ %>% filter(arm == !!arm)), row.names = FALSE)
}

cat("\n--- primary arm, one row per compartment pair ---\n")
for (d in DIRS) {
  hl <- summ %>% filter(arm == PRIMARY_ARM, direction == d)
  cat(sprintf("%s, clinical variants, %d features:\n", d, hl$arg_n))
  cat(sprintf("  ARG        SES %+.2f   containment %.3f\n", hl$arg_SES, hl$arg_obs))
  cat(sprintf("  background SES %+.2f ± %.2f   containment %.3f ± %.3f  (%d draws)\n",
              hl$bg_SES_mean, hl$bg_SES_sd, hl$bg_obs_mean, hl$bg_obs_sd, hl$n_draw))
  cat(sprintf("  ARG sits %+.2f sd from the background mean; empirical p = %.3f, BH-adjusted p = %.3f\n",
              hl$z_vs_bg, hl$p_emp, hl$q_BH))
}

save_tab(as.data.frame(summ), "T26a_background_containment.tsv")
save_tab(as.data.frame(per_draw), "T26b_background_per_draw.tsv")
cat("\nBackground containment complete.\n")
