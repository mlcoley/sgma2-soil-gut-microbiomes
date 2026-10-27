#!/usr/bin/env Rscript
# Within-household excess similarity restricted to rare features.
# Whole-profile household tests are dominated by cosmopolitan features that cannot
# individuate a household, so the test is repeated on feature spaces capped at a
# series of household prevalences (100%, 50%, 25%, 10%). Same statistic and null at
# every cap; only the feature space changes. Prevalence is measured over households
# (a feature is present in a household if it occurs in any of its compartments).
# For each layer, cap, metric (Jaccard on presence, Bray-Curtis on abundance) and
# compartment pair, the mean within-household similarity (diagonal of the household x
# household matrix) is compared with 999 permutations of the household labels. BH is
# applied across the whole sweep (120 tests). Restricting to rare features lowers the
# feature count and widens the null, so the minimum detectable effect (mde_pct) is
# reported for every cell. Every cell reaching BH-adjusted p < 0.05 gets an influence
# check (T13c): the leave-one-household-out test for the household supplying the most
# of the excess, and for H016. The null is rebuilt without the dropped household.
# Outputs: T13 (sweep), T13b (power trade-off), T13c (influence), Figures 5 and S3.
# Usage: Rscript analysis/07_rare_feature_sharing.R   (from the repository root)

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(vegan) })
source("analysis/00_load_features.R")
source("analysis/util_fig_rare_sharing.R")

set.seed(1)
N_PERM     <- 999
THR        <- 1                       # detection threshold (reads) for presence
PREV_CAPS  <- c(1.00, 0.50, 0.25, 0.10)
PAIR_TYPES <- list(c("gut", "dung"), c("gut", "soil"), c("dung", "soil"))

log_lines <- c(sprintf("Rare-feature household sharing - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

key    <- load_key()
GATING <- load_gating()
arg    <- gate_layer(get_layer("arg", key), GATING)
ko     <- get_layer("ko_all", key)
clin   <- clinical_args()

# ARG variant layers: gene-level exclusion list (never the symbol-level list)
GENE_GATING <- load_gating("control_arg_genes_blocklist.tsv", "gene-level")
gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")), stringsAsFactors = FALSE)
gab <- gab[gab$sample_id %in% key$sample_id, ]
g_samp <- sort(unique(gab$sample_id)); g_gene <- sort(unique(gab$gene_id))
ALLELE <- matrix(0, length(g_samp), length(g_gene), dimnames = list(g_samp, g_gene))
ALLELE[cbind(gab$sample_id, gab$gene_id)] <- gab$count
ALLELE <- ALLELE[, !colnames(ALLELE) %in% GENE_GATING$features, drop = FALSE]
hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")), stringsAsFactors = FALSE)
gene2sym <- setNames(hits$arg_symbol, hits$gene_id)
clin_genes <- names(gene2sym)[gene2sym %in% clin]

LAYERS <- list(
  "ARG (clinical)"        = arg$M[, colnames(arg$M) %in% clin, drop = FALSE],
  "ARG (all)"             = arg$M,
  "ARG allele (clinical)" = ALLELE[, colnames(ALLELE) %in% clin_genes, drop = FALSE],
  "ARG allele (all)"      = ALLELE,
  "KO (function)"         = ko$M)

## ---- similarity and permutation machinery ---------------------------------
jaccard <- function(P) {
  inter <- tcrossprod(P * 1); sz <- rowSums(P)
  J <- inter / pmax(outer(sz, sz, "+") - inter, 1)
  dimnames(J) <- list(rownames(P), rownames(P)); J
}
bray_sim <- function(M) {
  S <- 1 - as.matrix(vegdist(tss(M), "bray")); S[!is.finite(S)] <- 0; S
}
hh_matrix <- function(S, kk, c1, c2) {
  i1 <- which(kk$compartment == c1); i2 <- which(kk$compartment == c2)
  h1 <- as.character(kk$household[i1]); h2 <- as.character(kk$household[i2])
  hs <- sort(intersect(unique(h1), unique(h2)))
  if (length(hs) < 5) return(NULL)
  M <- matrix(NA_real_, length(hs), length(hs), dimnames = list(hs, hs))
  for (a in seq_along(hs)) {
    ra <- i1[h1 == hs[a]]
    for (b in seq_along(hs)) M[a, b] <- mean(S[ra, i2[h2 == hs[b]], drop = FALSE])
  }
  M
}
perm_test <- function(M, n_perm = N_PERM) {
  n <- nrow(M); obs <- mean(diag(M))
  nul <- replicate(n_perm, mean(M[cbind(seq_len(n), sample.int(n))]))
  list(n_households = n, observed = obs, null_mean = mean(nul), null_sd = sd(nul),
       z = (obs - mean(nul)) / sd(nul), ratio = obs / mean(nul),
       mde_pct = 100 * 1.645 * sd(nul) / mean(nul),
       p = (1 + sum(nul >= obs)) / (n_perm + 1))
}

# Household prevalence of each feature: in how many households does it appear in
# at least one compartment?
hh_prevalence <- function(M, kk) {
  P <- M >= THR
  hh <- as.character(kk$household)
  agg <- rowsum(P * 1, hh) > 0            # households x features
  colMeans(agg)
}

## ---- the sweep -----------------------------------------------------------
logmsg("== rare-feature sweep ==")
logmsg("prevalence measured over %d households; caps: %s",
       length(unique(key$household)), paste0(100 * PREV_CAPS, "%", collapse = ", "))

res <- list()
for (ln in names(LAYERS)) {
  Mx <- LAYERS[[ln]]
  kk <- key[match(rownames(Mx), key$sample_id), ]
  prev <- hh_prevalence(Mx, kk)
  logmsg("")
  logmsg("-- %s: %d features; household prevalence median %.2f, %.0f%% of features in >75%% of households",
         ln, ncol(Mx), median(prev), 100 * mean(prev > 0.75))
  for (cap in PREV_CAPS) {
    keepf <- names(prev)[prev <= cap]
    if (length(keepf) < 10) {
      logmsg("   cap <=%3.0f%%: only %d features — SKIPPED (below the 10-feature floor)",
             100 * cap, length(keepf)); next
    }
    Mk <- Mx[, keepf, drop = FALSE]
    # Samples emptied by the restriction would give undefined similarities.
    empt <- sum(rowSums(Mk >= THR) == 0)
    for (metric in c("Jaccard (presence)", "Bray-Curtis (abundance)")) {
      S <- if (metric == "Jaccard (presence)") jaccard(Mk >= THR) else bray_sim(Mk)
      for (pr in PAIR_TYPES) {
        HM <- hh_matrix(S, kk, pr[1], pr[2]); if (is.null(HM)) next
        pt <- perm_test(HM)
        res[[length(res) + 1]] <- data.frame(
          layer = ln, prev_cap = cap, n_features = length(keepf),
          pct_features_kept = round(100 * length(keepf) / ncol(Mx), 1),
          n_samples_emptied = empt, metric = metric,
          pair = paste(pr, collapse = "-"),
          n_households = pt$n_households, within = pt$observed,
          null_mean = pt$null_mean, z = pt$z, ratio = pt$ratio,
          mde_pct = pt$mde_pct, p = pt$p, stringsAsFactors = FALSE)
      }
    }
    logmsg("   cap <=%3.0f%%: %5d features (%4.1f%% kept), %d samples emptied",
           100 * cap, length(keepf), 100 * length(keepf) / ncol(Mx), empt)
  }
}
t13 <- do.call(rbind, res)
# BH across the whole sweep.
t13$q_BH <- p.adjust(t13$p, method = "BH")
t13$contamination_filter <- gating_col(GATING)
t13[] <- lapply(t13, function(x) if (is.numeric(x)) round(x, 4) else x)
save_tab(t13, "T13_rare_feature_sharing.tsv")

logmsg("")
logmsg("== RESULT ==")
logmsg("%d tests across %d layers x %d caps x 2 metrics x 3 pairs; BH across the whole sweep",
       nrow(t13), length(LAYERS), length(PREV_CAPS))
logmsg("reaching BH-adjusted p < 0.05: %d", sum(t13$q_BH < 0.05))
top <- head(t13[order(t13$p), ], 10)
for (i in seq_len(nrow(top)))
  logmsg("  %-21s cap<=%3.0f%% %-23s %-9s z=%+5.2f ratio=%.3f p=%.3f BH-adjusted p=%.3f (MDE %.1f%%, %d feat)",
         top$layer[i], 100 * top$prev_cap[i], top$metric[i], top$pair[i],
         top$z[i], top$ratio[i], top$p[i], top$q_BH[i], top$mde_pct[i], top$n_features[i])

logmsg("")
logmsg("MDE by cap (median across tests) — the cost of restricting:")
mde_by_cap <- t13 %>% group_by(prev_cap) %>%
  summarise(median_mde_pct = round(median(mde_pct), 2),
            median_abs_z = round(median(abs(z)), 2),
            max_z = round(max(z), 2), .groups = "drop") %>% as.data.frame()
for (i in seq_len(nrow(mde_by_cap)))
  logmsg("  cap <=%3.0f%%: median MDE %6.2f%% | median |z| %.2f | max z %+.2f",
         100 * mde_by_cap$prev_cap[i], mde_by_cap$median_mde_pct[i],
         mde_by_cap$median_abs_z[i], mde_by_cap$max_z[i])
save_tab(mde_by_cap, "T13b_rare_feature_power_tradeoff.tsv")

# ===========================================================================
# Influence check on every cell reaching BH-adjusted p < 0.05 -> T13c
# ===========================================================================
# H016's dung sample is the most soil-like in the cohort. A leave-one-out on a
# permutation statistic must drop the household from the null as well as from the
# observed statistic, so the household x household matrix is subset (not just the
# diagonal entry skipped).
logmsg("")
logmsg("== influence check on significant cells ==")
sig <- t13[t13$q_BH < 0.05, ]
infl <- list()
if (nrow(sig)) {
  for (i in seq_len(nrow(sig))) {
    ln <- sig$layer[i]; cap <- sig$prev_cap[i]
    Mx <- LAYERS[[ln]]; kk <- key[match(rownames(Mx), key$sample_id), ]
    prev <- hh_prevalence(Mx, kk)
    Mk <- Mx[, names(prev)[prev <= cap], drop = FALSE]
    S  <- if (sig$metric[i] == "Jaccard (presence)") jaccard(Mk >= THR) else bray_sim(Mk)
    pr <- strsplit(sig$pair[i], "-")[[1]]
    HM <- hh_matrix(S, kk, pr[1], pr[2])
    full <- perm_test(HM)
    # per-household share of the excess over the null
    exc <- diag(HM) - full$null_mean
    shr <- exc / sum(exc)
    worst <- names(sort(shr, decreasing = TRUE))[1]
    drop_one <- function(h) {
      keep <- setdiff(rownames(HM), h)
      perm_test(HM[keep, keep, drop = FALSE])
    }
    d16 <- if ("H016" %in% rownames(HM)) drop_one("H016") else NULL
    dw  <- drop_one(worst)
    infl[[length(infl) + 1]] <- data.frame(
      layer = ln, prev_cap = cap, metric = sig$metric[i], pair = sig$pair[i],
      z_full = round(full$z, 2), p_full = full$p,
      top_household = worst, top_share_of_excess = round(100 * max(shr), 1),
      z_drop_top = round(dw$z, 2), p_drop_top = dw$p,
      z_drop_H016 = if (is.null(d16)) NA_real_ else round(d16$z, 2),
      p_drop_H016 = if (is.null(d16)) NA_real_ else d16$p,
      stringsAsFactors = FALSE)
  }
  infl <- do.call(rbind, infl)
  save_tab(infl, "T13c_rare_feature_influence.tsv")
  for (i in seq_len(nrow(infl)))
    logmsg("  %-21s cap<=%3.0f%% %-23s %-9s z=%+.2f | top hh %s supplies %.1f%% of excess -> z=%+.2f (p=%.3f) | drop H016 -> z=%+.2f (p=%.3f)",
           infl$layer[i], 100 * infl$prev_cap[i], infl$metric[i], infl$pair[i],
           infl$z_full[i], infl$top_household[i], infl$top_share_of_excess[i],
           infl$z_drop_top[i], infl$p_drop_top[i], infl$z_drop_H016[i], infl$p_drop_H016[i])
} else logmsg("  (no cell reached BH-adjusted p < 0.05)")

save_fig(f17a_main(t13), "F17a_rare_feature_sharing.png", width = 9.4, height = 5.3)
save_fig(f17b_supplement(t13), "F17b_rare_feature_sharing_all_layers.png",
         width = 9.6, height = 8.4)

writeLines(unlist(log_lines), file.path(LOG_DIR, "07_rare_feature_sharing.log"))
cat("\nRare-feature household sharing complete.\n")
