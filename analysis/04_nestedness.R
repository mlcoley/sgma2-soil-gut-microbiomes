#!/usr/bin/env Rscript
# Directional containment of child-stool features in the environment (cattle dung +
# soil) at equal sampling effort, with the hypergeometric null.
#
# A compartment's feature set is a union over its samples and grows with sample
# number (soil 117 samples, stool 39), so every metric is computed on equal-sized
# subsamples (n = the smallest compartment) over B_SUB draws and reported as
# mean +/- sd. The null draws a source-sized set at random from the union of the
# three compartments; containment then follows a hypergeometric distribution,
# E[C(A->B)] = |B| / N, so the SES uses the exact moments rather than a simulation.
# Outputs: T8a (equal-effort richness), T8b (containment and SES), T8f (gut and
# environment feature counts), Figure 3.
# Usage: Rscript analysis/04_nestedness.R   (from the repository root)

suppressPackageStartupMessages({
  library(vegan); library(ggplot2); library(dplyr); library(tidyr); library(scales)
})
source("analysis/00_load_features.R")
source("analysis/util_fig_containment.R")

set.seed(1)
B_SUB      <- 100     # equal-effort subsamples per layer x threshold
THRESHOLDS <- c(1, 5) # presence = count >= threshold

log_lines <- c(sprintf("Containment - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines) + 1]] <<- m }

theme_sgma <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

key <- load_key()
ann <- load_arg_annotation()
clin <- clinical_args(ann)
arg <- get_layer("arg", key)
ko  <- get_layer("ko_all", key)

# ARG-family layers use the symbol-level exclusion list; gene-level (variant)
# layers use the gene-level list. Never apply one list to the other level.
GATING <- load_gating()
logmsg("filter: %s", GATING$note)
GENE_GATING <- load_gating("control_arg_genes_blocklist.tsv", "gene-level")

# ARG variant layers: catalog representative genes (95% identity clusters).
gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")), stringsAsFactors = FALSE)
gab <- gab[gab$sample_id %in% key$sample_id, ]
g_samp <- sort(unique(gab$sample_id)); g_gene <- sort(unique(gab$gene_id))
ALLELE <- matrix(0, nrow = length(g_samp), ncol = length(g_gene),
                 dimnames = list(g_samp, g_gene))
ALLELE[cbind(gab$sample_id, gab$gene_id)] <- gab$count
ALLELE <- ALLELE[, !colnames(ALLELE) %in% GENE_GATING$features, drop = FALSE]
hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")), stringsAsFactors = FALSE)
gene2sym <- setNames(hits$arg_symbol, hits$gene_id)
clin_genes <- names(gene2sym)[gene2sym %in% clin]
logmsg("variant layer: %d samples x %d catalog genes (%s); clinical variants %d",
       nrow(ALLELE), ncol(ALLELE), GENE_GATING$note,
       sum(colnames(ALLELE) %in% clin_genes))

LAYERS_OQ2 <- list(
  "ARG (clinical)"        = drop_gated(arg$M[, colnames(arg$M) %in% clin, drop = FALSE], GATING),
  "ARG (all)"             = drop_gated(arg$M, GATING),
  "ARG allele (clinical)" = ALLELE[, colnames(ALLELE) %in% clin_genes, drop = FALSE],
  "ARG allele (all)"      = ALLELE,
  "KO (function)"         = ko$M
)
stopifnot(identical(names(LAYERS_OQ2), LAYER_ORDER_OQ2))
logmsg("layers: %s", paste(sprintf("%s %d features", names(LAYERS_OQ2),
       vapply(LAYERS_OQ2, ncol, 1L)), collapse = " | "))

# C(A->B) = |A n B| / |A|
containment <- function(A, B) if (length(A)) length(intersect(A, B)) / length(A) else NA_real_
overlap_coef <- function(A, B) {
  m <- min(length(A), length(B)); if (m) length(intersect(A, B)) / m else NA_real_
}

# ===========================================================================
# Compartment feature sets under controlled sampling effort -> T8a
# ===========================================================================
# A compartment's feature set = features present (count >= thr) in >= 1 of its samples.
comp_sets <- function(M, kk, thr, samp_idx = NULL) {
  lapply(setNames(COMP_LEVELS, COMP_LEVELS), function(cc) {
    i <- if (is.null(samp_idx)) which(kk$compartment == cc) else samp_idx[[cc]]
    if (!length(i)) return(integer(0))
    which(colSums(M[i, , drop = FALSE] >= thr) > 0)
  })
}
draw_equal <- function(kk, n_eq) lapply(setNames(COMP_LEVELS, COMP_LEVELS),
                                       function(cc) sample(which(kk$compartment == cc), n_eq))

logmsg("== compartment feature sets ==")
t8a <- list()
for (ln in names(LAYERS_OQ2)) {
  M  <- LAYERS_OQ2[[ln]]
  kk <- key[match(rownames(M), key$sample_id), ]
  n_by <- table(kk$compartment)[COMP_LEVELS]
  n_eq <- min(n_by)
  for (thr in THRESHOLDS) {
    full <- comp_sets(M, kk, thr)
    rich <- replicate(B_SUB, vapply(comp_sets(M, kk, thr, draw_equal(kk, n_eq)), length, 1L))
    t8a[[length(t8a) + 1]] <- data.frame(
      layer = ln, threshold = thr, compartment = COMP_LEVELS,
      n_samples_all = as.integer(n_by), n_samples_equal = n_eq,
      richness_all = vapply(full, length, 1L)[COMP_LEVELS],
      richness_equal_mean = round(rowMeans(rich)[COMP_LEVELS], 1),
      richness_equal_sd = round(apply(rich, 1, sd)[COMP_LEVELS], 2), row.names = NULL)
  }
}
t8a <- do.call(rbind, t8a)
t8a$inflation_from_effort <- round(t8a$richness_all / pmax(t8a$richness_equal_mean, 1e-9), 3)
save_tab(t8a, "T8a_compartment_richness.tsv")

# ===========================================================================
# Directional containment and null SES -> T8b
# ===========================================================================
contain_with_null <- function(S, A, B, label) {
  N <- length(unique(unlist(S))); K <- length(B); n <- length(A)
  obs <- containment(A, B)
  mu  <- K / N
  sdC <- if (n > 0 && N > 1) sqrt(n * mu * (1 - mu) * (N - n) / (N - 1)) / n else NA_real_
  data.frame(direction = label, observed = obs, null_mean = mu,
             SES = if (is.na(sdC) || sdC == 0) NA_real_ else (obs - mu) / sdC)
}

logmsg("== directional containment ==")
t8b <- list()
for (ln in names(LAYERS_OQ2)) {
  M  <- LAYERS_OQ2[[ln]]
  kk <- key[match(rownames(M), key$sample_id), ]
  n_eq <- min(table(kk$compartment))
  for (thr in THRESHOLDS) {
    reps <- lapply(seq_len(B_SUB), function(b) {
      S   <- comp_sets(M, kk, thr, draw_equal(kk, n_eq))
      env <- union(S$dung, S$soil)
      rbind(
        contain_with_null(S, S$gut,  env,    "C(gut -> env)"),
        contain_with_null(S, env,    S$gut,  "C(env -> gut)"),
        contain_with_null(S, S$dung, S$soil, "C(dung -> soil)"),
        contain_with_null(S, S$gut,  S$soil, "C(gut -> soil)"),
        contain_with_null(S, S$gut,  S$dung, "C(gut -> dung)"),
        contain_with_null(S, S$soil, S$gut,  "C(soil -> gut)"),
        data.frame(direction = "overlap_coef(gut, env)",
                   observed = overlap_coef(S$gut, env),
                   null_mean = NA_real_, SES = NA_real_))
    })
    agg <- do.call(rbind, reps) %>% group_by(direction) %>%
      summarise(observed_mean = mean(observed), observed_sd = sd(observed),
                null_mean = mean(null_mean), SES_mean = mean(SES),
                SES_sd = sd(SES), .groups = "drop")
    agg$layer <- ln; agg$threshold <- thr; agg$n_equal <- n_eq
    t8b[[length(t8b) + 1]] <- agg
  }
}
t8b <- do.call(rbind, t8b) %>%
  select(layer, threshold, n_equal, direction, observed_mean, observed_sd,
         null_mean, SES_mean, SES_sd)
t8b[] <- lapply(t8b, function(x) if (is.numeric(x)) round(x, 4) else x)
t8b$contamination_filter <- gating_col(GATING)
save_tab(t8b, "T8b_directional_containment.tsv")
for (ln in names(LAYERS_OQ2)) {
  r <- t8b[t8b$layer == ln & t8b$threshold == 1, ]
  g2e <- r[r$direction == "C(gut -> env)", ]; e2g <- r[r$direction == "C(env -> gut)", ]
  logmsg("containment [%s] C(gut->env)=%.3f (null %.3f, SES %+.2f) vs C(env->gut)=%.3f (null %.3f, SES %+.2f)",
         ln, g2e$observed_mean, g2e$null_mean, g2e$SES_mean,
         e2g$observed_mean, e2g$null_mean, e2g$SES_mean)
}

# ---------------------------------------------------------------------------
# Counts of features shared by, and private to, stool and the environment -> T8f
# (equal effort, mean and sd over the same B_SUB draws)
# ---------------------------------------------------------------------------
logmsg("== stool / environment feature counts ==")
t8f <- list()
for (ln in names(LAYERS_OQ2)) {
  M  <- LAYERS_OQ2[[ln]]
  kk <- key[match(rownames(M), key$sample_id), ]
  n_eq <- min(table(kk$compartment))
  for (thr in THRESHOLDS) {
    reps <- do.call(rbind, lapply(seq_len(B_SUB), function(b) {
      S   <- comp_sets(M, kk, thr, draw_equal(kk, n_eq))
      env <- union(S$dung, S$soil)
      data.frame(gut_total = length(S$gut), env_total = length(env),
                 shared = length(intersect(S$gut, env)),
                 gut_private = length(setdiff(S$gut, env)),
                 env_private = length(setdiff(env, S$gut)))
    }))
    t8f[[length(t8f) + 1]] <- data.frame(
      layer = ln, threshold = thr, n_equal = n_eq,
      gut_total    = mean(reps$gut_total),    gut_total_sd    = sd(reps$gut_total),
      env_total    = mean(reps$env_total),    env_total_sd    = sd(reps$env_total),
      shared       = mean(reps$shared),       shared_sd       = sd(reps$shared),
      gut_private  = mean(reps$gut_private),  gut_private_sd  = sd(reps$gut_private),
      env_private  = mean(reps$env_private),  env_private_sd  = sd(reps$env_private),
      pct_env_private = 100 * mean(reps$env_private) / mean(reps$env_total))
  }
}
t8f <- do.call(rbind, t8f)
t8f[] <- lapply(t8f, function(x) if (is.numeric(x)) round(x, 1) else x)
t8f$contamination_filter <- gating_col(GATING)
save_tab(t8f, "T8f_private_repertoire_partition.tsv")
for (i in seq_len(nrow(t8f)))
  logmsg("feature counts [%s thr>=%d]: gut %.0f (%.0f private) | env %.0f (%.0f private, %.1f%%) | shared %.0f",
         t8f$layer[i], t8f$threshold[i], t8f$gut_total[i], t8f$gut_private[i],
         t8f$env_total[i], t8f$env_private[i], t8f$pct_env_private[i], t8f$shared[i])

# Figure 3
save_fig(f12a_containment(t8b, t8f), "F12a_containment.png", width = 10.5, height = 4.6)

stopifnot(abs(containment(1:5, 1:10) - 1) < 1e-12,
          abs(containment(1:10, 1:5) - 0.5) < 1e-12)

writeLines(log_lines, file.path(LOG_DIR, "04_nestedness.log"))
cat("\nContainment complete.\n")
