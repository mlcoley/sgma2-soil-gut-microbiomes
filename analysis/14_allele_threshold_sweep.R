#!/usr/bin/env Rscript
# ARG-variant containment and household sharing against the catalog identity threshold.
# Input: feature_tables/arg_reclust/arg_id{95,90,85,80}_cluster.tsv (representative, member).
# Outputs: T16a (merge summary), T16b (containment vs threshold), T16c (sharing vs
# threshold).
# Usage: Rscript analysis/14_allele_threshold_sweep.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(vegan) })
source("analysis/00_load_features.R")

set.seed(78)
N_PERM <- 999
THR    <- 1
B_SUB  <- 100     # equal-effort subsamples, as in 04_nestedness.R

log_lines <- c(sprintf("ARG variant threshold sweep - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

CLUST_DIR <- file.path(FT_DIR, "arg_reclust")
cf <- if (dir.exists(CLUST_DIR)) sort(list.files(CLUST_DIR, "^arg_id[0-9]+_cluster\\.tsv$",
                                                 full.names = TRUE)) else character(0)
if (!length(cf))
  stop("no cluster tables in ", CLUST_DIR,
       "; run pipeline/11_feature-tables/arg_recluster_sweep.sbatch")

key  <- load_key()
clin <- clinical_args()
GENE_GATING <- load_gating("control_arg_genes_exclusion_list.tsv", "gene-level")
gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")), stringsAsFactors = FALSE)
gab <- gab[gab$sample_id %in% key$sample_id, ]
gab <- gab[!gab$gene_id %in% GENE_GATING$features, ]
hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")), stringsAsFactors = FALSE)
gene2sym <- setNames(hits$arg_symbol, hits$gene_id)

## ---- machinery (as in 04_nestedness.R and 07_rare_feature_sharing.R) ------
containment <- function(A, B) if (length(A)) length(intersect(A, B)) / length(A) else NA_real_
contain_ses <- function(S, A, B) {
  N <- length(unique(unlist(S))); K <- length(B); n <- length(A)
  if (!n || !N) return(c(NA, NA, NA))
  mu <- K / N
  sdC <- if (N > 1) sqrt(n * mu * (1 - mu) * (N - n) / (N - 1)) / n else NA_real_
  obs <- containment(A, B)
  c(obs, mu, if (is.na(sdC) || sdC == 0) NA_real_ else (obs - mu) / sdC)
}
draw_equal <- function(kk, n_eq) unlist(lapply(COMP_LEVELS, function(cc) {
  i <- which(kk$compartment == cc); sample(i, min(n_eq, length(i)))
}))
comp_sets <- function(M, kk, idx) {
  lapply(setNames(COMP_LEVELS, COMP_LEVELS), function(cc) {
    j <- idx[kk$compartment[idx] == cc]
    if (!length(j)) character(0) else colnames(M)[colSums(M[j, , drop = FALSE] >= THR) > 0]
  })
}
jaccard <- function(P) {
  inter <- tcrossprod(P * 1); sz <- rowSums(P)
  J <- inter / pmax(outer(sz, sz, "+") - inter, 1)
  dimnames(J) <- list(rownames(P), rownames(P)); J
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
  list(observed = obs, null_mean = mean(nul), z = (obs - mean(nul)) / sd(nul),
       ratio = obs / mean(nul), mde_pct = 100 * 1.645 * sd(nul) / mean(nul),
       p = (1 + sum(nul >= obs)) / (n_perm + 1), n_households = n)
}

## ---- sweep ---------------------------------------------------------------
res_cont <- list(); res_share <- list(); res_merge <- list()
for (f in cf) {
  tag <- sub(".*arg_id([0-9]+)_cluster\\.tsv$", "\\1", f)
  idpct <- as.numeric(tag)
  cl <- read.delim(f, header = FALSE, col.names = c("rep", "member"), stringsAsFactors = FALSE)
  map <- setNames(cl$rep, cl$member)

  g <- gab
  g$feature <- map[g$gene_id]
  n_unmapped <- sum(is.na(g$feature))
  if (n_unmapped) {
    logmsg("id=%s%%: %d abundance rows have no cluster assignment; kept as their own feature",
           tag, n_unmapped)
    g$feature[is.na(g$feature)] <- g$gene_id[is.na(g$feature)]
  }
  agg <- g %>% group_by(sample_id, feature) %>%
    summarise(count = sum(count), .groups = "drop")
  samps <- sort(unique(agg$sample_id)); feats <- sort(unique(agg$feature))
  M <- matrix(0, length(samps), length(feats), dimnames = list(samps, feats))
  M[cbind(agg$sample_id, agg$feature)] <- agg$count
  kk <- key[match(rownames(M), key$sample_id), ]

  # clinical subset: a merged cluster is clinical if ANY member gene is
  clin_feats <- unique(map[names(gene2sym)[gene2sym %in% clin]])
  clin_feats <- clin_feats[!is.na(clin_feats)]

  # how much merging happened overall and per ARG symbol
  n_in  <- length(unique(gab$gene_id)); n_out <- ncol(M)
  res_merge[[length(res_merge) + 1]] <- data.frame(
    identity_pct = idpct, n_genes_in = n_in, n_features_out = n_out,
    merge_ratio = round(n_out / n_in, 3),
    n_clinical_features = length(intersect(clin_feats, colnames(M))))
  logmsg("id=%s%%: %d genes -> %d features (%.1f%% retained); clinical features %d",
         tag, n_in, n_out, 100 * n_out / n_in, length(intersect(clin_feats, colnames(M))))

  for (layer in c("clinical", "all")) {
    Mk <- if (layer == "clinical") M[, colnames(M) %in% clin_feats, drop = FALSE] else M
    if (ncol(Mk) < 20) next
    n_eq <- min(table(kk$compartment))

    # ---- containment, equal effort + exact hypergeometric null ------------
    reps <- lapply(seq_len(B_SUB), function(b) {
      S <- comp_sets(Mk, kk, draw_equal(kk, n_eq))
      rbind(data.frame(direction = "C(gut -> soil)",  t(contain_ses(S, S$gut, S$soil))),
            data.frame(direction = "C(dung -> soil)", t(contain_ses(S, S$dung, S$soil))),
            data.frame(direction = "C(gut -> dung)",  t(contain_ses(S, S$gut, S$dung))))
    })
    agg2 <- do.call(rbind, reps)
    names(agg2)[2:4] <- c("observed", "null_mean", "SES")
    agg2 <- agg2 %>% group_by(direction) %>%
      summarise(observed = mean(observed, na.rm = TRUE),
                null_mean = mean(null_mean, na.rm = TRUE),
                SES = mean(SES, na.rm = TRUE), SES_sd = sd(SES, na.rm = TRUE), .groups = "drop")
    agg2$identity_pct <- idpct; agg2$layer <- layer; agg2$n_features <- ncol(Mk)
    res_cont[[length(res_cont) + 1]] <- as.data.frame(agg2)

    # ---- household sharing --------------------------------------------------
    S <- jaccard(Mk >= THR)
    for (pr in list(c("dung", "soil"), c("gut", "soil"), c("gut", "dung"))) {
      HM <- hh_matrix(S, kk, pr[1], pr[2]); if (is.null(HM)) next
      pt <- perm_test(HM)
      res_share[[length(res_share) + 1]] <- data.frame(
        identity_pct = idpct, layer = layer, pair = paste(pr, collapse = "-"),
        n_features = ncol(Mk), n_households = pt$n_households,
        within = pt$observed, null_mean = pt$null_mean, z = pt$z, ratio = pt$ratio,
        mde_pct = pt$mde_pct, p = pt$p, stringsAsFactors = FALSE)
    }
  }
}

t16a <- do.call(rbind, res_merge); save_tab(t16a, "T16a_threshold_merge_summary.tsv")
t16b <- do.call(rbind, res_cont)
t16b[] <- lapply(t16b, function(x) if (is.numeric(x)) round(x, 4) else x)
save_tab(t16b, "T16b_containment_vs_threshold.tsv")
t16c <- do.call(rbind, res_share)
t16c$q_BH <- p.adjust(t16c$p, method = "BH")
t16c[] <- lapply(t16c, function(x) if (is.numeric(x)) round(x, 4) else x)
save_tab(t16c, "T16c_sharing_vs_threshold.tsv")

logmsg("")
logmsg("== containment vs clustering threshold (SES; the 95%% row is the control) ==")
for (ly in unique(t16b$layer)) for (d in unique(t16b$direction)) {
  r <- t16b[t16b$layer == ly & t16b$direction == d, ]
  r <- r[order(-r$identity_pct), ]
  logmsg("  %-8s %-16s %s", ly, d,
         paste(sprintf("id%d%%: SES %+.2f", r$identity_pct, r$SES), collapse = " | "))
}
logmsg("")
logmsg("== household dung-soil sharing vs clustering threshold (z) ==")
for (ly in unique(t16c$layer)) for (p in unique(t16c$pair)) {
  r <- t16c[t16c$layer == ly & t16c$pair == p, ]; r <- r[order(-r$identity_pct), ]
  logmsg("  %-8s %-10s %s", ly, p,
         paste(sprintf("id%d%%: z %+.2f (BH-adjusted p %.3f)", r$identity_pct, r$z, r$q_BH), collapse = " | "))
}
logmsg("")

writeLines(unlist(log_lines), file.path(LOG_DIR, "14_allele_threshold_sweep.log"))
cat("\nARG variant threshold sweep complete.\n")
