#!/usr/bin/env Rscript
# Shared utils for the same-household pair analyses (24_pair_classifier.R and
# 25_feast_source_tracking.R).
#   r8_layers()                 feature layers aligned to one sample key
#   r8_pair_index()             every unique pair, its stratum and household label
#   r8_encode()                 symmetric pair encodings (presence co-occurrence; minimum abundance)
#   r8_auc()                    rank-based AUC
#   r8_household_folds()        cross-validation folds assigned by household
#   r8_permute_households()     household labels shuffled within compartment
#   r8_permute_for_stratum()    stratum-aware permutation for cross-compartment strata

suppressPackageStartupMessages({ library(vegan) })

# ---------------------------------------------------------------------------
# Feature layers
# ---------------------------------------------------------------------------
r8_layers <- function(include_ko = FALSE) {
  key    <- load_key()
  gating <- load_gating()
  arg    <- gate_layer(get_layer("arg", key), gating)
  k      <- arg$key
  clin   <- clinical_args()

  gene_gating <- load_gating("control_arg_genes_exclusion_list.tsv", "gene-level")
  gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")),
                    stringsAsFactors = FALSE)
  gab <- gab[gab$sample_id %in% k$sample_id, ]
  gs  <- sort(unique(gab$sample_id)); gg <- sort(unique(gab$gene_id))
  A   <- matrix(0, length(gs), length(gg), dimnames = list(gs, gg))
  A[cbind(gab$sample_id, gab$gene_id)] <- gab$count
  A   <- A[, !colnames(A) %in% gene_gating$features, drop = FALSE]

  hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")),
                     stringsAsFactors = FALSE)
  gene2sym   <- setNames(hits$arg_symbol, hits$gene_id)
  clin_genes <- names(gene2sym)[gene2sym %in% clin]

  # every layer is aligned to the same sample order
  if (!all(k$sample_id %in% rownames(A)))
    stop("variant matrix is missing samples present in the key")
  A <- A[k$sample_id, , drop = FALSE]

  out <- list(
    `ARG variant (clinical)` = A[, colnames(A) %in% clin_genes, drop = FALSE],
    `ARG family (clinical)`  = arg$M[, colnames(arg$M) %in% clin, drop = FALSE]
  )
  if (include_ko) {
    ko <- get_layer("ko_all", key)
    out[["KO (function)"]] <- ko$M[k$sample_id, , drop = FALSE]
  }
  for (nm in names(out)) {
    out[[nm]] <- out[[nm]][, colSums(out[[nm]]) > 0, drop = FALSE]
    stopifnot(identical(rownames(out[[nm]]), k$sample_id))
  }
  list(layers = out, key = k)
}

# ---------------------------------------------------------------------------
# Pair index
# ---------------------------------------------------------------------------
# Returns one row per unique unordered pair: the two sample row indices, the
# compartment-pair stratum, and whether the two share a household.
r8_pair_index <- function(key) {
  n  <- nrow(key)
  ij <- t(combn(n, 2))
  cp <- paste(pmin(as.character(key$compartment[ij[, 1]]),
                   as.character(key$compartment[ij[, 2]])),
              pmax(as.character(key$compartment[ij[, 1]]),
                   as.character(key$compartment[ij[, 2]])), sep = "-")
  data.frame(
    i = ij[, 1], j = ij[, 2],
    stratum = cp,
    same_hh = key$household[ij[, 1]] == key$household[ij[, 2]],
    testable = !cp %in% c("gut-gut", "dung-dung"),
    stringsAsFactors = FALSE
  )
}

R8_STRATUM_ROLE <- c(
  "gut-soil"  = "primary",
  "dung-gut"  = "primary",
  "soil-soil" = "control",
  "dung-soil" = "control"
)
R8_STRATUM_LABEL <- c(
  "gut-soil"  = "gut to soil",
  "dung-gut"  = "gut to dung",
  "soil-soil" = "soil to soil",
  "dung-soil" = "dung to soil"
)

# ---------------------------------------------------------------------------
# Symmetric pair encodings
# ---------------------------------------------------------------------------
# presence  : 1 if the feature is detected in BOTH samples. Unweighted row sum is
#             shared richness. This is where the rare-feature sharing signal lives.
# abundance : elementwise minimum of the two within-sample relative abundances.
#             Unweighted row sum is Bray-Curtis similarity, the statistic measured
#             at AUC 0.500 on the dung-to-soil stratum.
r8_encode <- function(M, pairs, encoding = c("presence", "abundance"),
                      summaries = FALSE) {
  encoding <- match.arg(encoding)
  P <- (M > 0)
  if (encoding == "presence") {
    X <- P[pairs$i, , drop = FALSE] & P[pairs$j, , drop = FALSE]
    storage.mode(X) <- "integer"
  } else {
    R <- tss(M)
    X <- pmin(R[pairs$i, , drop = FALSE], R[pairs$j, , drop = FALSE])
  }
  dimnames(X) <- list(NULL, colnames(M))
  if (!summaries) return(X)

  prev <- colMeans(P)
  S <- lapply(c(1.00, 0.25, 0.10), function(cap) {
    keep <- prev <= cap
    if (sum(keep) < 10) return(NULL)
    Pk <- P[, keep, drop = FALSE]
    inter <- rowSums(Pk[pairs$i, , drop = FALSE] & Pk[pairs$j, , drop = FALSE])
    sz    <- rowSums(Pk)
    union <- sz[pairs$i] + sz[pairs$j] - inter
    m <- cbind(inter, inter / pmax(union, 1))
    colnames(m) <- sprintf(c("shared_count_prev%.2f", "jaccard_prev%.2f"), cap)
    m
  })
  cbind(X, do.call(cbind, S))
}

# ---------------------------------------------------------------------------
# AUC, rank-based
# ---------------------------------------------------------------------------
# Mann-Whitney identity
r8_auc <- function(score, label) {
  np <- sum(label); nn <- sum(!label)
  if (np == 0 || nn == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[label]) - np * (np + 1) / 2) / (np * nn)
}

# ---------------------------------------------------------------------------
# Grouped folds
# ---------------------------------------------------------------------------
r8_household_folds <- function(households, k = 5) {
  hh <- sort(unique(households))
  f  <- sample(rep_len(seq_len(k), length(hh)))
  setNames(f, hh)
}

# ---------------------------------------------------------------------------
# null
# ---------------------------------------------------------------------------
# Shuffle household labels within compartment
r8_permute_households <- function(key) {
  hh <- key$household
  for (cc in levels(key$compartment)) {
    idx <- which(key$compartment == cc)
    hh[idx] <- sample(hh[idx])
  }
  hh
}

# ---------------------------------------------------------------------------
# Stratum-aware permutation use this for cross-compartment strata.
# ---------------------------------------------------------------------------

r8_permute_for_stratum <- function(key, stratum) {
  parts <- strsplit(stratum, "-", fixed = TRUE)[[1]]
  hh <- key$household
  if (parts[1] == parts[2]) {                 # within-compartment: free
    idx <- which(key$compartment == parts[1])
    hh[idx] <- sample(hh[idx])
    return(hh)
  }
  # cross-compartment: permute the side with the fewest samples per household
  per_hh <- vapply(parts, function(cc)
    mean(table(key$household[key$compartment == cc])), numeric(1))
  side <- parts[which.min(per_hh)]
  idx  <- which(key$compartment == side)
  hh[idx] <- sample(hh[idx])
  hh
}
