#!/usr/bin/env Rscript
# ARG-variant containment: sensitivity to which sample supplied each catalog
# representative
# Output: T10.
# Usage: Rscript analysis/13_allele_source_sensitivity.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/00_load_features.R")

set.seed(1)
B_SUB      <- 100
THRESHOLDS <- c(1, 5)

log_lines <- c(sprintf("ARG variant source sensitivity - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

key  <- load_key()
ann  <- load_arg_annotation()
clin <- clinical_args(ann)

# --- source compartment of every catalog gene ------------------------------
# The gene id prefix is the sample its representative was assembled from
full  <- load_key(biological_only = FALSE, soil_primary_only = FALSE)
scomp <- setNames(as.character(full$compartment), full$sample_id)

GENE_GATING <- load_gating("control_arg_genes_exclusion_list.tsv", "gene-level")
gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")), stringsAsFactors = FALSE)
gab <- gab[gab$sample_id %in% key$sample_id & !gab$gene_id %in% GENE_GATING$features, ]

genes  <- sort(unique(gab$gene_id))
g_src  <- scomp[sub("_k141.*$", "", genes)]
stopifnot(!any(is.na(g_src)))
names(g_src) <- genes
logmsg("catalog genes in cohort: %d | source composition: %s", length(genes),
       paste(sprintf("%s %.1f%%", names(table(g_src)), 100 * prop.table(table(g_src))),
             collapse = " "))

samps <- sort(unique(gab$sample_id))
M <- matrix(0L, nrow = length(samps), ncol = length(genes), dimnames = list(samps, genes))
M[cbind(gab$sample_id, gab$gene_id)] <- gab$count
kk <- key[match(rownames(M), key$sample_id), ]

gene2sym <- with(read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")),
                            stringsAsFactors = FALSE), setNames(arg_symbol, gene_id))
clin_genes <- names(gene2sym)[gene2sym %in% clin]

# --- containment with the exact hypergeometric null ------------------------
# |A n B| ~ Hypergeom(N = |union|, K = |B|, n = |A|), so E[C(A->B)] = K/N and
# sd[C] = sqrt(n * mu * (1 - mu) * (N - n) / (N - 1)) / n.
contain_ses <- function(sets, A, B) {
  N <- length(unique(unlist(sets))); K <- length(B); n <- length(A)
  obs <- if (n) length(intersect(A, B)) / n else NA_real_
  mu  <- K / N
  sdC <- if (n > 0 && N > 1) sqrt(n * mu * (1 - mu) * (N - n) / (N - 1)) / n else NA_real_
  c(observed = obs, null = mu,
    SES = if (is.na(sdC) || sdC == 0) NA_real_ else (obs - mu) / sdC)
}

run_one <- function(cols, layer_label, subset_label, thr) {
  Ms <- M[, cols, drop = FALSE]
  if (ncol(Ms) < 20) return(NULL)
  n_eq <- min(table(kk$compartment)[COMP_LEVELS])
  reps <- replicate(B_SUB, {
    idx <- lapply(setNames(COMP_LEVELS, COMP_LEVELS),
                  function(cc) sample(which(kk$compartment == cc), n_eq))
    S <- lapply(idx, function(i) which(colSums(Ms[i, , drop = FALSE] >= thr) > 0))
    S$env <- union(S$dung, S$soil)
    vapply(list(`gut -> dung`  = c("gut", "dung"),
                `dung -> soil` = c("dung", "soil"),
                `gut -> soil`  = c("gut", "soil"),
                `gut -> env`   = c("gut", "env")),
           function(p) contain_ses(S[c("gut","dung","soil")], S[[p[1]]], S[[p[2]]]),
           numeric(3))
  }, simplify = "array")
  dirs <- dimnames(reps)[[2]]
  data.frame(
    layer = layer_label, gene_subset = subset_label, n_features = ncol(Ms),
    threshold = thr, direction = dirs,
    observed = round(apply(reps["observed", , ], 1, mean), 4),
    null     = round(apply(reps["null", , ], 1, mean), 4),
    SES      = round(apply(reps["SES", , ], 1, mean), 3),
    SES_sd   = round(apply(reps["SES", , ], 1, sd), 3),
    row.names = NULL)
}

SUBSETS <- list(
  "all sources"        = names(g_src),
  "non-gut-sourced"    = names(g_src)[g_src != "gut"],
  "gut-sourced only"   = names(g_src)[g_src == "gut"]
)

res <- list()
for (lab in c("ARG allele (clinical)", "ARG allele (all)")) {
  base <- if (lab == "ARG allele (clinical)") intersect(colnames(M), clin_genes) else colnames(M)
  for (sn in names(SUBSETS)) for (thr in THRESHOLDS)
    res[[length(res)+1]] <- run_one(intersect(base, SUBSETS[[sn]]), lab, sn, thr)
}
t10 <- do.call(rbind, res)

# Circular cells
t10$evaluable <- with(t10, ifelse(
  gene_subset == "non-gut-sourced" & direction == "gut -> env", "NO — circular (env holds nearly all non-gut-sourced genes)",
  ifelse(gene_subset == "gut-sourced only" & grepl("^gut ->", direction),
         "NO — circular (gut holds essentially all gut-sourced genes)", "yes")))
t10$contamination_filter <- gating_col(GENE_GATING)
save_tab(t10, "T10_allele_source_sensitivity.tsv")
logmsg("informative rows: %d of %d (rest circular by construction; see `evaluable`)",
       sum(t10$evaluable == "yes"), nrow(t10))

for (lab in unique(t10$layer)) {
  logmsg("== %s ==", lab)
  for (d in c("gut -> dung", "dung -> soil", "gut -> env")) {
    r <- t10[t10$layer == lab & t10$direction == d & t10$threshold == 1, ]
    logmsg("  %-13s %s", d, paste(sprintf("[%s, n=%d] obs %.3f null %.3f SES %+.2f",
           r$gene_subset, r$n_features, r$observed, r$null, r$SES), collapse = " | "))
  }
}

# gut -> dung SES with and without gut-sourced representatives
gd <- t10[t10$direction == "gut -> dung" & t10$threshold == 1, ]
for (lab in unique(gd$layer)) {
  a <- gd$SES[gd$layer == lab & gd$gene_subset == "all sources"]
  b <- gd$SES[gd$layer == lab & gd$gene_subset == "non-gut-sourced"]
  logmsg("[%s] gut->dung SES all sources %+.2f -> non-gut-sourced %+.2f", lab, a, b)
}

writeLines(log_lines, file.path(LOG_DIR, "13_allele_source_sensitivity.log"))
cat("\nARG variant source sensitivity complete.\n")
