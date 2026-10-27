#!/usr/bin/env Rscript
# Sensitivity of the headline results to the definition of the clinical ARG subset.
# clinical_args() is AMRFinderPlus element_type == "AMR" & scope == "core". Four
# definitions are run; virulence elements are never included.
#   A  AMR & core                 (used in the paper)
#   B  AMR, any scope             (adds 9 mostly chromosomal efflux/related symbols)
#   C  AMR & core + STRESS & core (adds 3 quaternary-ammonium resistance symbols)
#   D  AMR any scope + STRESS core
# n_symbols counts the definition's symbols that are columns of the wide ARG table
# (all 229 libraries); n_symbols_detected counts those detected in the 194 analysed
# samples (170 under A). One construction is held fixed and only the subset varies, so
# columns are comparable to each other; there is no equal-effort subsampling, so SES
# and z are not comparable to T8b. Read the spread across rows, not the level.
# Output: T20.
# Usage: Rscript analysis/10_clinical_subset.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
source("analysis/00_load_features.R")

set.seed(1)
N_PERM <- 9999

log_lines <- c(sprintf("Clinical-subset sensitivity - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

ann  <- load_arg_annotation()
key  <- load_key()
GAT  <- load_gating()
GGAT <- load_gating("control_arg_genes_blocklist.tsv", "gene-level")

DEFS <- list(
  A = list(label = "AMR & core (current)",        sel = ann$element_type == "AMR" & ann$scope == "core"),
  B = list(label = "AMR, any scope",              sel = ann$element_type == "AMR"),
  C = list(label = "AMR core + STRESS core",      sel = (ann$element_type == "AMR" & ann$scope == "core") |
                                                        (ann$element_type == "STRESS" & ann$scope == "core")),
  D = list(label = "AMR any + STRESS core",       sel = ann$element_type == "AMR" |
                                                        (ann$element_type == "STRESS" & ann$scope == "core"))
)
for (d in names(DEFS)) DEFS[[d]]$symbols <- ann$arg_symbol[DEFS[[d]]$sel]

## ---- symbol layer ---------------------------------------------------------
Marg <- load_matrix(file.path(FT_DIR, LAYERS$arg$file))
Marg <- Marg[rownames(Marg) %in% key$sample_id, , drop = FALSE]
kA   <- key[match(rownames(Marg), key$sample_id), ]
gf   <- genome_factor(kA)$factor

## ---- ARG variant layer ----------------------------------------------------
gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")), stringsAsFactors = FALSE)
gab <- gab[gab$sample_id %in% key$sample_id & !gab$gene_id %in% GGAT$features, ]
hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")), stringsAsFactors = FALSE)
g2s  <- setNames(hits$arg_symbol, hits$gene_id)

## ---- shared machinery -----------------------------------------------------
ses <- function(S, A, B) {
  N <- length(unique(unlist(S))); K <- length(B); n <- length(A)
  if (!n || !N) return(NA_real_)
  mu <- K / N
  sdv <- sqrt(n * mu * (1 - mu) * (N - n) / (N - 1)) / n
  if (is.na(sdv) || sdv == 0) NA_real_ else (length(intersect(A, B)) / n - mu) / sdv
}
comp_sets <- function(P, kk) lapply(setNames(COMP_LEVELS, COMP_LEVELS),
  function(c) { j <- which(kk$compartment == c)
                if (!length(j)) character(0) else colnames(P)[colSums(P[j, , drop = FALSE]) > 0] })
perm_z <- function(M, n_perm = N_PERM) {
  n <- nrow(M); obs <- mean(diag(M))
  nul <- replicate(n_perm, mean(M[cbind(seq_len(n), sample.int(n))]))
  c(z = (obs - mean(nul)) / sd(nul), p = (1 + sum(nul >= obs)) / (n_perm + 1))
}

out <- list()
for (d in names(DEFS)) {
  sy  <- setdiff(DEFS[[d]]$symbols, GAT$features)          # symbol-level exclusion list
  cols <- intersect(sy, colnames(Marg))
  Ms  <- Marg[, cols, drop = FALSE]
  Ps  <- Ms > MIN_COUNT

  # richness and per-genome abundance by compartment
  rich <- rowSums(Ps)
  cpg  <- rowSums(Ms) / gf[rownames(Ms)]
  med  <- function(v, cp) median(v[kA$compartment == cp], na.rm = TRUE)

  # containment, ARG family level
  Ssym <- comp_sets(Ps, kA)
  env  <- union(Ssym$dung, Ssym$soil)
  ses_sym <- ses(Ssym, Ssym$gut, env)

  # ARG variant layer under this definition. Subset membership is a property of the
  # definition; contamination exclusion at this layer is the gene-level list's job
  # alone (applied to `gab` above), never the symbol-level list.
  cg  <- names(g2s)[g2s %in% DEFS[[d]]$symbols]
  ga  <- gab[gab$gene_id %in% cg, ]
  ssm <- sort(unique(ga$sample_id)); ggn <- sort(unique(ga$gene_id))
  A <- matrix(0, length(ssm), length(ggn), dimnames = list(ssm, ggn))
  A[cbind(ga$sample_id, ga$gene_id)] <- ga$count
  Pa <- A >= 1
  ka <- key[match(rownames(Pa), key$sample_id), ]
  Sal <- comp_sets(Pa, ka)
  ses_all <- ses(Sal, Sal$gut, union(Sal$dung, Sal$soil))

  # dung-soil household signal, ARG variant level
  it <- tcrossprod(Pa * 1); szv <- rowSums(Pa)
  SJ <- it / pmax(outer(szv, szv, "+") - it, 1); dimnames(SJ) <- list(rownames(Pa), rownames(Pa))
  di <- which(ka$compartment == "dung"); si <- which(ka$compartment == "soil")
  hh <- sort(intersect(unique(ka$household[di]), unique(ka$household[si])))
  HM <- matrix(NA_real_, length(hh), length(hh), dimnames = list(hh, hh))
  for (a in seq_along(hh)) { ra <- di[ka$household[di] == hh[a]]
    for (b in seq_along(hh)) HM[a, b] <- mean(SJ[ra, si[ka$household[si] == hh[b]], drop = FALSE]) }
  pz <- perm_z(HM)

  out[[d]] <- data.frame(
    def = d, definition = DEFS[[d]]$label,
    n_symbols = length(cols), n_symbols_detected = sum(colSums(Ps) > 0),
    n_alleles = ncol(Pa),
    rich_gut = med(rich, "gut"), rich_dung = med(rich, "dung"), rich_soil = med(rich, "soil"),
    cpg_gut = med(cpg, "gut"), cpg_dung = med(cpg, "dung"), cpg_soil = med(cpg, "soil"),
    gut_over_soil = med(cpg, "gut") / med(cpg, "soil"),
    SES_gut_env_symbol = ses_sym, SES_gut_env_allele = ses_all,
    dungsoil_z = unname(pz["z"]), dungsoil_p = unname(pz["p"]))
}
t20 <- do.call(rbind, out)
save_tab(t20, "T20_clinical_subset_sensitivity.tsv")

logmsg("Clinical-subset definitions and what each changes:")
logmsg("")
logmsg("  %-3s %-28s %8s %9s %9s", "def", "definition", "symbols", "detected", "alleles")
for (i in seq_len(nrow(t20)))
  logmsg("  %-3s %-28s %8d %9d %9d", t20$def[i], t20$definition[i], t20$n_symbols[i],
         t20$n_symbols_detected[i], t20$n_alleles[i])
logmsg("")
logmsg("  %-3s %10s %10s %10s | %14s", "def", "cpg gut", "cpg soil", "gut/soil", "richness order")
for (i in seq_len(nrow(t20))) {
  r <- c(gut = t20$rich_gut[i], dung = t20$rich_dung[i], soil = t20$rich_soil[i])
  logmsg("  %-3s %10.3f %10.3f %10.2fx | %s", t20$def[i], t20$cpg_gut[i], t20$cpg_soil[i],
         t20$gut_over_soil[i], paste(names(sort(r, decreasing = TRUE)), collapse = " > "))
}
logmsg("")
logmsg("  %-3s %18s %18s %12s", "def", "SES gut->env sym", "SES gut->env all", "dung-soil z")
for (i in seq_len(nrow(t20)))
  logmsg("  %-3s %18.2f %18.2f %12s", t20$def[i], t20$SES_gut_env_symbol[i],
         t20$SES_gut_env_allele[i], sprintf("%+.2f (p %.4f)", t20$dungsoil_z[i], t20$dungsoil_p[i]))

logmsg("")
logmsg("SES and z above are not effort-controlled; read the spread across definitions, not the level.")
logmsg("")
logmsg("Checks across definitions:")
chk <- function(name, v, tol_rel = NULL) {
  rng <- range(v); same_sign <- all(sign(v) == sign(v[1]))
  msg <- if (!is.null(tol_rel)) {
    sprintf("%s: %.2f to %.2f (%.1f%% spread)", name, rng[1], rng[2],
            100 * (rng[2] - rng[1]) / abs(v[1]))
  } else sprintf("%s: %.2f to %.2f", name, rng[1], rng[2])
  logmsg("  %-46s %s", msg, if (same_sign) "sign STABLE" else "SIGN FLIPS")
}
logmsg("  richness ordering identical across definitions: %s",
       if (length(unique(apply(t20[, c("rich_gut","rich_dung","rich_soil")], 1,
           function(r) paste(order(r, decreasing = TRUE), collapse = "")))) == 1) "YES" else "NO")
chk("gut/soil per-genome ratio", t20$gut_over_soil, TRUE)
chk("SES gut->env (symbol)", t20$SES_gut_env_symbol)
chk("SES gut->env (allele)", t20$SES_gut_env_allele)
chk("dung-soil household z", t20$dungsoil_z)
logmsg("  dung-soil significant (p < 0.05) under all definitions: %s",
       if (all(t20$dungsoil_p < 0.05)) "YES" else "NO")

writeLines(unlist(log_lines), file.path(LOG_DIR, "10_clinical_subset.log"))
cat("\nClinical-subset sensitivity complete.\n")
