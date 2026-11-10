#!/usr/bin/env Rscript
# Sensitivity of presence definitions: read count >= 1 versus covered-fraction (breadth) >= 0.5.
# Inputs (built by pipeline/11_feature-tables/rebuild_breadth_filtered.sbatch):
#   feature_tables/ko_abundance_wide_cf50.tsv.gz
#   feature_tables/arg_abundance_wide_cf50.tsv.gz
#   feature_tables/arg_gene_abundance.tsv.gz (with a covered_fraction column)
# Outputs: T18a (median richness by rule), T18b (containment SES and dung-soil
# household-sharing z by rule), T18c (T18b plus containment SES for child stool in
# cattle dung and child stool in soil). T18c is computed inside the same equal-effort
# draws as T18b, so no extra random draw is made and every T18b column is unchanged.
# Usage: Rscript analysis/09_breadth_sensitivity.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(vegan) })
source("analysis/00_load_features.R")

set.seed(1)
N_PERM <- 999; B_SUB <- 100; CF <- 0.5

log_lines <- c(sprintf("Presence-rule sensitivity - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

ko_cf  <- file.path(FT_DIR, "ko_abundance_wide_cf50.tsv.gz")
arg_cf <- file.path(FT_DIR, "arg_abundance_wide_cf50.tsv.gz")
if (!file.exists(ko_cf) || !file.exists(arg_cf)) {
  stop("breadth-filtered tables not present in ", FT_DIR,
       "; run pipeline/11_feature-tables/rebuild_breadth_filtered.sbatch")
}

key   <- load_key()
clin  <- clinical_args()
GAT   <- load_gating()

## ---- two feature spaces ----------------------------------------------
align <- function(M) { i <- intersect(rownames(M), key$sample_id); M[i, colSums(M[i, , drop = FALSE]) > 0, drop = FALSE] }
KO_a  <- align(get_layer("ko_all", key)$M)
KO_b  <- align(load_matrix(ko_cf,  "woltka"))
AR_a  <- align(drop_gated(get_layer("arg", key)$M,  GAT))
AR_b  <- align(drop_gated(load_matrix(arg_cf, "woltka"), GAT))
logmsg("feature counts — KO: %d (count>=1) vs %d (breadth>=%.1f) | ARG: %d vs %d",
       ncol(KO_a), ncol(KO_b), CF, ncol(AR_a), ncol(AR_b))

# ARG variant layer
gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")), stringsAsFactors = FALSE)
has_cf <- "covered_fraction" %in% names(gab)
logmsg("allele table carries covered_fraction: %s", has_cf)
GG <- load_gating("control_arg_genes_exclusion_list.tsv", "gene-level")
gab <- gab[gab$sample_id %in% key$sample_id & !gab$gene_id %in% GG$features, ]
hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")), stringsAsFactors = FALSE)
g2s  <- setNames(hits$arg_symbol, hits$gene_id)
mk_allele <- function(d) {
  d <- d[d$gene_id %in% names(g2s)[g2s %in% clin], ]
  ss <- sort(unique(d$sample_id)); gg <- sort(unique(d$gene_id))
  M <- matrix(0, length(ss), length(gg), dimnames = list(ss, gg))
  M[cbind(d$sample_id, d$gene_id)] <- d$count; M
}
AL_a <- mk_allele(gab)
AL_b <- if (has_cf) mk_allele(gab[suppressWarnings(as.numeric(gab$covered_fraction)) >= CF, ]) else NULL

SPACES <- list("KO (function)" = list(a = KO_a, b = KO_b),
               "ARG (clinical)" = list(a = AR_a[, colnames(AR_a) %in% clin, drop = FALSE],
                                       b = AR_b[, colnames(AR_b) %in% clin, drop = FALSE]))
if (!is.null(AL_b)) SPACES[["ARG allele (clinical)"]] <- list(a = AL_a, b = AL_b)

## ---- richness under each presence def ----------------------------------
logmsg("")
logmsg("== richness under each presence rule ==")
t18a <- list()
for (ln in names(SPACES)) {
  for (v in c("a", "b")) {
    M <- SPACES[[ln]][[v]]; kk <- key[match(rownames(M), key$sample_id), ]
    r <- tapply(rowSums(M > 0), kk$compartment, median)
    t18a[[length(t18a) + 1]] <- data.frame(
      layer = ln, rule = if (v == "a") "count >= 1" else sprintf("breadth >= %.1f", CF),
      gut = r[["gut"]], dung = r[["dung"]], soil = r[["soil"]])
  }
}
t18a <- do.call(rbind, t18a)
save_tab(t18a, "T18a_richness_by_presence_rule.tsv")
for (ln in unique(t18a$layer)) {
  x <- t18a[t18a$layer == ln, ]
  logmsg("  %-21s %-16s gut %7.0f | dung %7.0f | soil %7.0f", ln, x$rule[1], x$gut[1], x$dung[1], x$soil[1])
  logmsg("  %-21s %-16s gut %7.0f | dung %7.0f | soil %7.0f", "", x$rule[2], x$gut[2], x$dung[2], x$soil[2])
  logmsg("  %-21s %-16s gut %+6.1f%% | dung %+6.1f%% | soil %+6.1f%%", "", "-> change",
         100*(x$gut[2]/x$gut[1]-1), 100*(x$dung[2]/x$dung[1]-1), 100*(x$soil[2]/x$soil[1]-1))
}
# Is soil still the richest compartment for KO under the breadth def
ko <- t18a[t18a$layer == "KO (function)", ]
ko_drop <- setNames(1 - as.numeric(ko[2, COMP_LEVELS]) / as.numeric(ko[1, COMP_LEVELS]), COMP_LEVELS)
logmsg("")
logmsg("KO richness loss under the breadth rule: %s",
       paste(sprintf("%s %.1f%%", COMP_LEVELS, 100*ko_drop), collapse = " | "))
logmsg("soil vs gut KO richness: count>=1 %.0f vs %.0f (soil %s); breadth %.0f vs %.0f (soil %s)",
       ko$soil[1], ko$gut[1], if (ko$soil[1] > ko$gut[1]) "higher" else "lower",
       ko$soil[2], ko$gut[2], if (ko$soil[2] > ko$gut[2]) "higher" else "lower")

## ---- containment and household sharing ------------------------------------------
draw_equal <- function(kk, n) unlist(lapply(COMP_LEVELS, function(c){ i <- which(kk$compartment==c); sample(i, min(n, length(i))) }))
csets <- function(M, kk, idx) lapply(setNames(COMP_LEVELS, COMP_LEVELS), function(c){
  j <- idx[kk$compartment[idx]==c]; if(!length(j)) character(0) else colnames(M)[colSums(M[j,,drop=FALSE]>0)>0] })
ses <- function(S, A, B) { N <- length(unique(unlist(S))); K <- length(B); n <- length(A)
  if(!n||!N) return(NA_real_); mu <- K/N
  sd <- sqrt(n*mu*(1-mu)*(N-n)/(N-1))/n; if(is.na(sd)||sd==0) NA_real_ else (length(intersect(A,B))/n - mu)/sd }
jac <- function(P){ it <- tcrossprod(P*1); sz <- rowSums(P); J <- it/pmax(outer(sz,sz,"+")-it,1); dimnames(J) <- list(rownames(P),rownames(P)); J }
hhm <- function(S, kk, c1, c2){ i1 <- which(kk$compartment==c1); i2 <- which(kk$compartment==c2)
  h1 <- as.character(kk$household[i1]); h2 <- as.character(kk$household[i2])
  hs <- sort(intersect(unique(h1),unique(h2))); if(length(hs)<5) return(NULL)
  M <- matrix(NA_real_,length(hs),length(hs),dimnames=list(hs,hs))
  for(a in seq_along(hs)){ ra <- i1[h1==hs[a]]; for(b in seq_along(hs)) M[a,b] <- mean(S[ra,i2[h2==hs[b]],drop=FALSE]) }; M }
pz <- function(M){ n <- nrow(M); o <- mean(diag(M)); nl <- replicate(N_PERM, mean(M[cbind(seq_len(n), sample.int(n))]))
  c(z=(o-mean(nl))/sd(nl), p=(1+sum(nl>=o))/(N_PERM+1)) }

logmsg("")
logmsg("== containment and household sharing under each rule ==")
t18b <- list()
for (ln in names(SPACES)) for (v in c("a","b")) {
  M <- SPACES[[ln]][[v]]; kk <- key[match(rownames(M), key$sample_id), ]
  rule <- if (v=="a") "count >= 1" else sprintf("breadth >= %.1f", CF)
  neq <- min(table(kk$compartment))
  cs <- replicate(B_SUB, { S <- csets(M, kk, draw_equal(kk, neq)); env <- union(S$dung,S$soil)
    c(ses(S,S$gut,env), ses(S,S$dung,S$soil), ses(S,S$gut,S$dung), ses(S,S$gut,S$soil)) })
  S <- jac(M > 0)
  sh <- sapply(list(c("dung","soil"),c("gut","soil"),c("gut","dung")), function(pr){
    H <- hhm(S,kk,pr[1],pr[2]); if(is.null(H)) c(z=NA,p=NA) else pz(H) })
  t18b[[length(t18b)+1]] <- data.frame(layer=ln, rule=rule, n_features=ncol(M),
    SES_gut_env=mean(cs[1,],na.rm=TRUE), SES_dung_soil=mean(cs[2,],na.rm=TRUE),
    SES_gut_dung=mean(cs[3,],na.rm=TRUE), SES_gut_soil=mean(cs[4,],na.rm=TRUE),
    z_dung_soil=sh[1,1], p_dung_soil=sh[2,1], z_gut_soil=sh[1,2], z_gut_dung=sh[1,3])
}
t18b <- do.call(rbind, t18b)
t18b[] <- lapply(t18b, function(x) if(is.numeric(x)) round(x,3) else x)
save_tab(t18b[, !names(t18b) %in% c("SES_gut_dung", "SES_gut_soil")],
         "T18b_containment_sharing_by_presence_rule.tsv")
save_tab(t18b, "T18c_containment_sharing_by_presence_rule_pairs.tsv")
for (ln in unique(t18b$layer)) { x <- t18b[t18b$layer==ln,]
  for (r in 1:2)
    logmsg("  %-21s %-16s nfeat %5d | SES dung->soil %+8.2f | gut->dung %+8.2f | gut->soil %+8.2f | dung-soil z %+5.2f (p %.3f)",
           if (r == 1) ln else "", x$rule[r], x$n_features[r], x$SES_dung_soil[r],
           x$SES_gut_dung[r], x$SES_gut_soil[r], x$z_dung_soil[r], x$p_dung_soil[r]) }

writeLines(unlist(log_lines), file.path(LOG_DIR, "09_breadth_sensitivity.log"))
cat("\nPresence-rule sensitivity complete.\n")
