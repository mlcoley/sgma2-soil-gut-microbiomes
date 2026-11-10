#!/usr/bin/env Rscript
# QC and batch effect investigation
#   Check 1  read depth by sample type: T17a
#   Check 2  separation of controls from biological samples, on Woltka genus: T17b
#   Check 3  sequencing-plate effects, tested overall and within each compartment
#            because plate is unbalanced across compartments: T17c; and whether the
#            dung-soil household signal survives a plate-restricted permutation: T17d
# Usage: Rscript analysis/08_qc_batch_audit.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(vegan) })
source("analysis/00_load_features.R")

set.seed(1)
N_PERM <- 999

log_lines <- c(sprintf("QC and batch audit - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

meta <- read.delim(file.path(QC_DIR, "decontam_metadata.tsv"),
                   stringsAsFactors = FALSE)
prep <- read.delim(file.path(QC_DIR, "646_metag-prep-info.txt"),
                   stringsAsFactors = FALSE, check.names = FALSE)
prep$sample_name <- trimws(prep$sample_name)
meta$plate <- trimws(prep$sample_plate[match(meta$field_name, prep$sample_name)])
meta$well  <- trimws(prep$well_id     [match(meta$field_name, prep$sample_name)])
meta$run   <- trimws(prep$run_id      [match(meta$field_name, prep$sample_name)])
logmsg("libraries with plate assigned: %d of %d | distinct plates: %d | distinct runs: %d",
       sum(!is.na(meta$plate) & meta$plate != ""), nrow(meta),
       length(unique(na.omit(meta$plate))), length(unique(na.omit(meta$run))))

# ===========================================================================
# Check 1: read depth by sample type -> T17a
# ===========================================================================
logmsg("")
logmsg("== check 1: read depth by sample type ==")
t17a <- meta %>% group_by(sample_type) %>%
  summarise(n = n(), min_reads = min(posthost_reads), median_reads = median(posthost_reads),
            max_reads = max(posthost_reads), .groups = "drop") %>% as.data.frame()
for (i in seq_len(nrow(t17a)))
  logmsg("  %-12s n=%3d  min %11.0f  median %11.0f  max %11.0f",
         t17a$sample_type[i], t17a$n[i], t17a$min_reads[i], t17a$median_reads[i], t17a$max_reads[i])

blanks <- meta[meta$sample_type %in% c("neg_blank", "pcr_blank"), ]
bio    <- meta[meta$sample_type == "true_sample", ]
analyzed <- rownames(get_layer("ko_all", load_key())$M)
bio_an <- bio[bio$sample_id %in% analyzed, ]
logmsg("shallowest 'biological' overall is %s at %.0f reads — the failed library already dropped",
       bio$sample_id[which.min(bio$posthost_reads)], min(bio$posthost_reads))
logmsg("SEPARATION (analyzed cohort only): deepest blank %.0f vs shallowest analyzed biological %.0f -> %.1fx gap",
       max(blanks$posthost_reads), min(bio_an$posthost_reads),
       min(bio_an$posthost_reads) / max(blanks$posthost_reads))
logmsg("  median blank %.0f vs median analyzed biological %.0f -> %.0fx",
       median(blanks$posthost_reads), median(bio_an$posthost_reads),
       median(bio_an$posthost_reads) / max(median(blanks$posthost_reads), 1))
flagged <- blanks[blanks$posthost_reads > 1e5, c("sample_id", "field_name", "sample_type", "posthost_reads")]
if (nrow(flagged)) {
  logmsg("blanks above 100k reads (%d):", nrow(flagged))
  for (i in seq_len(nrow(flagged)))
    logmsg("    %s (%s, %s): %.0f reads", flagged$sample_id[i], flagged$field_name[i],
           flagged$sample_type[i], flagged$posthost_reads[i])
} else logmsg("no blank exceeds 100k reads")
t17a$verdict <- ifelse(t17a$sample_type %in% c("neg_blank", "pcr_blank"),
                       "should be near-empty", "should be deep")
save_tab(t17a, "T17a_depth_by_sample_type.tsv")

# ===========================================================================
# Check 2: controls vs biological samples (Woltka genus) -> T17b
# ===========================================================================
logmsg("")
logmsg("== check 2: controls vs biological samples (Woltka genus) ==")
G <- load_matrix(file.path(FT_DIR, "cohort_genus.tsv"), "woltka")
G <- G[rowSums(G) > 0, , drop = FALSE]                      # drop zero-read libraries
mm <- meta[match(rownames(G), meta$sample_id), ]
keep <- !is.na(mm$sample_type)
G <- G[keep, , drop = FALSE]; mm <- mm[keep, ]
mm$class <- ifelse(mm$sample_type == "true_sample", "biological",
            ifelse(mm$sample_type %in% c("neg_blank", "pcr_blank"), "blank",
            ifelse(mm$sample_type == "pos_dup", "positive control", "mock")))
logmsg("libraries in ordination: %d (%s)", nrow(G),
       paste(sprintf("%s %d", names(table(mm$class)), as.integer(table(mm$class))), collapse = ", "))

D <- vegdist(tss(G), "bray")
a1 <- adonis2(D ~ mm$class, permutations = N_PERM, by = "terms")
r2null <- (length(unique(mm$class)) - 1) / (nrow(G) - 1)
logmsg("PERMANOVA composition ~ control class: R2 = %.3f (random-grouping expectation %.3f), p = %.3f",
       a1$R2[1], r2null, a1$`Pr(>F)`[1])

# distance of each library to its nearest biological sample
Dm <- as.matrix(D)
bio_i <- which(mm$class == "biological")
nn_to_bio <- vapply(seq_len(nrow(G)), function(i)
  min(Dm[i, setdiff(bio_i, i)]), numeric(1))
t17b <- data.frame(sample_id = rownames(G), field_name = mm$field_name, class = mm$class,
                   sample_type = mm$sample_type, plate = mm$plate,
                   posthost_reads = mm$posthost_reads,
                   nearest_biological_BC = round(nn_to_bio, 4), stringsAsFactors = FALSE)
save_tab(t17b, "T17b_control_separation.tsv")
for (cl in c("biological", "blank", "positive control", "mock"))
  if (any(mm$class == cl))
    logmsg("  median distance to nearest biological — %-17s %.3f", cl,
           median(nn_to_bio[mm$class == cl]))

# positive controls are aliquots of donor study samples
pos <- mm[mm$class == "positive control" & !is.na(mm$pos_dup_of), ]
if (nrow(pos)) {
  hit <- 0
  for (i in seq_len(nrow(pos))) {
    s <- pos$sample_id[i]; don <- pos$pos_dup_of[i]
    if (!don %in% rownames(Dm)) next
    nearest <- names(sort(Dm[s, setdiff(rownames(Dm), s)]))[1]
    same_comp <- mm$compartment[match(nearest, mm$sample_id)] == mm$compartment[match(s, mm$sample_id)]
    if (isTRUE(same_comp)) hit <- hit + 1
    logmsg("    %s (donor %s): nearest library = %s, rank of donor = %d",
           s, don, nearest, which(names(sort(Dm[s, setdiff(rownames(Dm), s)])) == don)[1])
  }
  logmsg("  positive controls whose nearest neighbour is the same compartment: %d/%d", hit, nrow(pos))
}

# ===========================================================================
# Check 3: plate / batch effects -> T17c, T17d
# ===========================================================================
logmsg("")
logmsg("== check 3: plate / batch ==")
key <- load_key()
key$plate <- meta$plate[match(key$sample_id, meta$sample_id)]
tab <- table(key$plate, key$compartment)
logmsg("plate x compartment (analyzed cohort):")
for (p in rownames(tab))
  logmsg("  %-10s gut=%2d dung=%2d soil=%3d", p, tab[p, "gut"], tab[p, "dung"], tab[p, "soil"])
logmsg("plate is unbalanced across compartments, so plate is also tested within each compartment")

LAY <- list("KO (function)" = get_layer("ko_all", key)$M,
            "ARG (filtered)" = gate_layer(get_layer("arg", key), load_gating())$M,
            "genus (taxonomy)" = G[rownames(G) %in% key$sample_id, , drop = FALSE])
t17c <- list()
for (ln in names(LAY)) {
  M <- LAY[[ln]]; kk <- key[match(rownames(M), key$sample_id), ]
  ok <- !is.na(kk$plate) & kk$plate != "" & rowSums(M) > 0
  M <- M[ok, , drop = FALSE]; kk <- kk[ok, ]
  d <- vegdist(tss(M), "bray")
  a <- adonis2(d ~ kk$plate, permutations = N_PERM, by = "terms")
  g <- length(unique(kk$plate))
  t17c[[length(t17c) + 1]] <- data.frame(
    layer = ln, scope = "all compartments", n = nrow(M), groups = g,
    R2 = round(a$R2[1], 4), R2_null_expected = round((g - 1) / (nrow(M) - 1), 4),
    R2_excess = round(a$R2[1] - (g - 1) / (nrow(M) - 1), 4), p = a$`Pr(>F)`[1])
  for (cc in COMP_LEVELS) {
    i <- kk$compartment == cc
    if (sum(i) < 15 || length(unique(kk$plate[i])) < 2) next
    dd <- vegdist(tss(M[i, , drop = FALSE]), "bray")
    aa <- adonis2(dd ~ kk$plate[i], permutations = N_PERM, by = "terms")
    gg <- length(unique(kk$plate[i]))
    t17c[[length(t17c) + 1]] <- data.frame(
      layer = ln, scope = paste("within", cc), n = sum(i), groups = gg,
      R2 = round(aa$R2[1], 4), R2_null_expected = round((gg - 1) / (sum(i) - 1), 4),
      R2_excess = round(aa$R2[1] - (gg - 1) / (sum(i) - 1), 4), p = aa$`Pr(>F)`[1])
  }
}
t17c <- do.call(rbind, t17c); save_tab(t17c, "T17c_plate_permanova.tsv")
for (i in seq_len(nrow(t17c)))
  logmsg("  %-17s %-18s n=%3d R2=%.4f (null %.4f, excess %+.4f) p=%.3f",
         t17c$layer[i], t17c$scope[i], t17c$n[i], t17c$R2[i],
         t17c$R2_null_expected[i], t17c$R2_excess[i], t17c$p[i])

## ---- is the dung-soil household signal a plate effect -------------------
logmsg("")
logmsg("== check 3b: dung-soil household signal vs plate ==")
GENE_GATING <- load_gating("control_arg_genes_exclusion_list.tsv", "gene-level")
gab <- read.delim(gzfile(file.path(FT_DIR, "arg_gene_abundance.tsv.gz")), stringsAsFactors = FALSE)
gab <- gab[gab$sample_id %in% key$sample_id & !gab$gene_id %in% GENE_GATING$features, ]
hits <- read.delim(gzfile(file.path(FT_DIR, "arg_hits_long.tsv.gz")), stringsAsFactors = FALSE)
g2s <- setNames(hits$arg_symbol, hits$gene_id)
clin_g <- names(g2s)[g2s %in% clinical_args()]
gab <- gab[gab$gene_id %in% clin_g, ]
ss <- sort(unique(gab$sample_id)); gg <- sort(unique(gab$gene_id))
A <- matrix(0, length(ss), length(gg), dimnames = list(ss, gg))
A[cbind(gab$sample_id, gab$gene_id)] <- gab$count
P <- A >= 1
kk <- key[match(rownames(P), key$sample_id), ]

inter <- tcrossprod(P * 1); sz <- rowSums(P)
S <- inter / pmax(outer(sz, sz, "+") - inter, 1); dimnames(S) <- list(rownames(P), rownames(P))

di <- which(kk$compartment == "dung"); si <- which(kk$compartment == "soil")
hh <- sort(intersect(unique(kk$household[di]), unique(kk$household[si])))
HM <- matrix(NA_real_, length(hh), length(hh), dimnames = list(hh, hh))
for (a in seq_along(hh)) {
  ra <- di[kk$household[di] == hh[a]]
  for (b in seq_along(hh)) HM[a, b] <- mean(S[ra, si[kk$household[si] == hh[b]], drop = FALSE])
}
dung_plate <- kk$plate[di][match(hh, kk$household[di])]

# plate sharing between own- and between-household dung-soil pairs
own <- 0; own_n <- 0; btw <- 0; btw_n <- 0
for (a in seq_along(hh)) {
  dp <- dung_plate[a]
  for (b in seq_along(hh)) {
    sp <- kk$plate[si][kk$household[si] == hh[b]]
    if (a == b) { own <- own + sum(sp == dp); own_n <- own_n + length(sp) }
    else        { btw <- btw + sum(sp == dp); btw_n <- btw_n + length(sp) }
  }
}
logmsg("  dung-soil pairs sharing a plate: own-household %d/%d = %.1f%% | between %d/%d = %.1f%% (%.2fx)",
       own, own_n, 100*own/own_n, btw, btw_n, 100*btw/btw_n, (own/own_n)/(btw/btw_n))

perm_z <- function(M, strata = NULL, n_perm = N_PERM) {
  n <- nrow(M); obs <- mean(diag(M))
  nul <- replicate(n_perm, {
    idx <- seq_len(n)
    if (is.null(strata)) idx <- sample(idx)
    else for (g in unique(strata)) { j <- which(strata == g); idx[j] <- j[sample(length(j))] }
    mean(M[cbind(seq_len(n), idx)])
  })
  c(obs = obs, null = mean(nul), z = (obs - mean(nul)) / sd(nul),
    p = (1 + sum(nul >= obs)) / (n_perm + 1))
}
free <- perm_z(HM)
rest <- perm_z(HM, strata = dung_plate)
logmsg("  UNRESTRICTED permutation : z = %+.2f, p = %.3f", free["z"], free["p"])
logmsg("  PLATE-RESTRICTED         : z = %+.2f, p = %.3f  (labels shuffled only within a dung plate)",
       rest["z"], rest["p"])
t17d <- data.frame(test = c("unrestricted", "plate-restricted"),
                   observed = c(free["obs"], rest["obs"]), null_mean = c(free["null"], rest["null"]),
                   z = c(free["z"], rest["z"]), p = c(free["p"], rest["p"]),
                   own_plate_share = round(100*own/own_n, 1),
                   between_plate_share = round(100*btw/btw_n, 1))
save_tab(t17d, "T17d_dungsoil_plate_restricted.tsv")

writeLines(unlist(log_lines), file.path(LOG_DIR, "08_qc_batch_audit.log"))
cat("\nQC and batch audit complete.\n")
