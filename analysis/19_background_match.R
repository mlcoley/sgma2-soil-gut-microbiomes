#!/usr/bin/env Rscript
# Builds the matched background gene sets.
# Inputs: T23a, T23b (18_background_target_profile.R), feature_tables/
# background_gene_summary.tsv.gz and arg_gene_lengths.tsv (pipeline/12_background-genes).
# Outputs: T25 set files (arm1_magnitude, arm1_shape, arm2_family), T25_arm2_family_map,
# and background_gene_ids.txt (the union of selected genes, input to
# pipeline/12_background-genes/extract_background_abundance.sbatch).
# Usage: Rscript analysis/19_background_match.R   (from the repository root)

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
source("analysis/00_load_features.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(191)
B_DRAW    <- 100    # matched sets per arm
N_BIN_MAG <- 8      # quantile bins, occupancy magnitude
N_BIN_ABN <- 8      # quantile bins, log10 abundance
N_BIN_LEN <- 4      # quantile bins, gene length
OG_POOL   <- 5      # arm2: sample the OG from this many closest candidates, so
                    # draws vary instead of collapsing to one deterministic set
PRIMARY   <- "clinical"

BG_FILE  <- file.path(FT_DIR, "background_gene_summary.tsv.gz")
ARG_LEN  <- file.path(FT_DIR, "arg_gene_lengths.tsv")
N_COMP   <- c(gut = 39, dung = 38, soil = 117)       

cat("\n=== build matched background sets ===\n\n")

## ---- inputs --------------------------------------------------------------
for (f in c(BG_FILE, ARG_LEN)) {
  if (!file.exists(f)) stop(sprintf(
    "missing input: %s\n  built by pipeline/12_background-genes/background_gene_summary.sbatch", f))
}

arg <- read.delim(file.path(TAB_DIR, "T23a_arg_target_profile.tsv"),
                  stringsAsFactors = FALSE)
if (PRIMARY == "clinical") arg <- arg[arg$is_clinical, ]
# arg_gene_lengths.tsv is headerless
alen <- read.delim(ARG_LEN, stringsAsFactors = FALSE, header = FALSE)
if (identical(as.character(alen[1, 1]), "gene_id")) alen <- alen[-1, , drop = FALSE]
names(alen)[1:2] <- c("gene_id", "gene_length")
alen$gene_length <- as.integer(alen$gene_length)
arg <- merge(arg, alen[, c("gene_id", "gene_length")], by = "gene_id")
cat(sprintf("ARG target genes (%s): %d  (lengths joined for %d)\n",
            PRIMARY, nrow(arg), sum(!is.na(arg$gene_length))))
arg <- arg[!is.na(arg$gene_length) & arg$total_count > 0, ]

# decompressed through a pipe (fread gz support needs R.utils)
bg <- fread(cmd = paste("gzip -dc", shQuote(BG_FILE)), sep = "\t",
            header = TRUE, showProgress = FALSE)
cat(sprintf("background pool as delivered: %s genes\n", format(nrow(bg), big.mark = ",")))

## ---- semantic exclusions ---------------------------------------------------
# AMRFinderPlus genes and unannotated genes were already dropped by HPC job
# remaining exclusions, matched on KO identifiers
MOBILITY_KO <- c(
  # transposase / integrase / recombinase / IS elements
  "K07497", "K07496", "K07495", "K07483", "K07484", "K07485", "K07486", "K07487",
  "K07488", "K07489", "K07491", "K07492", "K07493", "K07494", "K14059", "K06400",
  # plasmid conjugation / relaxase / T4SS
  "K03194", "K03195", "K03196", "K03197", "K03199", "K03200", "K03201", "K03202",
  "K03203", "K03204", "K03205", "K12056",
  # phage integrase / structural
  "K01356", "K06400", "K07729"
)
DEFENCE_KO <- c(
  # CRISPR-Cas
  "K09951", "K09952", "K19046", "K19075", "K19076", "K19077", "K19123", "K19124",
  "K19125", "K19126", "K07012", "K09002",
  # restriction-modification
  "K01153", "K01154", "K01155", "K03427", "K06223", "K07316", "K07317", "K07318",
  "K13581"
)
METAL_KO <- c(
  # mer / ars / czc / cop / sil / qac-adjacent efflux
  "K00520", "K08365", "K08366", "K19064", "K03325", "K03893", "K01551", "K03741",
  "K07235", "K07236", "K07237", "K07238", "K07239", "K07240", "K07787", "K07788",
  "K07789", "K17686", "K01533", "K17686"
)
EXCLUDE_KO <- unique(c(MOBILITY_KO, DEFENCE_KO, METAL_KO))

# restricting pool to KO-annotated genes is analogous to the clinical ARG subset 
# limited to well-characterised determinants (element_type == "AMR" & scope == "core").
REQUIRE_KO <- TRUE

n0 <- nrow(bg)
if (REQUIRE_KO) {
  bg <- bg[ko != "-"]
  cat(sprintf("REQUIRE_KO: pool restricted to KO-annotated genes, %s dropped (%.1f%%)\n",
              format(n0 - nrow(bg), big.mark = ","), 100 * (n0 - nrow(bg)) / n0))
}
n1 <- nrow(bg)                     # after REQUIRE_KO, before the semantic screen
bg <- bg[!(ko %in% EXCLUDE_KO)]
cat(sprintf("semantic KO exclusions (mobility/defence/metal): %s genes dropped (%.2f%% of the pre-screen pool)\n",
            format(n1 - nrow(bg), big.mark = ","), 100 * (n1 - nrow(bg)) / n1))
if (!REQUIRE_KO)
  cat(sprintf("  %.1f%% of the pool has no KO ('-') and is not screened\n",
              100 * mean(bg$ko == "-")))

bg <- bg[total_count > 0 & !is.na(gene_length) & gene_length > 0]
cat(sprintf("background pool after exclusions: %s genes\n", format(nrow(bg), big.mark = ",")))

## ---- axes ----------------------------------------------------------------
# magnitude as a rate per compartment so soil does not dominate a raw count
magnitude <- function(g, d, s) (g / N_COMP["gut"] + d / N_COMP["dung"] + s / N_COMP["soil"]) / 3

arg$mag <- magnitude(arg$prev_gut, arg$prev_dung, arg$prev_soil)
arg$abn <- log10(arg$total_count)
arg$len <- arg$gene_length
bg[, mag := magnitude(prev_gut, prev_dung, prev_soil)]
bg[, abn := log10(total_count)]
bg[, len := gene_length]

# shape: the free axis; used for reporting and for shape-matched arm
shape_pattern <- function(g, d, s)
  paste0(ifelse(g > 0, "G", "-"), ifelse(d > 0, "D", "-"), ifelse(s > 0, "S", "-"))
arg$shape <- shape_pattern(arg$prev_gut, arg$prev_dung, arg$prev_soil)
bg[, shape := shape_pattern(prev_gut, prev_dung, prev_soil)]

## ---- bins defined from the ARG distribution ------------------------------
mkbin <- function(x, edges) cut(x, breaks = edges, include.lowest = TRUE, labels = FALSE)
qedges <- function(x, n) unique(quantile(x, seq(0, 1, length.out = n + 1), na.rm = TRUE))

E_MAG <- qedges(arg$mag, N_BIN_MAG); E_ABN <- qedges(arg$abn, N_BIN_ABN)
E_LEN <- qedges(arg$len, N_BIN_LEN)
# open the outer edges so pool genes beyond the ARG range fall in the end cells
open_out <- function(e) { e[1] <- -Inf; e[length(e)] <- Inf; e }

cell_of <- function(mag, abn, len, with_shape = NULL) {
  k <- paste(mkbin(mag, open_out(E_MAG)), mkbin(abn, open_out(E_ABN)),
             mkbin(len, open_out(E_LEN)), sep = ".")
  if (!is.null(with_shape)) k <- paste(k, with_shape, sep = ".")
  k
}

arg$cell  <- cell_of(arg$mag, arg$abn, arg$len)
arg$cellS <- cell_of(arg$mag, arg$abn, arg$len, arg$shape)
bg[, cell  := cell_of(mag, abn, len)]
bg[, cellS := cell_of(mag, abn, len, shape)]

## ---- draw ----------------------------------------------------------------
# without replacement within a draw, with replacement across draws
draw_matched <- function(cell_arg, cell_bg, gene_bg, label) {
  pool <- split(gene_bg, cell_bg)
  need <- table(cell_arg)
  avail <- vapply(names(need), function(k) length(pool[[k]]) %||% 0L, 1L)
  short <- need[avail < need]
  sets <- vector("list", B_DRAW)
  for (b in seq_len(B_DRAW)) {
    sets[[b]] <- unlist(lapply(names(need), function(k) {
      p <- pool[[k]]
      if (is.null(p)) return(character(0))
      sample(p, min(need[[k]], length(p)))
    }), use.names = FALSE)
  }
  list(sets = sets, need = need, avail = avail, short = short, label = label)
}

cat(sprintf("\ndrawing %d matched sets per arm ...\n", B_DRAW))
arm1 <- draw_matched(arg$cell,  bg$cell,  bg$gene_id, "arm1_magnitude")
arm1S <- draw_matched(arg$cellS, bg$cellS, bg$gene_id, "arm1_shape")

## ---- coverage and balance reporting ---------------------------------------
report_arm <- function(a, arg_df, cellcol) {
  n_need <- sum(a$need); n_got <- length(a$sets[[1]])
  cat(sprintf("\n--- %s ---\n", a$label))
  cat(sprintf("  ARG genes to match : %d\n", n_need))
  cat(sprintf("  matched per draw   : %d  (%.1f%% matched, %.1f%% UNMATCHED)\n",
              n_got, 100 * n_got / n_need, 100 * (n_need - n_got) / n_need))
  if (n_got != n_need)
    cat(sprintf("  feature counts differ (%d vs %d): SES scales with sqrt(N), so compare via raw containment\n",
                n_got, n_need))
  if (length(a$short)) {
    cat(sprintf("  cells short of donors: %d of %d\n", length(a$short), length(a$need)))
    # unmatched composition by occupancy shape
    um <- arg_df[arg_df[[cellcol]] %in% names(a$short), ]
    if (nrow(um)) {
      cat("  unmatched ARG genes by occupancy shape:\n")
      tb <- sort(table(um$shape), decreasing = TRUE)
      for (i in seq_along(tb))
        cat(sprintf("    %-3s %5d  (%4.1f%% of unmatched)\n",
                    names(tb)[i], tb[i], 100 * tb[i] / nrow(um)))
    }
  } else cat("  all cells supplied\n")

  # achieved balance on the three matched axes, first draw
  m <- bg[gene_id %in% a$sets[[1]]]
  cmp <- function(nm, av, bv) cat(sprintf(
    "  %-9s ARG median %8.4f (IQR %8.4f-%8.4f) | BG median %8.4f (IQR %8.4f-%8.4f)\n",
    nm, median(av), quantile(av, .25), quantile(av, .75),
    median(bv), quantile(bv, .25), quantile(bv, .75)))
  cat("  achieved balance (matched axes):\n")
  cmp("magnitude", arg_df$mag, m$mag); cmp("abundance", arg_df$abn, m$abn)
  cmp("length",    arg_df$len, m$len)
  cat("  occupancy shape (free in arm1):\n")
  sa <- prop.table(table(arg_df$shape)); sb <- prop.table(table(m$shape))
  for (k in union(names(sa), names(sb)))
    cat(sprintf("    %-3s ARG %5.1f%%  BG %5.1f%%\n", k,
                100 * (sa[k] %||% 0), 100 * (sb[k] %||% 0)))
  invisible(NULL)
}
report_arm(arm1,  arg, "cell")
report_arm(arm1S, arg, "cellS")

## ---- arm 2: family granularity ------------------------------------------
# Matched on the distribution of members per group
# ARG family-size distribution is highly skewed (median 1, 70% of symbols singletons) 
# so a median target would be satisfied by almost any grouping
fam <- read.delim(file.path(TAB_DIR, "T23b_arg_family_granularity.tsv"),
                  stringsAsFactors = FALSE)
target_sizes <- sort(fam$n_genes_detected, decreasing = TRUE)

cat("\n--- arm2_family: choosing the OG level ---\n")
pick_level <- function(col) {
  sz <- bg[, .N, by = c(col)][get(col) != "-"]
  q <- quantile(sz$N, c(.5, .75, .9, .99))
  cat(sprintf("  %-9s groups %7s | members/group median %3.0f q75 %4.0f q90 %5.0f q99 %6.0f\n",
              col, format(nrow(sz), big.mark = ","), q[1], q[2], q[3], q[4]))
  sz
}
sz_narrow <- pick_level("og_narrow"); sz_bact <- pick_level("og_bact")
q_arg <- quantile(target_sizes, c(.5, .75, .9, .99))
cat(sprintf("  %-9s groups %7d | members/group median %3.0f q75 %4.0f q90 %5.0f q99 %6.0f  <-- TARGET\n",
            "ARG symbol", length(target_sizes), q_arg[1], q_arg[2], q_arg[3], q_arg[4]))

# choose the level whose size distribution is closer to the ARG target on the log-quantile profile
dist_to_target <- function(sz) sum((log1p(quantile(sz$N, c(.5, .75, .9, .99))) -
                                    log1p(q_arg))^2)
lev <- if (dist_to_target(sz_narrow) <= dist_to_target(sz_bact)) "og_narrow" else "og_bact"
cat(sprintf("  -> using %s (log-quantile distance %.3f vs %.3f)\n", lev,
            dist_to_target(sz_narrow), dist_to_target(sz_bact)))

# For each ARG family of size k, take an OG with at least k members and subsample k of them:
# exact per-family size, exact total, and the family-size distribution reproduced
bg[, lev_og := get(lev)]
og_stat <- bg[lev_og != "-", .(N = .N, mmag = median(mag), mabn = median(abn),
                               mlen = median(log10(len))), by = lev_og]
setnames(og_stat, "lev_og", "og")
setorder(og_stat, N)

# ARG family centroids on the same three axes
fam_cent <- as.data.table(arg)[, .(k = .N, mmag = median(mag), mabn = median(abn),
                                   mlen = median(log10(len))), by = arg_symbol]
setorder(fam_cent, -k)

# Standardise so the three axes contribute comparably to the distance
sds <- c(mag = sd(og_stat$mmag), abn = sd(og_stat$mabn), len = sd(og_stat$mlen))
bg_by_og <- split(bg$gene_id, bg$lev_og)
bg_idx   <- split(seq_len(nrow(bg)), bg$lev_og)

draw_families <- function() {
  used <- character(0); out <- list(); genes <- character(0); shortfall <- 0L
  for (r in seq_len(nrow(fam_cent))) {
    k <- fam_cent$k[r]
    cand <- og_stat[N >= k & !og %in% used]
    if (!nrow(cand)) {
      cand <- og_stat[!og %in% used]
      if (!nrow(cand)) break
      cand <- cand[which.max(N)]
      shortfall <- shortfall + (k - cand$N)
    } else {
      d <- ((cand$mmag - fam_cent$mmag[r]) / sds["mag"])^2 +
           ((cand$mabn - fam_cent$mabn[r]) / sds["abn"])^2 +
           ((cand$mlen - fam_cent$mlen[r]) / sds["len"])^2
      # sample among the OG_POOL closest candidates (not which.min), so draws vary
      cand <- cand[sample(order(d)[seq_len(min(OG_POOL, length(d)))], 1)]
    }
    g <- cand$og; used <- c(used, g)
    take <- min(k, cand$N)
    # members: sample k from the closest 2k
    ii <- bg_idx[[g]]
    dm <- ((bg$mag[ii] - fam_cent$mmag[r]) / sds["mag"])^2 +
          ((bg$abn[ii] - fam_cent$mabn[r]) / sds["abn"])^2 +
          ((log10(bg$len[ii]) - fam_cent$mlen[r]) / sds["len"])^2
    near <- ii[order(dm)[seq_len(min(2L * take, length(ii)))]]
    # `near` holds row indices into bg; index through bg$gene_id
    genes <- c(genes, bg$gene_id[if (length(near) > take) sample(near, take) else near])
    out[[length(out) + 1]] <- data.frame(og = g, og_size = cand$N,
                                         target = k, taken = take,
                                         arg_symbol = fam_cent$arg_symbol[r],
                                         stringsAsFactors = FALSE)
  }
  list(map = do.call(rbind, out), genes = genes, shortfall = shortfall)
}

fam_sets <- vector("list", B_DRAW)
fd <- NULL
for (b in seq_len(B_DRAW)) {
  d <- draw_families()
  fam_sets[[b]] <- d$genes
  if (b == 1) { fd <- d$map; short1 <- d$shortfall }
}
fam_draw <- fd
n_take <- sum(fam_draw$taken); n_target <- sum(target_sizes)
cat(sprintf("  matched %d of %d ARG families\n", nrow(fam_draw), length(target_sizes)))
cat(sprintf("  members drawn: %s vs %s ARG alleles  (%s)\n",
            format(n_take, big.mark = ","), format(n_target, big.mark = ","),
            if (n_take == n_target) "exact"
            else sprintf("short by %d; use raw containment", n_target - n_take)))
if (short1 > 0)
  cat(sprintf("  %d members short across families with no OG large enough at this level\n",
              short1))
cat(sprintf("  achieved family-size profile vs target: median %.0f/%.0f  q90 %.0f/%.0f  max %.0f/%.0f\n",
            median(fam_draw$taken), median(target_sizes),
            quantile(fam_draw$taken, .9), quantile(target_sizes, .9),
            max(fam_draw$taken), max(target_sizes)))

## ---- write ---------------------------------------------------------------
save_sets <- function(a) {
  dt <- rbindlist(lapply(seq_along(a$sets), function(b)
    data.table(draw = b, gene_id = a$sets[[b]])))
  f <- file.path(TAB_DIR, sprintf("T25_%s_sets.tsv.gz", a$label))
  fwrite(dt, f, sep = "\t")
  cat(sprintf("wrote %s (%s rows)\n", basename(f), format(nrow(dt), big.mark = ",")))
}
cat("\n")
save_sets(arm1); save_sets(arm1S)
save_sets(list(sets = fam_sets, label = "arm2_family"))
save_tab(fam_draw, "T25_arm2_family_map.tsv")


# Union of every selected gene: input to the cluster abundance-extraction job
ids <- sort(unique(c(unlist(arm1$sets), unlist(arm1S$sets), unlist(fam_sets))))
writeLines(ids, file.path(TAB_DIR, "background_gene_ids.txt"))
cat(sprintf("wrote background_gene_ids.txt (%s genes)\n", format(length(ids), big.mark = ",")))
cat("\nBackground sets built.\n")
