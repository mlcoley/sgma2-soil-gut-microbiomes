#!/usr/bin/env Rscript
# Matched background gene sets for the restriction to catalog genes assembled from a
# cattle dung or soil sample
# Inputs: as 19_background_match.R, and feature_tables/background_abundance*.tsv.gz.
# Output: T25_arm1_magnitude_nongut_sets.tsv.gz
# Usage: Rscript analysis/20_background_match_nongut.R   (from the repository root)

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
source("analysis/00_load_features.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(1)
B_DRAW    <- 100    # matched sets per arm
N_BIN_MAG <- 8      # quantile bins, occupancy magnitude
N_BIN_ABN <- 8      # quantile bins, log10 abundance
N_BIN_LEN <- 4      # quantile bins, gene length
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
# AMRFinderPlus genes and unannotated genes dropped by HPC job;
# remaining exclusions below matched on KO identifiers
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

# Only about half the pool carries a KO, so the screen above cannot see the rest.
# restriction to KO-annotated genes is analogous to the clinical ARG subset limited 
# to well-characterised determinants (element_type == "AMR" & scope == "core").
REQUIRE_KO <- TRUE

n0 <- nrow(bg)
if (REQUIRE_KO) {
  bg <- bg[ko != "-"]
  cat(sprintf("REQUIRE_KO: pool restricted to KO-annotated genes, %s dropped (%.1f%%)\n",
              format(n0 - nrow(bg), big.mark = ","), 100 * (n0 - nrow(bg)) / n0))
}
n1 <- nrow(bg)                     # after REQUIRE_KO, before semantic screen
bg <- bg[!(ko %in% EXCLUDE_KO)]
cat(sprintf("semantic KO exclusions (mobility/defence/metal): %s genes dropped (%.2f%% of the pre-screen pool)\n",
            format(n1 - nrow(bg), big.mark = ","), 100 * (n1 - nrow(bg)) / n1))
if (!REQUIRE_KO)
  cat(sprintf("  %.1f%% of the pool has no KO ('-') and is not screened\n",
              100 * mean(bg$ko == "-")))

bg <- bg[total_count > 0 & !is.na(gene_length) & gene_length > 0]

# Restrict ARG set and the donor pool to genes assembled from dung or soil.
full  <- load_key(biological_only = FALSE, soil_primary_only = FALSE)
scomp <- setNames(as.character(full$compartment), full$sample_id)
src_bg  <- scomp[sub("_k141.*$", "", bg$gene_id)]
src_arg <- scomp[sub("_k141.*$", "", arg$gene_id)]
stopifnot(!anyNA(src_bg), !anyNA(src_arg))
cat(sprintf("non-gut restriction: ARG set %d -> %d; pool %s -> %s\n", nrow(arg), sum(src_arg != "gut"),
            format(nrow(bg), big.mark = ","), format(sum(src_bg != "gut"), big.mark = ",")))
bg  <- bg[src_bg != "gut"]
arg <- arg[src_arg != "gut", ]
fs <- sort(Sys.glob(file.path(FT_DIR, "background_abundance*.tsv.gz")))
ids <- unique(unlist(lapply(fs, function(f) fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t",
                                                  select = "gene_id", showProgress = FALSE)$gene_id)))
n_before <- nrow(bg)
bg <- bg[gene_id %in% ids]
cat(sprintf("donor pool limited to genes with abundance rows: %s -> %s\n",
            format(n_before, big.mark = ","), format(nrow(bg), big.mark = ",")))
cat(sprintf("background pool after exclusions: %s genes\n", format(nrow(bg), big.mark = ",")))

## ---- axes ----------------------------------------------------------------
# magnitude as a rate per compartment, so soil samples do not dominate raw count
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

## ---- write ---------------------------------------------------------------
dt <- rbindlist(lapply(seq_along(arm1$sets), function(b) data.table(draw = b, gene_id = arm1$sets[[b]])))
f <- file.path(TAB_DIR, "T25_arm1_magnitude_nongut_sets.tsv.gz")
fwrite(dt, f, sep = "\t")
cat(sprintf("wrote %s (%s rows)\n", basename(f), format(nrow(dt), big.mark = ",")))
cat("\nNon-gut background sets built.\n")
