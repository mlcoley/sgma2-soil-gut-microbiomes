#!/usr/bin/env Rscript
# Per-genome normalization: MicrobeCensus genome equivalents (GE) against the
# single-copy marker proxy.
# Clinical ARG abundance is reported as copies per average genome. Two denominators
# exist: (A) the median count of 57 bacterial ribosomal-protein marker KOs, measured
# through the gene catalog, and (B) MicrobeCensus GE = total bases / average genome
# size, computed from reads and independent of the catalog. The ARG numerator is
# catalog-mapped under either choice, so catalog capture cancels in the ratio only
# for (A). Compartment rankings and the per-compartment tilt of the two factors are
# compared, and the tilt is related to a catalog-capture proxy (reads on KO-carrying
# genes / post-host reads).
# Outputs: T11a-T11d.
# Usage: Rscript analysis/06_normalization.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/00_load_features.R")

set.seed(1)
GATING <- load_gating()

log_lines <- c(sprintf("Normalization comparison (GE vs single-copy) - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

key <- load_key()
logmsg("cohort: %d samples (gut %d / dung %d / soil %d)", nrow(key),
       sum(key$compartment == "gut"), sum(key$compartment == "dung"),
       sum(key$compartment == "soil"))

## ---- the two factors -----------------------------------------------------
ge_f <- file.path(FT_DIR, "microbecensus_ge.tsv")
if (!file.exists(ge_f))
  stop("microbecensus_ge.tsv not present — nothing to compare")
ge <- read.delim(ge_f, stringsAsFactors = FALSE)

ko  <- get_layer("ko_all", key)                       # KO is not filtered
mk  <- intersect(load_single_copy_markers(), colnames(ko$M))
sc  <- apply(ko$M[, mk, drop = FALSE], 1, median)     # single-copy factor, per sample
logmsg("single-copy factor from %d of %d marker KOs present in the KO table",
       length(mk), length(load_single_copy_markers()))

miss <- setdiff(key$sample_id, ge$sample_id)
if (length(miss))
  logmsg("WARNING: %d cohort sample(s) absent from the GE table: %s",
         length(miss), paste(miss, collapse = ", "))

## ---- catalog-capture proxy (explanatory covariate only) -------------------
acc <- read.delim(file.path(QC_DIR, "sgma2_step1_read_accounting.tsv"), stringsAsFactors = FALSE)

d <- data.frame(sample_id = key$sample_id, compartment = key$compartment,
                stringsAsFactors = FALSE)
d$ags          <- ge$average_genome_size[match(d$sample_id, ge$sample_id)]
d$ge           <- ge$genome_equivalents[match(d$sample_id, ge$sample_id)]
d$single_copy  <- sc[d$sample_id]
d$posthost_reads <- acc$posthost_reads[match(d$sample_id, acc$sample_id)]
d$posthost_bases <- acc$posthost_bases[match(d$sample_id, acc$sample_id)]
d$ko_reads     <- rowSums(ko$M)[d$sample_id]
d$capture      <- d$ko_reads / d$posthost_reads
d <- d[!is.na(d$ge) & !is.na(d$single_copy), ]

## ---- independent validation of the GE table ------------------------------
# GE = total_bases / AGS, so GE * AGS must reproduce the library's base count.
d$implied_bases <- d$ge * d$ags
bases_r <- suppressWarnings(cor(log10(d$implied_bases), log10(d$posthost_bases),
                                use = "complete.obs"))
bases_dev <- median(d$implied_bases / d$posthost_bases, na.rm = TRUE)
logmsg("GE table validation: implied bases vs post-host bases  log10 r = %.4f, median ratio = %.3f",
       bases_r, bases_dev)
logmsg("  (ratio 1.000 = MicrobeCensus counted exactly the post-host libraries; independent of our pipeline)")

## ---- agreement between the two factors -----------------------------------
cor_row <- function(df, label) {
  data.frame(subset = label, n = nrow(df),
             spearman   = cor(df$single_copy, df$ge, method = "spearman"),
             pearson_log= cor(log10(df$single_copy), log10(df$ge)),
             stringsAsFactors = FALSE)
}
cors <- rbind(cor_row(d, "all"),
              do.call(rbind, lapply(COMP_LEVELS,
                                    function(c) cor_row(d[d$compartment == c, ], c))))
for (i in seq_len(nrow(cors)))
  logmsg("agreement [%-4s] n=%3d  Spearman %.3f  Pearson(log10) %.3f",
         cors$subset[i], cors$n[i], cors$spearman[i], cors$pearson_log[i])

## ---- the compartment-bias test -------------------------------------------
# The two factors are on different scales, so the ratio single-copy / GE is scaled
# to the cohort median (1.0 = normalized identically by both).
d$ratio     <- d$single_copy / d$ge
d$ratio_rel <- d$ratio / median(d$ratio)

kw <- kruskal.test(ratio_rel ~ compartment, data = d)
logmsg("")
logmsg("compartment tilt in single-copy/GE: Kruskal-Wallis chi2 = %.1f, df = %d, p = %.3g",
       kw$statistic, kw$parameter, kw$p.value)

comp_sum <- d %>% group_by(compartment) %>%
  summarise(n = n(),
            ags_median      = median(ags),
            ge_median       = median(ge),
            sc_median       = median(single_copy),
            capture_median  = median(capture),
            ratio_rel_median= median(ratio_rel),
            ratio_rel_iqr   = IQR(ratio_rel), .groups = "drop") %>%
  as.data.frame()
for (i in seq_len(nrow(comp_sum)))
  logmsg("  %-4s n=%3d  AGS %.2f Mbp  GE %8.0f  single-copy %8.0f  capture %.3f  ratio_rel %.2f",
         comp_sum$compartment[i], comp_sum$n[i], comp_sum$ags_median[i]/1e6,
         comp_sum$ge_median[i], comp_sum$sc_median[i], comp_sum$capture_median[i],
         comp_sum$ratio_rel_median[i])

pw <- pairwise.wilcox.test(d$ratio_rel, d$compartment, p.adjust.method = "BH")
logmsg("pairwise Wilcoxon on ratio_rel (BH):")
pwdf <- as.data.frame(as.table(pw$p.value)); names(pwdf) <- c("a", "b", "p_BH")
pwdf <- pwdf[!is.na(pwdf$p_BH), ]
for (i in seq_len(nrow(pwdf)))
  logmsg("  %s vs %s  BH-adjusted p = %.3g", pwdf$a[i], pwdf$b[i], pwdf$p_BH[i])

tilt <- max(comp_sum$ratio_rel_median) / min(comp_sum$ratio_rel_median)
logmsg("MAX COMPARTMENT TILT = %.2fx (%s vs %s)", tilt,
       comp_sum$compartment[which.max(comp_sum$ratio_rel_median)],
       comp_sum$compartment[which.min(comp_sum$ratio_rel_median)])

## ---- does the tilt track catalog capture? --------------------------------
cap_r  <- cor(log10(d$ratio_rel), log10(d$capture), method = "pearson")
cap_rs <- cor(d$ratio_rel, d$capture, method = "spearman")
logmsg("")
logmsg("ratio_rel vs catalog-capture proxy: Pearson(log10) %.3f, Spearman %.3f",
       cap_r, cap_rs)
for (cc in COMP_LEVELS) {
  s <- d[d$compartment == cc, ]
  logmsg("  within %-4s Spearman %.3f (n=%d)", cc,
         cor(s$ratio_rel, s$capture, method = "spearman"), nrow(s))
}

## ---- consequence for clinical ARG abundance --------------------------------
# Units differ between denominators, so compare the between-compartment ratios.
arg  <- gate_layer(get_layer("arg", key), GATING)
clin <- intersect(clinical_args(), colnames(arg$M))
clin_counts <- rowSums(arg$M[, clin, drop = FALSE])
common <- intersect(names(clin_counts), d$sample_id)
dd <- d[match(common, d$sample_id), ]
dd$clin_sc <- clin_counts[common] / dd$single_copy
dd$clin_ge <- clin_counts[common] / dd$ge

cons <- dd %>% group_by(compartment) %>%
  summarise(clinical_per_genome_singlecopy = median(clin_sc),
            clinical_per_genome_GE         = median(clin_ge), .groups = "drop") %>%
  as.data.frame()
logmsg("")
logmsg("clinical ARG copies per genome (median, %s):", gating_col(GATING))
for (i in seq_len(nrow(cons)))
  logmsg("  %-4s  single-copy %.3f   GE %.3f", cons$compartment[i],
         cons$clinical_per_genome_singlecopy[i], cons$clinical_per_genome_GE[i])

fold <- function(v, a, b) v[cons$compartment == a] / v[cons$compartment == b]
gs_sc <- fold(cons$clinical_per_genome_singlecopy, "gut", "soil")
gs_ge <- fold(cons$clinical_per_genome_GE,         "gut", "soil")
gd_sc <- fold(cons$clinical_per_genome_singlecopy, "gut", "dung")
gd_ge <- fold(cons$clinical_per_genome_GE,         "gut", "dung")
logmsg("gut/soil fold: single-copy %.2fx -> GE %.2fx   (change %.2fx)", gs_sc, gs_ge, gs_ge/gs_sc)
logmsg("gut/dung fold: single-copy %.2fx -> GE %.2fx   (change %.2fx)", gd_sc, gd_ge, gd_ge/gd_sc)
rank_sc <- cons$compartment[order(-cons$clinical_per_genome_singlecopy)]
rank_ge <- cons$compartment[order(-cons$clinical_per_genome_GE)]
logmsg("compartment ranking  single-copy: %s | GE: %s | preserved: %s",
       paste(rank_sc, collapse = ">"), paste(rank_ge, collapse = ">"),
       identical(as.character(rank_sc), as.character(rank_ge)))

## ---- tables --------------------------------------------------------------
save_tab(d[, c("sample_id", "compartment", "ags", "ge", "single_copy",
               "posthost_reads", "ko_reads", "capture", "ratio", "ratio_rel")],
         "T11a_normalization_per_sample.tsv")
save_tab(comp_sum, "T11b_normalization_by_compartment.tsv")
save_tab(cons,     "T11c_normalization_consequence.tsv")
save_tab(cors,     "T11d_normalization_agreement.tsv")

writeLines(unlist(log_lines), file.path(LOG_DIR, "06_normalization.log"))
cat("\nNormalization comparison complete.\n")
