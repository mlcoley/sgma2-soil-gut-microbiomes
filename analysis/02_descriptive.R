#!/usr/bin/env Rscript
# Descriptive overview: cohort inventory, clinical resistance-gene abundance by
# compartment (with tests and a contamination-filter sensitivity), and alpha
# diversity by compartment.
# Usage: Rscript analysis/02_descriptive.R   (from the repository root)

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(scales)
})
source("analysis/00_load_features.R")
set.seed(1)

log_lines <- c(sprintf("Descriptive overview - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines) + 1]] <<- m }

theme_sgma <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
key <- load_key()

# Negative-control filter, applied to every ARG-derived output here. Tables whose
# `layer` column mixes KO and ARG rows are labelled per row (KO is not filtered).
GATING <- load_gating()
stamp_gating <- function(df) {
  df$contamination_filter <- if ("layer" %in% names(df)) {
    ifelse(grepl("ARG", df$layer), gating_col(GATING), "n/a (KO not filtered)")
  } else gating_col(GATING)
  df
}

# ===========================================================================
# Inventory
# ===========================================================================
logmsg("== inventory ==")

# Key samples absent from the gene-catalog feature tables are excluded from all
# feature analyses.
ko_samps <- rownames(load_matrix(file.path(FT_DIR, LAYERS$ko_all$file), "woltka"))
missing  <- setdiff(key$sample_id, ko_samps)
if (length(missing))
  logmsg("NOTE: %d key sample(s) absent from the feature tables: %s",
         length(missing), paste(missing, collapse = ", "))
key_an <- key[!key$sample_id %in% missing, ]

# T1a: per-compartment counts and read depth (total post-host reads, R1 + R2).
# n_key = biological samples in the key; n_analyzed = those in the feature tables.
inv <- key %>% group_by(compartment) %>%
  summarise(n_key = n(), n_households_key = n_distinct(household), .groups = "drop") %>%
  left_join(key_an %>% group_by(compartment) %>%
              summarise(n_analyzed    = n(),
                        n_households  = n_distinct(household),
                        depth_min     = min(posthost_reads),
                        depth_mean    = mean(posthost_reads),
                        depth_median  = median(posthost_reads),
                        depth_max     = max(posthost_reads),
                        .groups = "drop"),
            by = "compartment") %>%
  select(compartment, n_key, n_analyzed, n_households,
         depth_min, depth_mean, depth_median, depth_max)
save_tab(inv, "T1a_compartment_inventory.tsv")
logmsg("compartments (key -> analyzed): %s",
       paste(sprintf("%s %d->%d", inv$compartment, inv$n_key, inv$n_analyzed), collapse = " | "))
logmsg("analyzed cohort total n = %d", sum(inv$n_analyzed))

# T1b: household x compartment completeness on the analyzed cohort.
key_only_by_hh <- key %>% filter(sample_id %in% missing) %>%
  group_by(household) %>%
  summarise(key_only = paste(sample_id, collapse = ","), .groups = "drop")
hh_comp <- key_an %>%
  count(household, compartment, .drop = FALSE) %>%
  pivot_wider(names_from = compartment, values_from = n, values_fill = 0) %>%
  mutate(complete_triad = (gut > 0) & (dung > 0) & (soil > 0)) %>%
  left_join(key_only_by_hh, by = "household") %>%
  mutate(key_only = ifelse(is.na(key_only), "", key_only)) %>%
  arrange(household)
save_tab(hh_comp, "T1b_household_completeness.tsv")
logmsg("households total = %d ; complete analyzable gut+dung+soil triad = %d",
       nrow(hh_comp), sum(hh_comp$complete_triad))

# ===========================================================================
# Clinical resistance-gene abundance by compartment
# ===========================================================================
logmsg("== resistome abundance ==")

ko  <- get_layer("ko_all", key)
arg <- gate_layer(get_layer("arg", key), GATING)
logmsg("ARG filter: %s%s", GATING$note,
       if (GATING$applied) sprintf(" -> %d ARG features retained", ncol(arg$M)) else "")
ann <- load_arg_annotation()

# Abundance in copies per average genome (denominator set by NORM_METHOD).
gf     <- genome_factor(key)
arg_pg <- per_genome(arg$M, gf)
catv   <- arg_category(ann)
res_syms  <- intersect(names(catv)[catv != "Virulence"], colnames(arg_pg))
clin_syms <- intersect(clinical_args(ann), colnames(arg_pg))
logmsg("resistome abundance normalization: %s", gf$method)

abund <- data.frame(sample_id = rownames(arg_pg), compartment = arg$key$compartment,
                    total = rowSums(arg_pg[, res_syms, drop = FALSE]),
                    clinical = rowSums(arg_pg[, clin_syms, drop = FALSE]))
save_tab(stamp_gating(abund), "T_arg_abundance.tsv")

# Probability of superiority P(a > b), ties counted as half.
p_superiority <- function(x, y) mean(outer(x, y, ">")) + 0.5 * mean(outer(x, y, "=="))

clin_ab <- data.frame(compartment = abund$compartment, cpg = abund$clinical)
clin_summary <- clin_ab %>%
  group_by(compartment) %>%
  summarise(n = n(),
            median = round(median(cpg), 4), q1 = round(quantile(cpg, .25), 4),
            q3 = round(quantile(cpg, .75), 4), iqr = round(IQR(cpg), 4),
            min = round(min(cpg), 4), max = round(max(cpg), 4), .groups = "drop") %>%
  mutate(quantity = "clinical ARG copies per average genome", denominator = gf$method)
save_tab(stamp_gating(clin_summary), "T_clinical_arg_abundance_summary.tsv")

# Kruskal-Wallis omnibus, then pairwise Wilcoxon rank-sum with BH adjustment.
kw_ab <- kruskal.test(cpg ~ compartment, data = clin_ab)
pw_ab <- pairwise.wilcox.test(clin_ab$cpg, clin_ab$compartment, p.adjust.method = "BH")
pairs_ab <- t(combn(levels(droplevels(clin_ab$compartment)), 2))
tests_ab <- do.call(rbind, lapply(seq_len(nrow(pairs_ab)), function(i) {
  a <- pairs_ab[i, 1]; b <- pairs_ab[i, 2]
  x <- clin_ab$cpg[clin_ab$compartment == a]; y <- clin_ab$cpg[clin_ab$compartment == b]
  w <- suppressWarnings(wilcox.test(x, y))
  data.frame(contrast = sprintf("%s vs %s", a, b), test = "Wilcoxon rank-sum",
             statistic = unname(w$statistic), p = w$p.value,
             p_adj_BH = pw_ab$p.value[b, a],
             median_ratio = round(median(x) / median(y), 3),
             p_superiority = round(p_superiority(x, y), 4))
}))
tests_ab <- bind_rows(
  data.frame(contrast = "compartment (omnibus)", test = "Kruskal-Wallis",
             statistic = unname(kw_ab$statistic), p = kw_ab$p.value,
             p_adj_BH = NA_real_, median_ratio = NA_real_, p_superiority = NA_real_),
  tests_ab)
tests_ab$quantity <- "clinical ARG copies per average genome"
tests_ab$denominator <- gf$method
save_tab(stamp_gating(tests_ab), "T_clinical_arg_abundance_tests.tsv")

for (i in seq_len(nrow(clin_summary)))
  logmsg("clinical ARG copies/genome [%-4s] median %.3f (IQR %.3f-%.3f, n=%d)",
         clin_summary$compartment[i], clin_summary$median[i],
         clin_summary$q1[i], clin_summary$q3[i], clin_summary$n[i])
logmsg("clinical ARG copies/genome: Kruskal-Wallis chi2=%.1f, df=%d, p=%.2e",
       kw_ab$statistic, kw_ab$parameter, kw_ab$p.value)
for (i in which(tests_ab$test == "Wilcoxon rank-sum"))
  logmsg("clinical ARG copies/genome [%s]: p=%.2e, BH-adjusted p=%.2e, P(superiority)=%.4f",
         tests_ab$contrast[i], tests_ab$p[i], tests_ab$p_adj_BH[i], tests_ab$p_superiority[i])

# ECDF of clinical ARG abundance by compartment (reads T_arg_abundance.tsv).
source("analysis/util_fig_ecdf.R")
ecdf_figure()

# Sensitivity to the contamination filter: abundance with and without it.
arg_un    <- get_layer("arg", key)
arg_pg_un <- per_genome(arg_un$M, gf)
res_un    <- intersect(names(catv)[catv != "Virulence"], colnames(arg_pg_un))
clin_all  <- clinical_args(ann)
clin_g    <- intersect(clin_all, colnames(arg_pg))
clin_u    <- intersect(clin_all, colnames(arg_pg_un))

gate_sens <- data.frame(
  sample_id   = rownames(arg_pg),
  compartment = arg$key$compartment,
  total_filtered    = rowSums(arg_pg[,    res_un[res_un %in% colnames(arg_pg)], drop = FALSE]),
  total_unfiltered  = rowSums(arg_pg_un[, res_un,  drop = FALSE])[rownames(arg_pg)],
  clinical_filtered   = rowSums(arg_pg[,    clin_g, drop = FALSE]),
  clinical_unfiltered = rowSums(arg_pg_un[, clin_u, drop = FALSE])[rownames(arg_pg)]
)
gate_sum <- gate_sens %>% group_by(compartment) %>%
  summarise(n = n(),
            total_filtered      = median(total_filtered),
            total_unfiltered    = median(total_unfiltered),
            clinical_filtered   = median(clinical_filtered),
            clinical_unfiltered = median(clinical_unfiltered),
            clinical_ratio   = median(clinical_unfiltered) / median(clinical_filtered),
            .groups = "drop") %>% as.data.frame()
save_tab(gate_sens, "T_arg_abundance_filter_sensitivity_per_sample.tsv")
save_tab(gate_sum,  "T_arg_abundance_filter_sensitivity.tsv")
for (i in seq_len(nrow(gate_sum)))
  logmsg("filter sensitivity [%-4s] clinical copies/genome: filtered %.3f vs unfiltered %.3f (%.2fx)",
         gate_sum$compartment[i], gate_sum$clinical_filtered[i],
         gate_sum$clinical_unfiltered[i], gate_sum$clinical_ratio[i])
gsf <- function(v) v[gate_sum$compartment == "gut"] / v[gate_sum$compartment == "soil"]
logmsg("filter sensitivity: gut/soil clinical contrast filtered %.2fx vs unfiltered %.2fx",
       gsf(gate_sum$clinical_filtered), gsf(gate_sum$clinical_unfiltered))

# ===========================================================================
# Alpha diversity by compartment
# ===========================================================================
logmsg("== alpha diversity ==")

alpha_of <- function(layer, label) {
  M <- layer$M
  data.frame(sample_id  = rownames(M),
             compartment = layer$key$compartment,
             household   = layer$key$household,
             layer       = label,
             richness    = rowSums(M > 0),
             shannon     = diversity(tss(M), "shannon"))
}
al <- rbind(alpha_of(ko, "KO"), alpha_of(arg, "ARG"))
save_tab(stamp_gating(al), "T_alpha_diversity.tsv")

# Kruskal-Wallis omnibus and BH-adjusted pairwise Wilcoxon per layer and metric,
# with a Fligner-Killeen test of spread (and per-compartment MAD and IQR) beside it.
kw_lines <- c(); kw_stats <- list(); alpha_disp <- list()
for (lay in unique(al$layer)) for (metric in c("richness", "shannon")) {
  d <- al[al$layer == lay, ]
  grp <- droplevels(d$compartment)
  kw <- kruskal.test(d[[metric]] ~ d$compartment)
  pw <- pairwise.wilcox.test(d[[metric]], d$compartment, p.adjust.method = "BH")
  fk <- fligner.test(d[[metric]], grp)
  pwm <- pw$p.value
  gp <- function(a, b) {
    v <- tryCatch(pwm[a, b], error = function(e) NA_real_)
    if (is.na(v)) v <- tryCatch(pwm[b, a], error = function(e) NA_real_)
    unname(v)
  }
  kw_stats[[length(kw_stats) + 1]] <- data.frame(
    layer = lay, metric = metric, chi2 = unname(kw$statistic),
    df = unname(kw$parameter), p = kw$p.value,
    p_gut_dung = gp("dung", "gut"), p_gut_soil = gp("soil", "gut"),
    p_dung_soil = gp("soil", "dung"))
  alpha_disp[[length(alpha_disp) + 1]] <- data.frame(
    layer = lay, metric = metric, test = "Fligner-Killeen (median-centred)",
    chi2 = unname(fk$statistic), df = unname(fk$parameter), p = fk$p.value,
    setNames(as.data.frame(t(round(tapply(d[[metric]], grp, mad), 4))),
             paste0("mad_", levels(grp))),
    setNames(as.data.frame(t(round(tapply(d[[metric]], grp, IQR), 4))),
             paste0("iqr_", levels(grp))),
    row.names = NULL)
  kw_lines <- c(kw_lines, sprintf("[%s | %s] Kruskal-Wallis chi2=%.1f, p=%.2e",
                                  lay, metric, kw$statistic, kw$p.value),
                capture.output(print(pw$p.value)),
                sprintf("[%s | %s] Fligner-Killeen (spread) chi2=%.1f, p=%.2e",
                        lay, metric, fk$statistic, fk$p.value), "")
}
kw_stats <- bind_rows(kw_stats)
alpha_disp <- bind_rows(alpha_disp)
writeLines(kw_lines, file.path(TAB_DIR, "T_alpha_tests.txt"))
save_tab(stamp_gating(kw_stats), "T_alpha_location_tests.tsv")
save_tab(stamp_gating(alpha_disp), "T_alpha_dispersion.tsv")
logmsg("alpha spread (Fligner-Killeen) differs by compartment in %d of %d layer x metric tests",
       sum(alpha_disp$p < 0.05), nrow(alpha_disp))

# Figure 2: four panels (KO and ARG family; richness and Shannon) with the
# Kruskal-Wallis omnibus read from `kw_stats`.
LAY_TITLE <- setNames(layer_display(c("KO", "ARG"), wrap = FALSE), c("KO", "ARG"))
alpha_panel <- function(lay, metric, ylab) {
  d  <- data.frame(compartment = al$compartment[al$layer == lay],
                   value       = al[[metric]][al$layer == lay])
  st <- kw_stats[kw_stats$layer == lay & kw_stats$metric == metric, ]
  ggplot(d, aes(compartment, value, color = compartment)) +
    geom_boxplot(width = .55, outlier.shape = NA, fill = NA) +
    geom_jitter(width = .15, alpha = .4, size = .9) +
    annotate("text", x = Inf, y = Inf, hjust = 1.04, vjust = 1.5, size = 3, color = "grey25",
             label = sprintf("KW χ² = %.1f, %s", st$chi2[1], fmt_p(st$p[1]))) +
    scale_color_manual(values = COMP_COLS, guide = "none") +
    scale_x_discrete(labels = compartment_display) +
    scale_y_continuous(expand = expansion(mult = c(.05, .15))) +
    labs(title = LAY_TITLE[[lay]], x = NULL, y = ylab) +
    theme_sgma +
    theme(plot.title = element_text(face = "plain", size = 10))
}

f5 <- patchwork::wrap_plots(
  list(alpha_panel("KO",  "richness", "Observed richness"),
       alpha_panel("KO",  "shannon",  "Shannon"),
       alpha_panel("ARG", "richness", "Observed richness"),
       alpha_panel("ARG", "shannon",  "Shannon")),
  nrow = 2)
save_fig(f5, "F5_alpha_by_compartment.png", width = 8, height = 5.2)

writeLines(log_lines, file.path(LOG_DIR, "02_descriptive.log"))
cat("\nDescriptive overview complete.\n")
