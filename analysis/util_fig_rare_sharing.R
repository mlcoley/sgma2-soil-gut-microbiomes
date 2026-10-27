#!/usr/bin/env Rscript
# Figure 5 (main) and Figure S3 (supplement): within-household excess similarity for
# the rare-feature sweep. Computes no statistic; draws from T13.
#
# A point is one whole test (the mean within-household similarity across all 38-39
# households, in standard deviations from the mean of 999 household-label
# permutations), not one sample pair. The prevalence cap is in facet columns, since it
# is four categories, and the compartment pair is on the x axis. BH is applied across
# the full 120-test sweep; the main figure shows only the two ARG-variant layers
# (48 tests) and uses the full-sweep adjusted p-values, not a recomputation on the
# subset.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
if (!exists("logmsg")) logmsg <- function(...) cat(sprintf(...), "\n")

F17_LAYER_ORDER <- c("ARG (all)", "ARG (clinical)", "ARG allele (all)",
                     "ARG allele (clinical)", "KO (function)")
F17_VARIANT     <- c("ARG allele (all)", "ARG allele (clinical)")
F17_PAIR_SHORT  <- c(`dung-soil` = "cattle dung\nsoil",
                     `gut-dung` = "child stool\ncattle dung",
                     `gut-soil` = "child stool\nsoil")

f17_prep <- function(t13) {
  t13 %>% mutate(
    layer  = factor(layer_display(layer), levels = layer_display(F17_LAYER_ORDER)),
    caplab = factor(sprintf("features in ≤%d%% of households", 100 * prev_cap),
                    levels = sprintf("features in ≤%d%% of households", c(100, 50, 25, 10))),
    pairx  = factor(F17_PAIR_SHORT[pair],
                    levels = F17_PAIR_SHORT[c("dung-soil", "gut-dung", "gut-soil")]),
    metric = factor(ifelse(grepl("Jaccard", metric), "Jaccard (presence)",
                           "Bray-Curtis (abundance)"),
                    levels = c("Jaccard (presence)", "Bray-Curtis (abundance)")),
    BH     = ifelse(q_BH < 0.05, "survives BH", "no"))
}

# n_tests is passed in because the main figure draws a subset while BH was applied
# across the full sweep.
f17_draw <- function(d, n_tests) {
  nf <- d %>% distinct(layer, caplab, n_features)
  ggplot(d, aes(pairx, z)) +
    geom_hline(yintercept = 0, linewidth = .25, colour = "grey35") +
    geom_hline(yintercept = 1.96, linetype = 2, linewidth = .3, colour = "grey60") +
    geom_point(aes(fill = metric, shape = BH), colour = "grey55", stroke = .25,
               size = 2.2, position = position_dodge(width = .55)) +
    geom_text(data = nf, aes(x = 0.55, y = 6.6,
                             label = sprintf("%s features", format(n_features, big.mark = ","))),
              inherit.aes = FALSE, hjust = 0, size = 2.75, colour = "grey45") +
    geom_text(data = transform(nf[nf$caplab == levels(nf$caplab)[1], , drop = FALSE],
                               xx = 0.55, yy = 1.96, lab = "1.96"),
              aes(x = xx, y = yy, label = lab), inherit.aes = FALSE,
              hjust = 0, vjust = -0.45, size = 2.5, colour = "grey55") +
    facet_grid(layer ~ caplab,
               labeller = labeller(layer = function(x) gsub(" ", "\n", x))) +
    scale_fill_manual(values = c("Jaccard (presence)" = "#3F6B57",
                                 "Bray-Curtis (abundance)" = "#C08A2A"), name = NULL) +
    scale_shape_manual(values = c("survives BH" = 24, "no" = 21), name = NULL,
                       breaks = "survives BH",
                       labels = sprintf("survives BH correction across %d tests", n_tests)) +
    coord_cartesian(ylim = c(-2.6, 7.1)) +
    labs(x = "compartment pair",
         y = "standard deviations from the\nhousehold-permuted expectation (z)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.spacing.y = unit(0, "lines"),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.text.x = element_text(size = 9.5),
          strip.text.y = element_text(size = 9, angle = 0),
          axis.text.x = element_text(size = 8.5, lineheight = .9),
          axis.title.x = element_text(margin = margin(t = 9))) +
    guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, colour = "grey55", stroke = .25, size = 2.8)),
           shape = guide_legend(order = 2, override.aes = list(fill = "grey55", colour = "grey55", stroke = .25, size = 2.8)))
}

# Figure 5: the two ARG-variant layers.
f17a_main <- function(t13) {
  d <- f17_prep(t13[t13$layer %in% F17_VARIANT, , drop = FALSE])
  f17_draw(d, n_tests = nrow(t13))
}

# Figure S3: all five layers. At tight caps the feature space is small enough that z
# and the empirical p disagree at the 5% level; the check below stops the figure if
# such a cell ever falls in a variant layer or survives BH, since the figure does not
# mark them.
f17b_supplement <- function(t13) {
  n_clear <- sum(t13$z > 1.96)
  n_bh    <- sum(t13$q_BH < 0.05)
  flagged <- t13[(t13$z > 1.96) != (t13$p < 0.025), , drop = FALSE]
  stopifnot(!any(flagged$layer %in% F17_VARIANT), !any(flagged$q_BH < 0.05))
  logmsg(paste0(
    "Fig S3 check: %d of %d points clear z = 1.96, %d survive BH; ",
    "%d calibration-flagged cells (BH-adjusted p %.3f-%.3f), none in an ARG ",
    "variant layer and none a BH survivor"),
    n_clear, nrow(t13), n_bh, nrow(flagged),
    min(flagged$q_BH), max(flagged$q_BH))
  f17_draw(f17_prep(t13), n_tests = nrow(t13))
}
