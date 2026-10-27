#!/usr/bin/env Rscript
# Figure builders for the pair classifier (Figure S7) and its SHAP attribution
# (Figure S8). They read tables and compute no statistic, so a caption or label fix
# never re-runs the permutation. Every panel plots a value against the null it must be
# read against.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

R8_ROLE_COLS <- c(primary = "#B5462F", control = "#3F6B57")
R8_STRAT_ORDER <- c("gut to soil", "gut to dung", "dung to soil", "soil to soil")

# Strip text only: the `stratum` values in T40c-T40e are join keys.
R8_STRAT_STRIP <- c("gut to soil"  = "child stool to soil",
                    "gut to dung"  = "child stool to cattle dung",
                    "dung to soil" = "cattle dung to soil",
                    "soil to soil" = "soil to soil")
r8_strat_labeller <- function() ggplot2::as_labeller(R8_STRAT_STRIP)

r8_theme <- function() {
  theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold", size = 8.5),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8.7, colour = "grey30"),
          plot.caption = element_text(size = 7.4, colour = "grey35", hjust = 0))
}

r8_figure_classifier <- function(t40d) {
  # specification on the y axis, where the labels have the full panel width
  spec_lab <- c(primary        = "PRIMARY\npresence + summaries",
                sens_nosummary = "sensitivity\npresence only",
                sens_abundance = "sensitivity\nabundance + summaries")
  keep <- intersect(names(spec_lab), unique(t40d$spec))

  obs <- t40d %>%
    filter(spec %in% keep) %>%
    mutate(stratum = factor(stratum, levels = intersect(R8_STRAT_ORDER, unique(stratum))),
           # reversed so the primary specification reads at the top of each panel
           spec = factor(spec, levels = rev(keep)),
           y    = as.numeric(spec),
           lab  = sprintf("AUC %.3f, p = %.3f", AUC, p),
           # flip the label left of its point when the point sits far right
           lab_h = ifelse(AUC > 0.62, 1.06, -0.06))

  p <- ggplot(obs, aes(x = AUC, y = y)) +
    geom_rect(aes(xmin = null_mean - 1.96 * null_sd, xmax = null_mean + 1.96 * null_sd,
                  ymin = y - 0.30, ymax = y + 0.30),
              fill = "grey85", colour = NA) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    # smallest AUC the design can detect at 80% power
    geom_segment(aes(x = detectable_AUC_80pct, xend = detectable_AUC_80pct,
                     y = y - 0.17, yend = y + 0.17),
                 colour = "grey30", linewidth = 0.45, linetype = "solid") +
    geom_point(aes(colour = role), size = 2.3) +
    geom_text(aes(label = lab, colour = role, hjust = lab_h), y = obs$y + 0.42,
              size = 2.4, show.legend = FALSE) +
    facet_wrap(~ stratum, nrow = 1, labeller = r8_strat_labeller()) +
    scale_y_continuous(breaks = seq_along(keep),
                       labels = spec_lab[rev(keep)],
                       limits = c(0.45, length(keep) + 0.70)) +
    scale_colour_manual(values = R8_ROLE_COLS, name = NULL,
                        labels = c(primary = "primary test", control = "positive control")) +
    expand_limits(x = c(0.42, 0.77)) +
    labs(x = "Mean AUC over repeated household-grouped cross-validation", y = NULL) +
    r8_theme() + theme(legend.position = "top",
                       axis.text.y = element_text(size = 7.6, lineheight = 0.95),
                       panel.grid.major.y = element_blank())
  save_fig(p, "F43_pair_classifier_vs_null.png", width = 12, height = 3.9)
  p
}


r8_figure_shap <- function(t40e, top_n = 15) {
  d <- t40e %>%
    group_by(stratum) %>%
    slice_max(mean_abs_shap, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(panel = stratum,
           panel = factor(panel, levels = unique(panel[order(
             match(stratum, R8_STRAT_ORDER))])),
           stratum = factor(stratum, levels = intersect(R8_STRAT_ORDER, unique(stratum))),
           feature = tidytext_reorder(feature, mean_abs_shap, stratum),
           agree = case_when(
             is.na(univariate_z) ~ "summary column (no univariate analogue)",
             abs(univariate_z) >= 2 ~ "agrees with univariate test",
             TRUE ~ "model only"))

  p <- ggplot(d, aes(x = mean_abs_shap, y = feature)) +
    geom_col(aes(fill = agree), width = 0.72) +
    geom_point(aes(x = null_shap_p95), shape = 18, size = 3.1, colour = "black") +
    facet_wrap(~ panel, scales = "free", ncol = 2,
               labeller = r8_strat_labeller()) +
    scale_fill_manual(values = c("agrees with univariate test" = "#3F6B57",
                                 "model only" = "#C08A2A",
                                 "summary column (no univariate analogue)" = "#6C7A89"),
                      name = NULL) +
    # shorten gene ids (P026_k141_85061_2 -> P026_2) and summary-column names
    scale_y_discrete(labels = function(v) {
      x <- sub("___.*$", "", v)
      x <- sub("^([^_]+)_[^_]+_[^_]+_([^_]+)$", "\\1_\\2", x)
      sub("^shared_count_", "shared_", x)
    }) +
    labs(x = "mean absolute SHAP contribution", y = NULL) +
    r8_theme() + theme(legend.position = "top",
                       axis.text.y = element_text(size = 8))
  save_fig(p, "F44_shap_attribution.png", width = 11, height = 8.5)
  p
}


# Small helper so features can be ordered within each facet without pulling in
# the tidytext dependency.
tidytext_reorder <- function(x, by, within) {
  factor(paste(x, within, sep = "___"),
         levels = unique(paste(x, within, sep = "___")[order(within, by)]))
}

