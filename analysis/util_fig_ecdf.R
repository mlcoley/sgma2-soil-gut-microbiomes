#!/usr/bin/env Rscript
# ECDF of clinical ARG abundance by compartment. Computes no statistic beyond the
# separation check: it reads T_arg_abundance.tsv, so a caption or axis fix never
# re-runs 02_descriptive.R. The log axis is for the dung and soil distributions; it
# does not make the gut/environment gap (a 1.3% ratio) any more visible.
# Usage: Rscript analysis/util_fig_ecdf.R

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
if (!exists("TAB_DIR")) source("analysis/00_load_features.R")
if (!exists("logmsg")) logmsg <- function(...) cat(sprintf(...), "\n")

ECDF_BRK <- c(0.15, 0.3, 0.6, 1.2, 2.5, 5)

ecdf_theme <- function() theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ecdf_figure <- function(tab = file.path(TAB_DIR, "T_arg_abundance.tsv")) {
  abund <- read.delim(tab, stringsAsFactors = FALSE) %>%
    mutate(compartment = factor(compartment, levels = c("gut", "dung", "soil")),
           env = ifelse(compartment == "gut", "gut", "environment"))
  stopifnot(!any(is.na(abund$compartment)), nrow(abund) == 194, all(abund$clinical > 0))

  g   <- abund$clinical[abund$env == "gut"]
  e   <- abund$clinical[abund$env == "environment"]
  gap <- min(g) - max(e)
  pos <- mean(outer(g, e, ">"))
  logmsg("ECDF: gut n=%d min=%.3f | env n=%d max=%.3f | gap=%.4f | P(sup)=%.4f",
         length(g), min(g), length(e), max(e), gap, pos)

  d <- abund
  d$compartment <- factor(compartment_display(d$compartment),
                          levels = unname(COMPARTMENT_DISPLAY))
  cols <- setNames(unname(COMP_COLS[names(COMPARTMENT_DISPLAY)]),
                   unname(COMPARTMENT_DISPLAY))
  p <- ggplot(d, aes(clinical, color = compartment)) +
    annotate("rect", xmin = max(e), xmax = min(g), ymin = -Inf, ymax = Inf,
             fill = "grey40", alpha = 0.12) +
    stat_ecdf(geom = "step", linewidth = 0.8) +
    scale_x_log10(breaks = ECDF_BRK, labels = format(ECDF_BRK, drop0trailing = TRUE)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = cols, name = NULL) +
    labs(x = "Clinical ARG abundance per average genome (log scale)",
         y = "Samples at or below this abundance") +
    ecdf_theme() + theme(legend.position = "top")
  save_fig(p, "F3b_clinical_arg_abundance_ecdf.png", width = 7.2, height = 4.8)
  invisible(c(gap = gap, pos = pos))
}

if (identical(environment(), globalenv()) && sys.nframe() == 0) ecdf_figure()
