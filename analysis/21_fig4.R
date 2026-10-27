#!/usr/bin/env Rscript
# Figure 4: clinical ARG variants against the matched non-ARG background, drawn from
# T26a (20_background_containment.R). Computes no statistic.
# Observed containment (left) and distance from the chance expectation for a set of
# that size (right) are drawn side by side because neither is readable alone:
# C(gut -> dung) is +4.81 on the standardized measure and -0.16 on observed
# containment. C(gut -> env) is dropped from the figure: it shares a numerator with
# C(env -> gut), so its standardized score is algebraically identical, and its observed
# value differs only by the set-size ratio (C(A->B) / C(B->A) = |B| / |A|).
# The background is drawn as a mean with an error bar (the 100 matched draws are
# independent draws, not a continuous envelope). Zero is marked on the standardized
# panel: it is the chance expectation, and only C(gut -> dung) sits above it.
# Usage: Rscript analysis/21_fig4.R   (from the repository root)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales); library(patchwork)
})
source("analysis/00_load_features.R")

theme_sgma <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

read_tab <- function(f) read.delim(file.path(TAB_DIR, f), stringsAsFactors = FALSE,
                                   check.names = FALSE)
t26a <- read_tab("T26a_background_containment.tsv")

PRIMARY_ARM <- "arm1_magnitude"
DROP_DIR <- "C(gut -> env)"
DIR_LABELS <- c("C(env -> gut)"   = "environment in child stool",
                "C(gut -> dung)"  = "child stool in cattle dung",
                "C(gut -> soil)"  = "child stool in soil",
                "C(dung -> soil)" = "cattle dung in soil")

bg <- t26a %>%
  filter(arm == PRIMARY_ARM, direction != DROP_DIR) %>%
  mutate(direction = factor(DIR_LABELS[direction],
                            levels = rev(DIR_LABELS)),
         raw_delta = arg_obs - bg_obs_mean,
         raw_z     = raw_delta / bg_obs_sd)
stopifnot(nrow(bg) == 4, !any(is.na(bg$direction)))

mk_bg_panel <- function(arg_col, mean_col, sd_col, xlab, title, ann,
                        drop_y = FALSE, as_pct = FALSE, zero_rule = FALSE) {
  d <- bg %>% transmute(direction,
                        arg = .data[[arg_col]],
                        m   = .data[[mean_col]],
                        s   = .data[[sd_col]],
                        ann = ann)
  rng <- range(c(d$arg, d$m - d$s, d$m + d$s))
  pad <- diff(rng) * 0.18
  g <- ggplot(d, aes(y = direction))
  if (zero_rule)
    g <- g + geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.45)
  g <- g +
    geom_errorbar(aes(xmin = m - s, xmax = m + s, colour = "matched non-ARG background"),
                  orientation = "y", width = .18, linewidth = .55) +
    geom_point(aes(x = m, colour = "matched non-ARG background"), size = 2.1) +
    geom_point(aes(x = arg, colour = "clinical ARG variants"), size = 3.1) +
    # anchored at a fixed x: centred text is clipped by the panel edge
    geom_text(aes(x = rng[1] - pad, label = ann), hjust = 0, vjust = -1.9,
              size = 3.1, colour = "grey30") +
    scale_colour_manual(values = c("clinical ARG variants" = "#B5462F",
                                   "matched non-ARG background" = "grey25"),
                        name = NULL,
                        breaks = c("clinical ARG variants", "matched non-ARG background")) +
    scale_x_continuous(limits = c(rng[1] - pad, rng[2] + pad),
                       labels = if (as_pct) scales::percent_format(accuracy = 1)
                                else waiver()) +
    labs(title = title, x = xlab, y = NULL) +
    theme_sgma +
    theme(legend.position = "bottom",
          axis.text = element_text(size = 10.5),
          legend.text = element_text(size = 10.5),
          axis.title.x = element_text(size = 10.5, lineheight = 1.15,
                                      margin = margin(t = 8)))
  if (drop_y) g <- g + theme(axis.text.y = element_blank())
  g
}

# The annotation is the z of the ARG value against the background draws on each
# panel's own scale. Adjusted p-values are in the caption.
ann_ses <- sprintf("z = %+.2f", bg$z_vs_bg)
ann_raw <- sprintf("z = %+.2f", bg$raw_z)

p30_raw <- mk_bg_panel("arg_obs", "bg_obs_mean", "bg_obs_sd",
                       "% of the source compartment's features\ndetected in the target",
                       "Observed containment", ann_raw, as_pct = TRUE)
p30_ses <- mk_bg_panel("arg_SES", "bg_SES_mean", "bg_SES_sd",
                       "standard deviations from the expectation\nfor a set of that size",
                       "Distance from chance expectation", ann_ses, drop_y = TRUE,
                       zero_rule = TRUE)
p30_div <- ggplot() +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.4) +
  theme_void()

f30 <- (p30_raw | p30_div | p30_ses) +
  patchwork::plot_layout(widths = c(1, 0.02, 1), guides = "collect") &
  theme(legend.position = "bottom")
save_fig(f30, "F30_background_containment_summary.png", width = 10.5, height = 4.4)
cat("\nFigure 4 complete.\n")
