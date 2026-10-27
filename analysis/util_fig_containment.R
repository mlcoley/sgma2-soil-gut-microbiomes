#!/usr/bin/env Rscript
# Figure 3 (containment of child-stool features in the environment) builder for
# 04_nestedness.R. Computes no statistic: it takes T8b and T8f and returns a ggplot.
#
# One direction is drawn because C(A->B) / C(B->A) = |B| / |A| exactly, so the
# mirrored direction only restates set sizes (both are printed under each column).
# Open circle = size-matched expectation, filled circle = observed, segment = the
# shortfall; the SES above each column is the same quantity standardized.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})
if (!exists("TAB_DIR")) source("analysis/00_load_features.R")

# Layers in figure order; 04_nestedness.R asserts it builds them in this order.
LAYER_ORDER_OQ2 <- c("ARG (clinical)", "ARG (all)",
                     "ARG allele (clinical)", "ARG allele (all)",
                     "KO (function)")

if (!exists("theme_sgma"))
  theme_sgma <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

# Blue = ARG family, violet = ARG variant, neutral = KO.
F12A_RUNG_COLS <- c("ARG family" = "#2a78d6", "ARG variant" = "#4a3aa7",
                    "KO" = "#4A4A4A")

f12a_containment <- function(t8b, t8f) {
  fd <- t8b %>%
    filter(threshold == 1, direction == "C(gut -> env)") %>%
    mutate(layer = factor(layer, levels = LAYER_ORDER_OQ2)) %>%
    arrange(layer)

  sz <- t8f %>%
    filter(threshold == 1) %>%
    mutate(layer = factor(layer, levels = LAYER_ORDER_OQ2)) %>%
    select(layer, gut_total, env_total)

  fd <- left_join(fd, sz, by = "layer")
  stopifnot(nrow(fd) == length(LAYER_ORDER_OQ2), !any(is.na(fd$gut_total)))

  fd <- fd %>% mutate(
    rung = ifelse(grepl("^KO", layer), "KO",
           ifelse(grepl("allele", layer), "ARG variant", "ARG family")),
    disp = layer_display(as.character(layer), wrap = FALSE),
    # split on the parenthesis so the resolution word is kept; four lines because the
    # set sizes on one line are wider than the column
    xlab = sprintf("%s\n(%s)\nchild stool %s\nenvironment %s",
                   sub(" *\\(.*", "", disp),
                   gsub("^.*\\((.*)\\)$", "\\1", disp),
                   formatC(round(gut_total), format = "d", big.mark = ","),
                   formatC(round(env_total), format = "d", big.mark = ",")),
    xlab = factor(xlab, levels = xlab))

  ggplot(fd, aes(xlab)) +
    geom_segment(aes(xend = xlab, y = null_mean, yend = observed_mean, colour = rung),
                 linewidth = 0.9) +
    geom_errorbar(aes(ymin = pmax(observed_mean - observed_sd, 0),
                      ymax = observed_mean + observed_sd, colour = rung),
                  width = .12, linewidth = .5) +
    geom_point(aes(y = null_mean, colour = rung, shape = "size-matched expectation"),
               size = 3, fill = "white", stroke = 1.1) +
    geom_point(aes(y = observed_mean, colour = rung, shape = "observed"), size = 3) +
    geom_text(aes(y = 1.0, label = sprintf("SES %.2f", SES_mean)),
              size = 3.4, colour = "grey20", vjust = -0.6) +
    scale_shape_manual(values = c("size-matched expectation" = 21, "observed" = 19),
                       breaks = c("size-matched expectation", "observed"),
                       name = NULL) +
    scale_colour_manual(values = F12A_RUNG_COLS, guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       breaks = seq(0, 1, .25),
                       limits = c(0, 1), expand = expansion(mult = c(.02, .12))) +
    labs(x = NULL, y = "Child stool features detected in the environment") +
    theme_sgma +
    theme(legend.position = "top",
          plot.title = element_blank(),
          axis.text.x = element_text(size = 9, lineheight = 1.15),
          panel.grid.major.x = element_blank())
}
