#!/usr/bin/env Rscript
# Figure 3 (containment of each compartment's features in another compartment) builder
# for 04_nestedness.R. Returns a ggplot.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})
if (!exists("TAB_DIR")) source("analysis/00_load_features.R")

# Layers in figure order
LAYER_ORDER_OQ2 <- c("ARG (clinical)", "ARG (all)",
                     "ARG allele (clinical)", "ARG allele (all)",
                     "KO (function)")

if (!exists("theme_sgma"))
  theme_sgma <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

# Blue = ARG family, violet = ARG variant, neutral = KO.
F12A_RUNG_COLS <- c("ARG family" = "#2a78d6", "ARG variant" = "#4a3aa7",
                    "KO" = "#4A4A4A")

f12a_containment <- function(t8b, t8g) {
  # One panel per compartment pair, one column per resolution within it
  DIR_DISPLAY <- c("C(gut -> dung)"  = "Child stool in cattle dung",
                   "C(gut -> soil)"  = "Child stool in soil",
                   "C(dung -> soil)" = "Cattle dung in soil")
  DIR_PAIR    <- c("C(gut -> dung)" = "gut-dung", "C(gut -> soil)" = "gut-soil",
                   "C(dung -> soil)" = "dung-soil")

  fd <- t8b %>%
    filter(threshold == 1, direction %in% names(DIR_DISPLAY)) %>%
    mutate(layer = factor(layer, levels = LAYER_ORDER_OQ2),
           pair  = unname(DIR_PAIR[direction]),
           panel = factor(unname(DIR_DISPLAY[direction]), levels = unname(DIR_DISPLAY))) %>%
    arrange(panel, layer)

  # a_total / b_total are the source and target set sizes of the pair, from the
  # same equal-effort procedure as T8b (T8g)
  sz <- t8g %>%
    filter(threshold == 1) %>%
    mutate(layer = factor(layer, levels = LAYER_ORDER_OQ2)) %>%
    select(layer, pair, a_total, b_total)

  fd <- left_join(fd, sz, by = c("layer", "pair"))
  stopifnot(nrow(fd) == 3 * length(LAYER_ORDER_OQ2), !any(is.na(fd$a_total)),
            !any(is.na(fd$b_total)))

  fd <- fd %>% mutate(
    rung = ifelse(grepl("^KO", layer), "KO",
           ifelse(grepl("allele", layer), "ARG variant", "ARG family")),
    disp = layer_display(as.character(layer), wrap = FALSE),
    clin = grepl("\\(clinical\\)$", disp),
    xlab = sprintf("%s\n%s\n%s / %s",
                   ifelse(clin, "Clinical ARG",
                          ifelse(grepl("\\(", disp), sub(" *\\(.*", "", disp),
                                 sub(" .*", "", disp))),
                   ifelse(clin, sub("^ARG ([a-z]+) .*$", "\\1", disp),
                          ifelse(grepl("\\(", disp),
                                 paste0("(", gsub("^.*\\((.*)\\)$", "\\1", disp), ")"),
                                 sub("^[^ ]+ ", "", disp))),
                   formatC(round(a_total), format = "d", big.mark = ","),
                   formatC(round(b_total), format = "d", big.mark = ",")),
    xlab = factor(xlab, levels = unique(xlab)))

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
              size = 3.2, colour = "grey20", vjust = -0.6) +
    facet_wrap(~ panel, nrow = 1, scales = "free_x") +
    scale_shape_manual(values = c("size-matched expectation" = 21, "observed" = 19),
                       breaks = c("size-matched expectation", "observed"),
                       name = NULL) +
    scale_colour_manual(values = F12A_RUNG_COLS, guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       breaks = seq(0, 1, .25),
                       limits = c(0, 1), expand = expansion(mult = c(.02, .12))) +
    labs(x = NULL, y = "Source compartment features\ndetected in the target") +
    theme_sgma +
    theme(legend.position = "top",
          plot.title = element_blank(),
          strip.text = element_text(size = 10.5, face = "bold"),
          axis.text.x = element_text(size = 8, lineheight = 1.15),
          panel.grid.major.x = element_blank())
}
