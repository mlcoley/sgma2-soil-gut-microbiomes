#!/usr/bin/env Rscript
# Figure S6: own-household minus decoy-household FEAST source proportion against its
# sign-flip null. Reuses the theme of util_fig_pairs.R.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("analysis/util_fig_pairs.R")

r5_figure_contrast <- function(t42b, t42d) {
  ord <- intersect(c("gut to soil", "gut to dung", "dung to soil", "soil to soil"),
                   unique(t42b$stratum))

  R5_STRIP <- c("gut to soil"   = "child stool to soil",
                "gut to dung"   = "child stool to cattle dung",
                "dung to soil"  = "cattle dung to soil\n(operative control)",
                "soil to soil"  = "soil to soil\n(positive control)")
  lab <- setNames(R5_STRIP[ord], ord)

  nd  <- t42d  %>% filter(stratum %in% ord) %>% mutate(stratum = factor(stratum, levels = ord))
  obs <- t42b  %>% filter(stratum %in% ord) %>%
    mutate(stratum = factor(stratum, levels = ord),
           # flip the annotation left of its line when the line sits far right
           lab   = sprintf("contrast %+.4f\np = %.4f\n%d/%d households > 0",
                           contrast, p, n_positive, n_household),
           lab_h = ifelse(contrast > 0, 1.06, -0.06))

  p <- ggplot(nd, aes(x = contrast)) +
    geom_histogram(bins = 45, fill = "grey82", colour = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_vline(data = obs, aes(xintercept = contrast), colour = "#B5462F", linewidth = 0.85) +
    geom_text(data = obs, aes(x = contrast, label = lab, hjust = lab_h),
              y = Inf, vjust = 1.2, size = 3.1, lineheight = 0.95, colour = "grey25") +
    facet_wrap(~ stratum, nrow = 1, scales = "free",
               labeller = labeller(stratum = lab)) +
    labs(x = "own-household minus decoy-household source proportion",
         y = "sign-flip draws") +
    r8_theme() +
    theme(legend.position = "none",
          strip.text  = element_text(face = "bold", size = 11),
          axis.text   = element_text(size = 9.5),
          axis.title  = element_text(size = 11),
          axis.title.x = element_text(margin = margin(t = 7)))
  save_fig(p, "F47_feast_household_contrast.png", width = 12, height = 4.2)
  p
}
