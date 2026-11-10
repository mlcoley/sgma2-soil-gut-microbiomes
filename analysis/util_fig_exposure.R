#!/usr/bin/env Rscript
# Builders for the exposure-association figures (forest plot of the 54 primary tests
# and the soil physico-chemical tests)
# Requires 00_load_features.R to be sourced first.

suppressPackageStartupMessages({ library(ggplot2) })

f11_read <- function(f) {
  p <- file.path(TAB_DIR, f)
  if (!file.exists(p)) stop("missing source table: ", p, call. = FALSE)
  read.delim(p, stringsAsFactors = FALSE)
}

# The shaded band is the smallest detectable effect at 80% power
build_f11 <- function(t7b = f11_read("T7b_exposure_associations.tsv")) {
  d <- t7b[t7b$family == "primary", ]
  stopifnot(nrow(d) == 54, length(unique(d$mde_r_at_80pct)) == 1)
  mde  <- d$mde_r_at_80pct[1]
  n_hh <- max(d$n_households)

  d$response    <- factor(response_display(d$response),
                          levels = unname(RESPONSE_DISPLAY))
  d$compartment <- factor(compartment_display(d$compartment),
                          levels = unname(COMPARTMENT_DISPLAY))
  d$exposure    <- exposure_display(d$exposure)
  cols <- setNames(unname(COMP_COLS[names(COMPARTMENT_DISPLAY)]),
                   unname(COMPARTMENT_DISPLAY))

  ggplot(d, aes(estimate_SD, reorder(exposure, estimate_SD), colour = compartment)) +
    annotate("rect", xmin = -mde, xmax = mde, ymin = -Inf, ymax = Inf,
             alpha = .07, fill = "grey20") +
    geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                   position = position_dodge(width = .6), linewidth = .4) +
    geom_point(position = position_dodge(width = .6), size = 1.6) +
    facet_wrap(~response) +
    scale_colour_manual(values = cols, name = NULL) +
    labs(x = "Standardized effect (SD of response per unit or per level of exposure)",
         y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
}

# Soil physico-chemical tests.
build_f11b <- function(t7c = f11_read("T7c_soil_drivers.tsv")) {
  stopifnot(nrow(t7c) == 24)
  t7c$response  <- factor(response_display(t7c$response),
                          levels = unname(RESPONSE_DISPLAY))
  t7c$survives  <- ifelse(t7c$q_BH < 0.05, "< 0.05", "≥ 0.05")

  ggplot(t7c, aes(estimate_SD, reorder(driver_label, estimate_SD), colour = survives)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = .4) +
    geom_point(size = 1.8) +
    facet_wrap(~response) +
    scale_colour_manual(values = c(`< 0.05` = "#B5462F", `≥ 0.05` = "grey55"),
                        name = "BH-adjusted p", breaks = c("< 0.05", "≥ 0.05")) +
    labs(x = "Standardized effect (SD of response per SD of soil property)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
}

F11_FILES <- c(f11  = "F11_exposure_associations.png",
               f11b = "F11b_soil_properties.png")
