#!/usr/bin/env Rscript
# Child covariates (age, sex, anthropometry, stool consistency) against the stool
# responses.
# 05_exposure.R tests household exposures; this tests the child's own
# characteristics. Gut only: these are properties of the child, so crossing them with
# dung and soil would add tests that carry no hypothesis. Same feasibility screen as
# 05_exposure.R (binary: minority class >= 8 households; continuous: >= 5 distinct
# values; both: >= 30 non-missing); one covariate at a time against clinical ARG
# abundance, clinical ARG richness and KO Shannon diversity (standardized), BH within
# the family.
# Child age is the age in whole months at stool collection, from the date of birth on the
# screening form. The date of birth in the household-demographics module disagrees for
# some households (some by exactly one year, consistent with a transcription error), and
# neither source can be shown correct from the data, so the identical family is re-run
# with that source (T39c). The deposited metadata carries both ages (child_age_months in
# sgma2_screening.csv and in the demographics table), not dates of birth. The disagreement
# is counted in the log.
# Outputs: T39a (screen), T39b (associations), T39c (age-source sensitivity).
# Household metadata is read from $SGMA2_HOUSEHOLD.
# Usage: Rscript analysis/16_child_covariates.R   (from the repository root)

suppressPackageStartupMessages({
  library(dplyr); library(vegan)
})
source("analysis/00_load_features.R")
source("analysis/util_exposure_panel.R")

set.seed(1)
GATING <- load_gating()

log_lines <- c(sprintf("Child covariates - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

rd <- exposure_rd

key <- load_key()
logmsg("cohort: %d samples across %d households", nrow(key), length(unique(key$household)))

# ===========================================================================
# Responses: as in 05_exposure.R, restricted to stool
# ===========================================================================
arg  <- gate_layer(get_layer("arg", key), GATING)
ko   <- get_layer("ko_all", key)
gf   <- genome_factor(key)
clin <- intersect(clinical_args(), colnames(arg$M))
arg_pg <- per_genome(arg$M, gf)

resp <- data.frame(
  sample_id     = rownames(arg$M),
  hh_id         = arg$key$household,
  compartment   = arg$key$compartment,
  arg_abundance = rowSums(arg_pg[, clin, drop = FALSE]),
  arg_richness  = rowSums(arg$M[, clin, drop = FALSE] > MIN_COUNT),
  stringsAsFactors = FALSE
)
kom <- ko$M[resp$sample_id, , drop = FALSE]
resp$ko_shannon <- diversity(kom, index = "shannon")
resp <- resp[resp$compartment == "gut", ]

# Stool is one sample per household, so a fixed-effects fit is exact.
stopifnot(!any(duplicated(resp$hh_id)))
logmsg("gut responses: %d samples, %d households (1:1 confirmed)",
       nrow(resp), length(unique(resp$hh_id)))
logmsg("responses: clinical ARG abundance (%s), clinical ARG richness, KO Shannon", gf$method)

# ===========================================================================
# Covariates
# ===========================================================================
stool <- rd("sgma2_stoolmetadata.csv")
anth  <- rd("sgma2_childanth.csv")
scr   <- rd("sgma2_screening.csv")
dem   <- rd("sgma2_hq_module1_hhdem_agestructure.csv")

# --- age-source reconciliation, recorded before either version is used ------
age_scr <- setNames(scr$child_age_months, scr$hh_id)
age_dem <- setNames(dem$child_age_months, dem$hh_id)
stopifnot(!any(is.na(age_scr)), !any(is.na(age_dem)), !anyNA(stool$hh_id))

hh_all  <- intersect(stool$hh_id, intersect(names(age_scr), names(age_dem)))
age_scr <- age_scr[hh_all]
age_dem <- age_dem[hh_all]

diff_mo   <- age_scr - age_dem
disagrees <- diff_mo != 0
year_typo <- abs(diff_mo) == 12

logmsg("")
logmsg("== age source reconciliation (whole months at stool collection) ==")
logmsg("households where the two sources disagree: %d of %d", sum(disagrees), length(hh_all))
logmsg("  ... consistent with a one-year transcription error (12 months apart): %d", sum(year_typo))
logmsg("Screening is primary (the form collected DOB for eligibility); the demographics-module")
logmsg("version runs as a full sensitivity (T39c).")

# --- the covariate panel -----------------------------------------------------
cov_df <- data.frame(hh_id = hh_all, stringsAsFactors = FALSE)
cov_df$child_age_months <- age_scr[match(cov_df$hh_id, hh_all)]
cov_df$bsc_sample  <- stool$bsc_sample[match(cov_df$hh_id, stool$hh_id)]
cov_df$bsc_typical <- stool$bsc_typical[match(cov_df$hh_id, stool$hh_id)]
cov_df$child_weight_kg <- anth$child_weight_kg[match(cov_df$hh_id, anth$hh_id)]
cov_df$child_height_cm <- anth$child_height_cm[match(cov_df$hh_id, anth$hh_id)]
# binary, coded 1 = Female
cov_df$child_female <- ifelse(dem$child_gender[match(cov_df$hh_id, dem$hh_id)] == "Female", 1L, 0L)
# expected to fail the screen; constructed so they appear in T39a with a reason
cov_df$breastfeeding   <- ifelse(scr$currently_breastfeeding[match(cov_df$hh_id, scr$hh_id)] == "Yes", 1L, 0L)
cov_df$daytime_at_school <- ifelse(scr$child_daytimelocation[match(cov_df$hh_id, scr$hh_id)] == "At school", 1L, 0L)
cov_df$residence_gt_6mo <- ifelse(scr$hh_residencetime[match(cov_df$hh_id, scr$hh_id)] == "Greater than 6 months", 1L, 0L)

COVARS <- setdiff(names(cov_df), "hh_id")

# assert the id join resolved
n_unmatched <- sum(!resp$hh_id %in% cov_df$hh_id)
if (n_unmatched > 0) stop("covariate join failed for ", n_unmatched, " gut samples")
for (v in COVARS)
  if (all(is.na(cov_df[[v]]))) stop("covariate is entirely NA after join: ", v)
logmsg("covariate join resolved for all %d gut households", nrow(resp))

# Weight and height at 20-59 months partly restate age.
logmsg("")
logmsg("== anthropometry vs age ==")
for (v in c("child_weight_kg", "child_height_cm"))
  logmsg("  %-16s vs child_age_months: Pearson r = %+.2f", v,
         cor(cov_df[[v]], cov_df$child_age_months, use = "complete.obs"))
logmsg("bsc_sample vs bsc_typical: Spearman %+.3f, identical in %d of %d households",
       cor(cov_df$bsc_sample, cov_df$bsc_typical, method = "spearman"),
       sum(cov_df$bsc_sample == cov_df$bsc_typical), nrow(cov_df))

# ===========================================================================
# Feasibility screen (thresholds as in 05_exposure.R)
# ===========================================================================
MIN_NONMISS  <- 30
MIN_MINORITY <- 8
MIN_LEVELS   <- 5

screen_one <- function(nm) {
  x  <- cov_df[[nm]]
  nn <- sum(!is.na(x)); u <- length(unique(x[!is.na(x)]))
  is_bin <- u <= 2
  minority <- if (is_bin) min(table(x[!is.na(x)])) else NA_integer_
  ok <- if (is_bin) (nn >= MIN_NONMISS && minority >= MIN_MINORITY)
        else        (nn >= MIN_NONMISS && u >= MIN_LEVELS)
  reason <- if (ok) ""
    else if (nn < MIN_NONMISS) sprintf("only %d non-missing (<%d)", nn, MIN_NONMISS)
    else if (is_bin) sprintf("minority class %d (<%d) — cannot support a test", minority, MIN_MINORITY)
    else sprintf("only %d distinct values (<%d)", u, MIN_LEVELS)
  data.frame(variable = nm, type = if (is_bin) "binary" else "continuous",
             n_nonmissing = nn, n_distinct = u,
             minority_n = minority, usable = ok, reason = reason,
             stringsAsFactors = FALSE)
}
t39a <- do.call(rbind, lapply(COVARS, screen_one))
t39a <- t39a[order(!t39a$usable, t39a$variable), ]
save_tab(t39a, "T39a_child_covariate_feasibility_screen.tsv")

usable <- t39a$variable[t39a$usable]
logmsg("")
logmsg("== covariate feasibility screen ==")
for (i in seq_len(nrow(t39a)))
  logmsg("  %-18s %-10s n=%2d %-6s %s", t39a$variable[i], t39a$type[i], t39a$n_nonmissing[i],
         if (t39a$usable[i]) "USABLE" else "EXCL", t39a$reason[i])
logmsg("usable: %d of %d", length(usable), nrow(t39a))

# ===========================================================================
# Power
# ===========================================================================
mde_r <- function(n, alpha = 0.05, power = 0.8) {
  z <- qnorm(1 - alpha / 2) + qnorm(power)
  tanh(z / sqrt(n - 3))
}
n_hh  <- nrow(resp)
MDE_R <- mde_r(n_hh)
MDE_D <- (qnorm(0.975) + qnorm(0.8)) * sqrt(4 / n_hh)
logmsg("")
logmsg("== power ==")
logmsg("n = %d stool households: smallest detectable |r| = %.2f (r^2 = %.0f%% of variance); smallest binary contrast = %.2f SD",
       n_hh, MDE_R, 100 * MDE_R^2, MDE_D)

# ===========================================================================
# The tests
# ===========================================================================
RESPONSES <- c(arg_abundance = "clinical ARG abundance (per genome)",
               arg_richness  = "clinical ARG richness",
               ko_shannon    = "KO Shannon diversity")

fit_one <- function(rsp, xv, cov_source = cov_df, label = "primary") {
  d <- resp
  d$x <- cov_source[[xv]][match(d$hh_id, cov_source$hh_id)]
  d$y <- d[[rsp]]
  d <- d[!is.na(d$x) & !is.na(d$y), ]
  if (nrow(d) < MIN_NONMISS) return(NULL)
  d$y <- as.numeric(scale(d$y))
  # binary stays 0/1 (group contrast in SD units); continuous is standardized
  if (length(unique(d$x)) > 2) d$x <- as.numeric(scale(d$x))
  m  <- lm(y ~ x, data = d)
  cf <- summary(m)$coefficients
  data.frame(analysis = label, response = rsp, covariate = xv,
             n_households = nrow(d),
             estimate_SD = cf["x", "Estimate"], se = cf["x", "Std. Error"],
             p = cf["x", "Pr(>|t|)"],
             ci_lo = cf["x", "Estimate"] - 1.96 * cf["x", "Std. Error"],
             ci_hi = cf["x", "Estimate"] + 1.96 * cf["x", "Std. Error"],
             stringsAsFactors = FALSE)
}

grid <- expand.grid(rsp = names(RESPONSES), xv = usable, stringsAsFactors = FALSE)
t39b <- do.call(rbind, lapply(seq_len(nrow(grid)),
                              function(i) fit_one(grid$rsp[i], grid$xv[i])))
t39b$q_BH <- p.adjust(t39b$p, method = "BH")
t39b$mde_r_at_80pct      <- round(MDE_R, 3)
t39b$contamination_filter <- gating_col(GATING)
t39b$normalization       <- gf$method
t39b$dob_source          <- "screening"
t39b <- t39b[order(t39b$p), ]
save_tab(t39b, "T39b_child_covariate_associations.tsv")

logmsg("")
logmsg("== associations (stool only) ==")
logmsg("family: %d tests (%d covariates x %d responses x 1 compartment), BH within",
       nrow(t39b), length(usable), length(RESPONSES))
for (i in seq_len(min(6, nrow(t39b))))
  logmsg("  %-14s ~ %-18s beta=%+.2f SD [%+.2f,%+.2f] p=%.3f BH-adjusted p=%.3f",
         t39b$response[i], t39b$covariate[i], t39b$estimate_SD[i],
         t39b$ci_lo[i], t39b$ci_hi[i], t39b$p[i], t39b$q_BH[i])
nsig <- sum(t39b$q_BH < 0.05, na.rm = TRUE)
logmsg("reaching BH-adjusted p < 0.05: %d of %d", nsig, nrow(t39b))

# ===========================================================================
# Date-of-birth source sensitivity: the identical family under the other source
# ===========================================================================
cov_alt <- cov_df
cov_alt$child_age_months <- age_dem[match(cov_alt$hh_id, hh_all)]
t39c <- do.call(rbind, lapply(seq_len(nrow(grid)),
                              function(i) fit_one(grid$rsp[i], grid$xv[i],
                                                  cov_source = cov_alt, label = "dob_dem_module")))
t39c$q_BH <- p.adjust(t39c$p, method = "BH")
t39c$dob_source <- "dem_module"
t39c <- t39c[order(t39c$p), ]
save_tab(t39c, "T39c_child_covariate_dob_sensitivity.tsv")

age_p <- t39b[t39b$covariate == "child_age_months", c("response", "estimate_SD", "p", "q_BH")]
age_a <- t39c[t39c$covariate == "child_age_months", c("response", "estimate_SD", "p", "q_BH")]
cmp   <- merge(age_p, age_a, by = "response", suffixes = c("_screening", "_dem"))
logmsg("")
logmsg("== age source sensitivity ==")
for (i in seq_len(nrow(cmp)))
  logmsg("  %-14s screening beta=%+.2f BH-adjusted p=%.3f  |  dem beta=%+.2f BH-adjusted p=%.3f", cmp$response[i],
         cmp$estimate_SD_screening[i], cmp$q_BH_screening[i],
         cmp$estimate_SD_dem[i], cmp$q_BH_dem[i])
logmsg("verdict identical under both sources: %s",
       if (all((cmp$q_BH_screening < 0.05) == (cmp$q_BH_dem < 0.05))) "yes" else "no")

writeLines(unlist(log_lines), file.path(LOG_DIR, "16_child_covariates.log"))
cat("\nChild covariates complete.\n")
