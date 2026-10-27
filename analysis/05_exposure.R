#!/usr/bin/env Rscript
# Exposure and soil-property associations with the compartment responses.
#
# Every household keeps cattle, so exposure terms measure intensity among keepers,
# not presence. The unit is the household (39), so every model has an effective n of
# 38-39 and power is the binding constraint. Exposures are tested one at a time
# (a multivariable model is not supportable at n = 39), with BH adjustment within
# each family: primary, exploratory, and four household socioeconomic variables added
# after the primary results were seen (kept out of the primary correction).
#   * soil: lme4::lmer with a household random intercept (3 microsites per household)
#   * stool, dung: lm on 39 / 38 rows
#   * Gaussian models; clinical ARG richness has no zeros (asserted below).
# A feasibility screen (T7a) runs before any model: binary variables need a minority
# class of >= 8 households, continuous variables >= 5 distinct values, and both need
# >= 30 non-missing. Variables that fail, or that are judged not validly measured,
# are excluded with the reason recorded.
# Outputs: T7a (screen), T7b (exposure associations), T7c (soil properties), and the
# forest plots for T7b and T7c.
# Household metadata is read from $SGMA2_HOUSEHOLD; no coordinate is read or written.
# Usage: Rscript analysis/05_exposure.R   (from the repository root)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(vegan); library(lme4)
})
source("analysis/00_load_features.R")
source("analysis/util_fig_exposure.R")
source("analysis/util_exposure_panel.R")

set.seed(1)
GATING <- load_gating()

log_lines <- c(sprintf("Exposure associations - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

rd <- exposure_rd

key <- load_key()
logmsg("cohort: %d samples across %d households", nrow(key), length(unique(key$household)))

# ===========================================================================
# Responses: per sample, then aligned to households
# ===========================================================================
arg  <- gate_layer(get_layer("arg", key), GATING)
ko   <- get_layer("ko_all", key)
gf   <- genome_factor(key)
clin <- intersect(clinical_args(), colnames(arg$M))
arg_pg <- per_genome(arg$M, gf)

resp <- data.frame(
  sample_id   = rownames(arg$M),
  hh_id       = arg$key$household,
  compartment = arg$key$compartment,
  arg_abundance = rowSums(arg_pg[, clin, drop = FALSE]),
  arg_richness  = rowSums(arg$M[, clin, drop = FALSE] > MIN_COUNT),
  stringsAsFactors = FALSE
)
kom <- ko$M[resp$sample_id, , drop = FALSE]
resp$ko_shannon <- diversity(kom, index = "shannon")

stopifnot(min(resp$arg_richness) > 0)
logmsg("clinical ARG richness range %d-%d, zero counts: %d",
       min(resp$arg_richness), max(resp$arg_richness), sum(resp$arg_richness == 0))
logmsg("responses: clinical ARG abundance (%s), clinical ARG richness, KO Shannon", gf$method)

# ===========================================================================
# Exposures: household level
# ===========================================================================
husb <- rd("sgma2_hq_module4_hhanimalhusbandry.csv")
n_cows_hh <- sum(husb$hhlivestock_cows == 1, na.rm = TRUE)
logmsg("households keeping cattle: %d of %d", n_cows_hh, nrow(husb))

num <- function(x) suppressWarnings(as.numeric(x))
exp_df <- exposure_panel()

PRIMARY     <- EXPOSURE_PRIMARY
EXPLORATORY <- EXPOSURE_EXPLORATORY
POSTHOC_SES <- EXPOSURE_POSTHOC_SES

# ===========================================================================
# Covariate feasibility screen -> T7a
# ===========================================================================
MIN_NONMISS <- 30
MIN_MINORITY <- 8
MIN_LEVELS   <- 5

# Construct-validity exclusions: variables that pass the statistical screen but do
# not measure what the name says. Reported in T7a with the reason and excluded from
# modelling.
INVALID_MEASURE <- c(
  dung_to_crops = paste0(
    "not validly measured — the source item anim_confined_dungcrops_Cows sits ",
    "inside the cow-confinement module, which was SKIPPED for 14 households ",
    "(every field in that module is NA for the same 14). The NA->0 recode ",
    "therefore scores them as not applying dung to crops. Confinement and dung ",
    "application are unrelated in practice (field knowledge), so at the ",
    "asked-household rate (19/25 = 76%) roughly 11 of ",
    "the 14 are misclassified as unexposed — ~28% of the cohort, attenuating ",
    "toward the null. Restricting to the 25 who were asked gives 19/6, which ",
    "fails MIN_MINORITY, so no valid version of this variable exists here."))

fam_of <- exposure_family_of
screen_one <- function(v, nm) {
  x  <- v[!is.na(v)]
  nn <- length(x); u <- length(unique(x))
  is_bin <- all(x %in% c(0, 1))
  minority <- if (is_bin) min(sum(x == 1), sum(x == 0)) else NA_integer_
  ok <- if (is_bin) (nn >= MIN_NONMISS && minority >= MIN_MINORITY)
        else        (nn >= MIN_NONMISS && u >= MIN_LEVELS)
  reason <- if (ok) "usable"
    else if (nn < MIN_NONMISS) sprintf("only %d non-missing (<%d)", nn, MIN_NONMISS)
    else if (is_bin) sprintf("minority class %d (<%d) — cannot support a test", minority, MIN_MINORITY)
    else sprintf("only %d distinct values (<%d)", u, MIN_LEVELS)
  excl_type <- if (ok) NA_character_ else "statistical_screen"
  if (nm %in% names(INVALID_MEASURE)) {
    ok <- FALSE
    reason <- unname(INVALID_MEASURE[[nm]])
    excl_type <- "construct_validity"
  }
  data.frame(variable = nm,
             family = fam_of(nm),
             type = if (is_bin) "binary" else "continuous",
             n_nonmissing = nn, n_missing = length(v) - nn,
             n_distinct = u, minority_n = minority,
             usable = ok, exclusion_type = excl_type, reason = reason,
             stringsAsFactors = FALSE)
}
cand <- c(PRIMARY, EXPLORATORY, POSTHOC_SES)
t7a  <- do.call(rbind, lapply(cand, function(nm) screen_one(exp_df[[nm]], nm)))
t7a  <- t7a[order(t7a$family != "primary", !t7a$usable, t7a$variable), ]
save_tab(t7a, "T7a_covariate_feasibility_screen.tsv")

logmsg("== covariate feasibility screen ==")
for (i in seq_len(nrow(t7a)))
  logmsg("  [%-11s] %-22s %-10s n=%2d  %s%s", t7a$family[i], t7a$variable[i], t7a$type[i],
         t7a$n_nonmissing[i], if (t7a$usable[i]) "USABLE" else "EXCLUDED — ",
         if (t7a$usable[i]) "" else t7a$reason[i])
usable <- t7a$variable[t7a$usable]
prim_usable <- intersect(PRIMARY, usable)
logmsg("-> %d of %d candidates usable (%d of %d primary)",
       length(usable), nrow(t7a), length(prim_usable), length(PRIMARY))

# ===========================================================================
# Power: smallest detectable effect at 80% power, alpha 0.05 (two-sided)
# ===========================================================================
mde_r <- function(n, alpha = 0.05, power = 0.8) {
  z <- qnorm(1 - alpha / 2) + qnorm(power)
  tanh(z / sqrt(n - 3))                       # Fisher-z inversion
}
n_hh <- length(unique(key$household))
MDE_R <- mde_r(n_hh)
MDE_D <- (qnorm(0.975) + qnorm(0.8)) * sqrt(4 / n_hh)   # binary exposure, balanced split
logmsg("n = %d households: smallest detectable |r| = %.2f (r^2 = %.0f%% of variance); smallest binary contrast = %.2f SD",
       n_hh, MDE_R, 100 * MDE_R^2, MDE_D)

# ===========================================================================
# Exposure -> response models -> T7b
# ===========================================================================
# Responses are standardized, so coefficients are in SD units.
RESPONSES <- c(arg_abundance = "clinical ARG abundance (per genome)",
               arg_richness  = "clinical ARG richness",
               ko_shannon    = "KO Shannon diversity")

fit_one <- function(cmp, rsp, xv) {
  d <- resp[resp$compartment == cmp, ]
  d$x <- exp_df[[xv]][match(d$hh_id, exp_df$hh_id)]
  d$y <- d[[rsp]]
  d <- d[!is.na(d$x) & !is.na(d$y), ]
  n_h <- length(unique(d$hh_id))
  if (n_h < MIN_NONMISS) return(NULL)
  d$y <- as.numeric(scale(d$y))
  xs  <- if (length(unique(d$x)) > 2) as.numeric(scale(d$x)) else d$x
  d$x <- xs
  est <- tryCatch({
    if (cmp == "soil") {
      m  <- lmer(y ~ x + (1 | hh_id), data = d, REML = FALSE)
      if (isSingular(m)) SINGULAR_FITS <<- SINGULAR_FITS + 1L
      cf <- summary(m)$coefficients
      se <- cf["x", "Std. Error"]
      # lmerMod gives no p-value: normal approximation
      c(cf["x", "Estimate"], se, 2 * pnorm(-abs(cf["x", "Estimate"] / se)))
    } else {
      m  <- lm(y ~ x, data = d)
      cf <- summary(m)$coefficients
      c(cf["x", "Estimate"], cf["x", "Std. Error"], cf["x", "Pr(>|t|)"])
    }
  }, error = function(e) rep(NA_real_, 3))
  data.frame(compartment = cmp, response = rsp, exposure = xv,
             family = fam_of(xv),
             n_samples = nrow(d), n_households = n_h,
             estimate_SD = est[1], se = est[2], p = est[3],
             ci_lo = est[1] - 1.96 * est[2], ci_hi = est[1] + 1.96 * est[2],
             stringsAsFactors = FALSE)
}

logmsg("== exposure associations ==")
SINGULAR_FITS <- 0L
grid <- expand.grid(cmp = COMP_LEVELS, rsp = names(RESPONSES), xv = usable,
                    stringsAsFactors = FALSE)
t7b <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  fit_one(grid$cmp[i], grid$rsp[i], grid$xv[i])))

# BH within each family separately, never pooled.
t7b$q_BH <- ave(t7b$p, t7b$family, FUN = function(p) p.adjust(p, method = "BH"))
ip <- t7b$family == "primary"
t7b$mde_r_at_80pct <- round(MDE_R, 3)
t7b$contamination_filter <- gating_col(GATING)
t7b$normalization <- gf$method
t7b <- t7b[order(t7b$family != "primary", t7b$p), ]
save_tab(t7b, "T7b_exposure_associations.tsv")

logmsg("soil models with a singular fit (zero between-household variance): %d of %d",
       SINGULAR_FITS, sum(t7b$compartment == "soil"))
np <- sum(ip); nsig <- sum(t7b$q_BH < 0.05, na.rm = TRUE)
logmsg("primary family: %d tests (%d exposures x %d responses x %d compartments); reaching BH-adjusted p < 0.05: %d",
       np, length(prim_usable), length(RESPONSES), length(COMP_LEVELS), nsig)
top <- head(t7b[ip, ], 6)
for (i in seq_len(nrow(top)))
  logmsg("  %-5s %-14s ~ %-22s beta=%+.2f SD [%+.2f,%+.2f] p=%.3f BH-adjusted p=%.3f (n_hh=%d)",
         top$compartment[i], top$response[i], top$exposure[i], top$estimate_SD[i],
         top$ci_lo[i], top$ci_hi[i], top$p[i], top$q_BH[i], top$n_households[i])

save_fig(build_f11(t7b), F11_FILES[["f11"]], width = 11, height = 5)

# ===========================================================================
# Soil properties -> T7c
# ===========================================================================
# Soil chemistry values are ICRAF FTIR-model predictions, not wet chemistry. Only
# properties whose predictions are reliable against the variation in these soils
# (1 - (RMSE / SD)^2: pH 87%, SOC 88%, K 95%, Ca 91%, B 87%, Mg 78%, CEC 77%, TN 69%)
# are tested; Mn (20%), Cu and Al (13%), and Fe, PSI, S and exchangeable acidity
# (<= 0%) are excluded as unreliable.
SOIL_VARS <- c(soil_ph = "pH", soil_soc = "soil organic carbon", soil_tn = "total N",
               soil_cec = "CEC", soil_ext_calc = "extractable Ca",
               soil_ext_mg = "extractable Mg", soil_potassium = "K", soil_boron = "B")

ed <- rd("sgma2_soiledaphicproperties.csv")
soil <- resp[resp$compartment == "soil", ]
ed_m <- ed[match(soil$sample_id, ed$sample_id), ]
if (all(is.na(ed_m$sample_id))) {
  # the edaphic table is keyed by field_name (H###-S0#), not P###
  ed_m <- ed[match(key$field_name[match(soil$sample_id, key$sample_id)], ed$sample_id), ]
}
logmsg("== soil properties ==")
logmsg("per-microsite chemistry joined for %d of %d analyzed soils",
       sum(!is.na(ed_m$sample_id)), nrow(soil))

t7c <- do.call(rbind, lapply(names(SOIL_VARS), function(v) {
  d <- data.frame(hh_id = soil$hh_id, x = num(ed_m[[v]]))
  out <- lapply(names(RESPONSES), function(rsp) {
    dd <- cbind(d, y = soil[[rsp]])
    dd <- dd[!is.na(dd$x) & !is.na(dd$y), ]
    if (nrow(dd) < 30) return(NULL)
    dd$y <- as.numeric(scale(dd$y)); dd$x <- as.numeric(scale(dd$x))
    m  <- lmer(y ~ x + (1 | hh_id), data = dd, REML = FALSE)
    cf <- summary(m)$coefficients
    se <- cf["x", "Std. Error"]
    data.frame(driver = v, driver_label = unname(SOIL_VARS[v]), response = rsp,
               n_soils = nrow(dd), n_households = length(unique(dd$hh_id)),
               estimate_SD = cf["x", "Estimate"], se = se,
               p = 2 * pnorm(-abs(cf["x", "Estimate"] / se)), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}))
t7c$q_BH <- p.adjust(t7c$p, method = "BH")
t7c$ci_lo <- t7c$estimate_SD - 1.96 * t7c$se
t7c$ci_hi <- t7c$estimate_SD + 1.96 * t7c$se
t7c$chemistry_source <- "ICRAF FTIR prediction (not wet chemistry)"
t7c$contamination_filter <- gating_col(GATING)
t7c <- t7c[order(t7c$p), ]
save_tab(t7c, "T7c_soil_drivers.tsv")

logmsg("soil properties reaching BH-adjusted p < 0.05: %d of %d", sum(t7c$q_BH < 0.05), nrow(t7c))
for (i in seq_len(min(6, nrow(t7c))))
  logmsg("  %-14s ~ %-22s beta=%+.2f SD p=%.4f BH-adjusted p=%.4f (n=%d soils, %d hh)",
         t7c$response[i], t7c$driver_label[i], t7c$estimate_SD[i],
         t7c$p[i], t7c$q_BH[i], t7c$n_soils[i], t7c$n_households[i])

save_fig(build_f11b(t7c), F11_FILES[["f11b"]], width = 11, height = 4.6)

writeLines(unlist(log_lines), file.path(LOG_DIR, "05_exposure.log"))
cat("\nExposure associations complete.\n")
