#!/usr/bin/env Rscript
# Per-analyte reliability of the ICRAF spectral predictions, and a screen to
# select soil properties to be tested in 05_exposure.R.
# Two statistics are reported for each analyte:
#   ICRAF R2   the model's performance on a 30% hold-out of ICRAF's own global library
#              (criterion R2 >= 0.6)
#   reliability = 1 - (RMSE_val / SD)^2, the share of the variation present in these
#              152 soils that a prediction with that error can resolve; SD is computed
#              from the per-sample predictions here (the report rounds SD to one decimal)
# Analytes with reliability >= 0.69 were retained. Input: ICRAF report table in
# data/metadata/icraf_J62_report_table2_validation.tsv.
# Output: T7e.
# Usage: Rscript analysis/12_soil_analyte_reliability.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/00_load_features.R")

set.seed(1)

ICRAF_TAB <- file.path(DATA_DIR, "metadata", "icraf_J62_report_table2_validation.tsv")

log_lines <- c(sprintf("Soil analyte reliability - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

# 8 retained, 5 excluded
SOIL_VARS <- c(soil_ph = "pH", soil_soc = "soil organic carbon", soil_tn = "total N",
               soil_cec = "CEC", soil_ext_calc = "extractable Ca",
               soil_ext_mg = "extractable Mg", soil_potassium = "K", soil_boron = "B")
SOIL_EXCLUDED_VARS <- c(soil_ext_alum = "extractable Al", soil_iron = "Fe",
                        soil_ext_mn = "extractable Mn", soil_psi = "PSI",
                        soil_exch_acid = "exchangeable acidity")
RELIABILITY_MIN <- 0.69

stopifnot(file.exists(ICRAF_TAB))
icraf <- read.delim(ICRAF_TAB, stringsAsFactors = FALSE)

ed_path <- file.path(META_DIR, "sgma2_soiledaphicproperties.csv")
stopifnot(file.exists(ed_path))
ed <- read.csv(ed_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
               na.strings = c("NA", ""))
logmsg("ICRAF prediction file: %d soil samples", nrow(ed))

# The analysed cohort is a subset of the soils ICRAF ran. SD is taken over all of
# them (the variation the prediction had to resolve)
key <- load_key()
soil_ids <- key$field_name[key$compartment == "soil"]
in_cohort <- ed$sample_id %in% soil_ids
logmsg("of these, in the analysed cohort: %d", sum(in_cohort))

sd_of <- function(v, idx = TRUE) { x <- suppressWarnings(as.numeric(v[idx]))
                                   sd(x, na.rm = TRUE) }

rows <- lapply(seq_len(nrow(icraf)), function(i) {
  v <- icraf$our_variable[i]
  present <- !is.na(v) && nzchar(v) && v %in% names(ed)
  sd_all <- if (present) sd_of(ed[[v]]) else NA_real_
  sd_coh <- if (present) sd_of(ed[[v]], in_cohort) else NA_real_
  rel    <- if (present) 1 - (icraf$rmse_val[i] / sd_all)^2 else NA_real_
  screen <- if (!present) "not measured in this cohort"
            else if (v %in% names(SOIL_VARS)) "retained"
            else if (v %in% names(SOIL_EXCLUDED_VARS)) "excluded"
            else "not tested (soil texture)"
  data.frame(
    icraf_property   = icraf$icraf_property[i],
    our_variable     = ifelse(present, v, NA_character_),
    analyte          = ifelse(present,
                              unname(c(SOIL_VARS, SOIL_EXCLUDED_VARS)[v]), NA_character_),
    model            = icraf$model[i],
    icraf_r2_val     = icraf$r2_val[i],
    icraf_rmse_val   = icraf$rmse_val[i],
    icraf_bias_val   = icraf$bias_val[i],
    icraf_rpiq_val   = icraf$rpiq_val[i],
    sd_this_cohort   = sd_all,
    sd_analysed_only = sd_coh,
    reliability      = rel,
    screen_outcome   = screen,
    stringsAsFactors = FALSE)
})
t7d <- do.call(rbind, rows)
t7d <- t7d[order(-replace(t7d$reliability, is.na(t7d$reliability), -Inf)), ]
rownames(t7d) <- NULL

## ---- table should reproduce the screen applied in 05_exposure.R
ret <- t7d$our_variable[t7d$screen_outcome == "retained"]
exc <- t7d$our_variable[t7d$screen_outcome == "excluded"]
stopifnot(setequal(ret, names(SOIL_VARS)), setequal(exc, names(SOIL_EXCLUDED_VARS)))

# retained analytes clear threshold and the excluded fail it
r_ret <- t7d$reliability[t7d$screen_outcome == "retained"]
r_exc <- t7d$reliability[t7d$screen_outcome == "excluded"]
stopifnot(all(r_ret >= RELIABILITY_MIN), all(r_exc < RELIABILITY_MIN))
logmsg("screen reproduces: %d retained (all >= %.0f%%), %d excluded (all < %.0f%%)",
       length(r_ret), 100 * RELIABILITY_MIN, length(r_exc), 100 * RELIABILITY_MIN)
logmsg("lowest retained: %.1f%%   highest excluded: %.1f%%",
       100 * min(r_ret), 100 * max(r_exc))

save_tab(t7d, "T7e_soil_analyte_reliability.tsv")

logmsg("")
logmsg("  %-16s %-22s %7s %8s %8s  %s", "analyte", "ICRAF property", "R2", "RMSE", "ours", "screen")
for (i in seq_len(nrow(t7d)))
  logmsg("  %-16s %-22s %7.3f %8.3f %7s  %s",
         ifelse(is.na(t7d$analyte[i]), "—", t7d$analyte[i]),
         t7d$icraf_property[i], t7d$icraf_r2_val[i], t7d$icraf_rmse_val[i],
         ifelse(is.na(t7d$reliability[i]), "—", sprintf("%.0f%%", 100 * t7d$reliability[i])),
         t7d$screen_outcome[i])
logmsg("")

writeLines(unlist(log_lines), file.path(LOG_DIR, "12_soil_analyte_reliability.log"))
cat("\nSoil analyte reliability complete.\n")
