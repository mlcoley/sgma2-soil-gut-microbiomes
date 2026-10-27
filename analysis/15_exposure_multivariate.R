#!/usr/bin/env Rscript
# Exposure associations against compositional (multivariate) responses.
# Every scalar-response test regresses exposures on a summary (clinical ARG abundance,
# richness, KO Shannon); composition can shift without moving any of them, so each
# usable exposure (T7a) is tested against the Bray-Curtis (primary) and Aitchison
# (sensitivity) composition of each compartment's KO and ARG-family layers.
# Output: T31a. BH is applied within family (primary, exploratory, post hoc SES), on
# Bray-Curtis only, as in 05_exposure.R; dung_to_crops is excluded there.
# Methods:
#  * Soil has 117 samples from 39 households and every exposure is household-level,
#    so permutation shuffles the 39 household values and re-expands them. Stool and
#    dung have one sample per household and permute freely (the same operation with a
#    group size of 1).
#  * The R2 floor is the mean R2 across the permutation distribution, correct under
#    restricted permutation; for stool and dung the closed form 1/(n-1) is reported
#    beside it as a check.
#  * R2 of a 1-df term is computed in closed form from the Gower matrix, which makes
#    999 vector permutations cheap; it is validated against vegan::adonis2 below and
#    the script stops if they disagree.
#  * betadisper is run on every binary exposure and reported NA for continuous ones.
# Requires T7a (run 05_exposure.R first).
# Usage: Rscript analysis/15_exposure_multivariate.R   (from the repository root)

suppressPackageStartupMessages({ library(dplyr); library(vegan) })
source("analysis/00_load_features.R")
source("analysis/util_exposure_panel.R")

set.seed(1)
GATING <- load_gating()
N_PERM <- 999

log_lines <- c(sprintf("Exposure vs multivariate responses - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

# ===========================================================================
# Inputs
# ===========================================================================
key <- load_key()
arg <- gate_layer(get_layer("arg", key), GATING)
ko  <- get_layer("ko_all", key)

exp_df <- exposure_panel()

t7a_path <- file.path(TAB_DIR, "T7a_covariate_feasibility_screen.tsv")
if (!file.exists(t7a_path))
  stop("T7a not found; run 05_exposure.R first (it defines the usable variables)")
t7a <- read.delim(t7a_path, stringsAsFactors = FALSE)
USABLE <- t7a$variable[t7a$usable]
FAMILY <- setNames(t7a$family, t7a$variable)

logmsg("exposure panel: %d usable of %d screened (%s)", length(USABLE), nrow(t7a),
       paste(sprintf("%d %s", table(FAMILY[USABLE]), names(table(FAMILY[USABLE]))), collapse = ", "))
for (i in which(!t7a$usable))
  logmsg("  EXCLUDED [%s] %-22s %s", t7a$exclusion_type[i], t7a$variable[i],
         substr(t7a$reason[i], 1, 90))
stopifnot(!"dung_to_crops" %in% USABLE)

EXPOSURE_LABELS <- c(
  cows_head              = "cattle head count",
  total_livestock_head   = "total livestock head count",
  cows_in_compound       = "cattle penned/stabled in compound",
  child_touch_feces      = "child touches livestock faeces",
  child_antibiotics_2mo  = "child antibiotics, past 2 mo",
  vet_meds_cows          = "veterinary medications given to cattle",
  chickens_head          = "chicken head count",
  child_eats_soil        = "child eats soil",
  hh_leftover_abx        = "household uses leftover antibiotics",
  soap_available         = "soap available at handwashing station",
  water_untreated        = "drinking water untreated",
  ses_income_crops       = "income source: crops",
  ses_income_ownbusiness = "income source: own business",
  ses_kitchen_garden     = "kitchen garden in compound",
  ses_land_over_1acre    = "owns over 1 acre")
stopifnot(all(USABLE %in% names(EXPOSURE_LABELS)))

# ===========================================================================
# PERMANOVA machinery
# ===========================================================================
clr <- function(M) { L <- log(M + 1); sweep(L, 1, rowMeans(L), "-") }

#' Gower double-centred matrix of a distance object.
gower <- function(d) {
  D2 <- as.matrix(d)^2
  n  <- nrow(D2)
  J  <- diag(n) - matrix(1 / n, n, n)
  -0.5 * J %*% D2 %*% J
}

#' R2 of a 1-df term, in closed form. Identical to adonis2's R2; asserted below.
r2_1df <- function(G, x, trG) {
  xc <- x - mean(x)
  ss <- sum(xc^2)
  if (ss <= 0) return(NA_real_)
  as.numeric(crossprod(xc, G %*% xc)) / (ss * trG)
}

#' PERMANOVA for one (layer x exposure), with the permutation restricted to
#' households. `hh` gives each row's household; shuffling operates on the
#' UNIQUE households and re-expands, so microsites of one household always move
#' together. With one sample per household this reduces to free permutation.
permanova_hh <- function(G, trG, x, hh, nperm = N_PERM) {
  obs <- r2_1df(G, x, trG)
  uh  <- unique(hh)
  # household-level value; asserted constant within household
  hv  <- vapply(uh, function(h) {
    v <- unique(x[hh == h]); stopifnot(length(v) == 1); v }, numeric(1))
  idx <- match(hh, uh)
  perm <- vapply(seq_len(nperm), function(i) r2_1df(G, hv[sample.int(length(uh))][idx], trG),
                 numeric(1))
  list(R2 = obs,
       R2_null = mean(perm, na.rm = TRUE),
       R2_null_sd = sd(perm, na.rm = TRUE),
       # +1 in numerator and denominator: the observed value is one of the
       # possible arrangements, so p can never be 0.
       p = (sum(perm >= obs - 1e-12, na.rm = TRUE) + 1) / (nperm + 1),
       n_units = length(uh))
}

# Layers per compartment. Reference taxonomy is usable in stool only, and nothing else
# in the analysis works in taxonomic feature space, so no taxonomic layer is tested.
LAYERS_BY_COMP <- list(
  gut  = c("function (KO)" = "ko", "resistome (ARG)" = "arg"),
  dung = c("function (KO)" = "ko", "resistome (ARG)" = "arg"),
  soil = c("function (KO)" = "ko", "resistome (ARG)" = "arg"))
MAT <- list(ko = ko$M, arg = arg$M)

# ---------------------------------------------------------------------------
# Validation against vegan::adonis2 before any result is computed, on a real layer
# with one binary and one continuous predictor.
# ---------------------------------------------------------------------------
logmsg("")
logmsg("== validating the closed-form R2 against vegan::adonis2 ==")
{
  vk <- key[key$compartment == "gut", ]
  Mv <- MAT$ko[vk$sample_id, , drop = FALSE]; Mv <- Mv[, colSums(Mv) > 0, drop = FALSE]
  dv <- vegdist(tss(Mv), method = "bray")
  Gv <- gower(dv); trGv <- sum(diag(Gv))
  for (tv in c("cows_in_compound", "cows_head")) {
    xv <- exp_df[[tv]][match(vk$household, exp_df$hh_id)]
    keepv <- !is.na(xv)
    Gs <- gower(as.dist(as.matrix(dv)[keepv, keepv])); trGs <- sum(diag(Gs))
    mine <- r2_1df(Gs, xv[keepv], trGs)
    theirs <- vegan::adonis2(as.dist(as.matrix(dv)[keepv, keepv]) ~ xv[keepv],
                             permutations = 99, by = "terms")$R2[1]
    logmsg("  %-18s closed form %.6f vs adonis2 %.6f  (diff %.2e)", tv, mine, theirs,
           abs(mine - theirs))
    if (abs(mine - theirs) > 1e-8)
      stop("closed-form R2 disagrees with adonis2 — do not trust the results below")
  }
  logmsg("  -> agreement exact; closed form is adonis2's R2")
}

# ===========================================================================
# PERMANOVA across the exposure panel -> T31a
# ===========================================================================
logmsg("")
logmsg("== exposure associations against compositional responses ==")

rows <- list()
for (cmp in COMP_LEVELS) {
  ck <- key[key$compartment == cmp, ]
  for (lname in names(LAYERS_BY_COMP[[cmp]])) {
    M <- MAT[[LAYERS_BY_COMP[[cmp]][[lname]]]]
    M <- M[intersect(rownames(M), ck$sample_id), , drop = FALSE]
    M <- M[, colSums(M) > 0, drop = FALSE]
    ids <- rownames(M); hh <- key$household[match(ids, key$sample_id)]
    for (dn in c("Bray-Curtis", "Aitchison")) {
      d_full <- if (dn == "Bray-Curtis") vegdist(tss(M), method = "bray") else dist(clr(M))
      Dm <- as.matrix(d_full)
      for (xv in USABLE) {
        x <- exp_df[[xv]][match(hh, exp_df$hh_id)]
        keep <- !is.na(x)
        if (length(unique(x[keep])) < 2) next
        G <- gower(as.dist(Dm[keep, keep])); trG <- sum(diag(G))
        r <- permanova_hh(G, trG, x[keep], hh[keep])
        is_bin <- all(x[keep] %in% c(0, 1))
        # dispersion companion (undefined for a continuous predictor)
        disp_p <- NA_real_
        if (is_bin) {
          bd <- try(vegan::permutest(
            vegan::betadisper(as.dist(Dm[keep, keep]), factor(x[keep])),
            permutations = N_PERM), silent = TRUE)
          if (!inherits(bd, "try-error")) disp_p <- bd$tab$`Pr(>F)`[1]
        }
        rows[[length(rows) + 1]] <- data.frame(
          compartment = cmp, layer = lname, distance = dn,
          exposure = xv, exposure_label = unname(EXPOSURE_LABELS[xv]),
          family = unname(FAMILY[xv]),
          type = if (is_bin) "binary" else "continuous",
          n_samples = sum(keep), n_households = r$n_units,
          n_exposed = if (is_bin) sum(x[keep] == 1) else NA_integer_,
          R2 = round(r$R2, 5),
          R2_null_expected = round(r$R2_null, 5),
          R2_null_sd = round(r$R2_null_sd, 5),
          R2_excess = round(r$R2 - r$R2_null, 5),
          # closed form is valid only under free permutation of independent units
          R2_null_closed_form = if (cmp == "soil") NA_real_ else round(1 / (r$n_units - 1), 5),
          permutation = if (cmp == "soil") "restricted to households" else "free (1 sample/household)",
          p = round(r$p, 4), dispersion_p = disp_p, stringsAsFactors = FALSE)
      }
    }
  }
  logmsg("  %s: %d layers x %d exposures x 2 distances", cmp,
         length(LAYERS_BY_COMP[[cmp]]), length(USABLE))
}
t31a <- do.call(rbind, rows)

# BH within family, on the primary distance only.
t31a$q_BH <- NA_real_
bc <- t31a$distance == "Bray-Curtis"
t31a$q_BH[bc] <- ave(t31a$p[bc], t31a$family[bc], FUN = function(p) p.adjust(p, method = "BH"))
t31a$contamination_filter <- gating_col(GATING)
t31a <- t31a[order(t31a$distance != "Bray-Curtis", t31a$family != "primary", t31a$p), ]
save_tab(t31a, "T31a_exposure_compositional_associations.tsv")

# recompute `bc` after the reorder
n_comp_layer <- sum(lengths(LAYERS_BY_COMP))
bc <- t31a$distance == "Bray-Curtis"
stopifnot(sum(bc) == nrow(t31a) / 2,
          all(table(t31a$family[bc]) ==
                c(exploratory = 5, posthoc_SES = 4, primary = 6) * n_comp_layer))

# closed form vs empirical floor where 1/(n-1) is valid
chk <- t31a[!is.na(t31a$R2_null_closed_form), ]
logmsg("")
logmsg("closed-form vs empirical R2 floor (gut/dung, where 1/(n-1) is valid):")
logmsg("  max |empirical - closed form| = %.5f over %d tests",
       max(abs(chk$R2_null_expected - chk$R2_null_closed_form)), nrow(chk))
soil_floor <- t31a[t31a$compartment == "soil", ]
logmsg("soil floor under RESTRICTED permutation: %.5f-%.5f (1/(n-1) would have been %.5f)",
       min(soil_floor$R2_null_expected), max(soil_floor$R2_null_expected), 1 / 38)

for (fam in c("primary", "exploratory", "posthoc_SES")) {
  sel <- bc & t31a$family == fam
  logmsg("")
  logmsg("[%s] %d Bray-Curtis tests; reaching BH-adjusted p < 0.05: %d", fam, sum(sel),
         sum(t31a$q_BH[sel] < 0.05, na.rm = TRUE))
  logmsg("  R2 at or below its floor: %d of %d", sum(t31a$R2_excess[sel] <= 0), sum(sel))
  top <- head(t31a[sel, ], 4)
  for (i in seq_len(nrow(top)))
    logmsg("  %-5s %-17s %-38s R2=%.4f (floor %.4f, excess %+.4f) p=%.3f BH-adjusted p=%.3f",
           top$compartment[i], top$layer[i], top$exposure_label[i], top$R2[i],
           top$R2_null_expected[i], top$R2_excess[i], top$p[i], top$q_BH[i])
}

nsig <- sum(t31a$q_BH[bc] < 0.05, na.rm = TRUE)
logmsg("")
logmsg("%d Bray-Curtis test(s) reach BH-adjusted p < 0.05; %d of %d primary tests have R2 at or below the floor",
       nsig, sum(t31a$R2_excess[bc & t31a$family == "primary"] <= 0), sum(bc & t31a$family == "primary"))

writeLines(unlist(log_lines), file.path(LOG_DIR, "15_exposure_multivariate.log"))
cat("\nExposure vs multivariate responses complete.\n")
