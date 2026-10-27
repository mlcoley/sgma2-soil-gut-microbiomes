#!/usr/bin/env Rscript
# Same-household pair classifier and SHAP attribution.
# Given two samples, can a model trained on other households tell whether they came from
# the same household? Fixed similarity statistics (Bray-Curtis) read AUC 0.500 on the
# dung-to-soil stratum, where a household signal is established with rare-feature
# sharing (07_rare_feature_sharing.R); a learner can reweight features and in principle
# see what a fixed statistic averages away.
# Two positive controls are scored beside the two primary tests (child stool to soil,
# child stool to cattle dung): soil to soil, and dung to soil, the operative one, since a
# signal is known to be present there that the fixed statistic misses. If the classifier
# recovers it, a null in the primary strata is informative; if it reads about 0.500 there,
# the primary nulls are uninterpretable and are reported that way.
# The deliverables are the AUC against its permutation null, the smallest AUC detectable
# at 80% power from the measured null width, and the attribution conditional on a control
# firing. XGBoost (max_depth 3, eta 0.1, 200 rounds, subsample 0.8, colsample_bytree 0.5;
# no tuning), five-fold cross-validation with folds assigned to households, mean AUC over
# repeats; the null permutes household labels within compartment and re-runs the whole
# pipeline, including fold construction. Three specifications: the primary (presence
# co-occurrence plus summary columns) and two sensitivities (presence only; minimum
# relative abundance plus summaries).
# Outputs: T40c (calibration of the fixed statistics), T40d (classifier), T40e (SHAP),
# Figures S7 and S8.
# Memory-bound: run with R8_NCORE=6 and nothing else in R alongside it. Environment
# overrides: R8_NPERM (199), R8_NREPEAT (5), R8_NTHREAD (1), R8_NCORE (12), R8_LAYER.
# Usage: Rscript analysis/22_pair_classifier.R   (from the repository root)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(vegan); library(xgboost); library(parallel)
})
source("analysis/00_load_features.R")
source("analysis/util_pair_features.R")
source("analysis/util_fig_pairs.R")

PRIMARY_LAYER <- Sys.getenv("R8_LAYER", "ARG variant (clinical)")
N_PERM    <- as.integer(Sys.getenv("R8_NPERM",    "199"))
N_REPEAT  <- as.integer(Sys.getenv("R8_NREPEAT",  "5"))
N_THREAD  <- as.integer(Sys.getenv("R8_NTHREAD",  "1"))
N_CORE    <- as.integer(Sys.getenv("R8_NCORE",    "12"))
N_FOLD    <- 5
PRIMARY_ENCODING <- "presence"

# ---------------------------------------------------------------------------
# Specifications
# ---------------------------------------------------------------------------
# The primary specification supplies the summary columns (see util_pair_features.R). The
# first sensitivity omits them; its contrast with the primary measures how much of the
# result comes from supplying the sums. The second swaps presence for abundance.
SPECS <- list(
  list(id = "primary",       encoding = "presence",  summaries = TRUE,
       label = "presence co-occurrence + summary statistics"),
  list(id = "sens_nosummary", encoding = "presence",  summaries = FALSE,
       label = "presence co-occurrence only (as originally scoped)"),
  list(id = "sens_abundance", encoding = "abundance", summaries = TRUE,
       label = "minimum relative abundance + summary statistics")
)

# Fixed hyperparameters; no tuning search.
XGB_PARAMS <- list(objective = "binary:logistic", max_depth = 3, eta = 0.1,
                   subsample = 0.8, colsample_bytree = 0.5, min_child_weight = 1,
                   nthread = N_THREAD, eval_metric = "auc")
XGB_ROUNDS <- 200

set.seed(1)
log_lines <- c(sprintf("Same-household pair classifier - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines) + 1]] <<- m }

# ===========================================================================
# Inputs
# ===========================================================================
LAY   <- r8_layers(include_ko = FALSE)
key   <- LAY$key
pairs <- r8_pair_index(key)

logmsg("cohort: %d samples, %d households", nrow(key), length(unique(key$household)))
logmsg("pairs: %d total, %d same-household (%.2f%%)",
       nrow(pairs), sum(pairs$same_hh), 100 * mean(pairs$same_hh))
logmsg("")
logmsg("stratum          total  same_hh  role")
for (s in names(R8_STRATUM_LABEL)) {
  sel <- pairs$stratum == s
  logmsg("  %-14s %5d  %5d    %s", R8_STRATUM_LABEL[s], sum(sel),
         sum(pairs$same_hh[sel]), R8_STRATUM_ROLE[s])
}
for (s in c("gut-gut", "dung-dung"))
  logmsg("  %-14s %5d  %5d    dropped (structurally zero positives)",
         s, sum(pairs$stratum == s), sum(pairs$same_hh[pairs$stratum == s]))
for (nm in names(LAY$layers))
  logmsg("layer %-24s %5d features", nm, ncol(LAY$layers[[nm]]))

STRATA <- names(R8_STRATUM_LABEL)

# ===========================================================================
# T40c: calibration of the fixed statistics (what each can detect, per stratum and layer)
# ===========================================================================
logmsg("")
logmsg("== T40c: calibration of fixed statistics (permutation null, B = %d) ==",
       min(N_PERM * 5, 999))

calib_stats <- list()
for (lname in names(LAY$layers)) {
  M <- LAY$layers[[lname]]
  P <- M > 0
  calib_stats[[paste0(lname, " | Bray-Curtis similarity")]] <-
    (1 - as.matrix(vegdist(tss(M), "bray")))[cbind(pairs$i, pairs$j)]
  calib_stats[[paste0(lname, " | shared richness")]] <-
    tcrossprod(P * 1)[cbind(pairs$i, pairs$j)]
  prev <- colMeans(P)
  for (cap in c(0.25, 0.10)) {
    Pr <- P[, prev <= cap, drop = FALSE]
    if (ncol(Pr) < 10) next
    inter <- tcrossprod(Pr * 1); sz <- rowSums(Pr)
    calib_stats[[sprintf("%s | rare co-presence (prev<=%.2f)", lname, cap)]] <-
      inter[cbind(pairs$i, pairs$j)]
    calib_stats[[sprintf("%s | rare Jaccard (prev<=%.2f)", lname, cap)]] <-
      (inter / pmax(outer(sz, sz, "+") - inter, 1))[cbind(pairs$i, pairs$j)]
  }
  # confound for the unnormalised co-presence count: feature-rich samples co-detect more
  # features with no sharing, so joint richness is scored on its own
  szr <- rowSums(P[, prev <= 0.10, drop = FALSE])
  calib_stats[[paste0(lname, " | joint rare richness (confound probe)")]] <-
    szr[pairs$i] + szr[pairs$j]
}

B_CAL <- min(N_PERM * 5, 999)
perm_labels <- replicate(B_CAL, {
  hh <- r8_permute_households(key); hh[pairs$i] == hh[pairs$j]
}, simplify = FALSE)

t40c <- do.call(rbind, lapply(names(calib_stats), function(nm) {
  do.call(rbind, lapply(STRATA, function(s) {
    sel <- pairs$stratum == s
    sc  <- calib_stats[[nm]][sel]
    obs <- r8_auc(sc, pairs$same_hh[sel])
    nul <- vapply(perm_labels, function(L) r8_auc(sc, L[sel]), numeric(1))
    data.frame(layer = sub(" \\|.*$", "", nm), statistic = sub("^.*\\| ", "", nm),
               stratum = R8_STRATUM_LABEL[s], role = R8_STRATUM_ROLE[s],
               n_pos = sum(pairs$same_hh[sel]), n_pairs = sum(sel),
               AUC = obs, null_mean = mean(nul), null_sd = sd(nul),
               z = (obs - mean(nul)) / sd(nul),
               p = (1 + sum(nul >= obs)) / (length(nul) + 1),
               stringsAsFactors = FALSE)
  }))
}))
t40c$q <- NA_real_
for (st in unique(t40c$statistic)) {
  s <- t40c$statistic == st
  t40c$q[s] <- p.adjust(t40c$p[s], "BH")
}
save_tab(t40c, "T40c_fixed_statistic_calibration.tsv")

bc <- t40c[t40c$layer == PRIMARY_LAYER & t40c$statistic == "Bray-Curtis similarity", ]
for (r in seq_len(nrow(bc)))
  logmsg("  Bray-Curtis %-13s AUC=%.3f null=%.3f+/-%.3f z=%+.2f p=%.3f",
         bc$stratum[r], bc$AUC[r], bc$null_mean[r], bc$null_sd[r], bc$z[r], bc$p[r])

# ===========================================================================
# The classifier
# ===========================================================================
# One repeat: assign households to folds; per fold, train on pairs whose both members are
# training households and score pairs whose both members are held out. Returns the pooled
# AUC over held-out pairs and optionally the per-fold SHAP contributions.
run_one_repeat <- function(X, y, pair_hh_i, pair_hh_j, households,
                           want_shap = FALSE) {
  folds <- r8_household_folds(households, N_FOLD)
  fi <- folds[pair_hh_i]; fj <- folds[pair_hh_j]
  score <- rep(NA_real_, length(y))
  shap_acc <- NULL; shap_n <- 0

  for (f in seq_len(N_FOLD)) {
    tr <- which(fi != f & fj != f)      # both members in training households
    te <- which(fi == f & fj == f)      # both members held out
    if (length(te) == 0 || sum(y[tr]) < 2) next
    npos <- sum(y[tr]); nneg <- length(tr) - npos
    par  <- XGB_PARAMS
    par$scale_pos_weight <- if (npos > 0) nneg / npos else 1
    dtr <- xgb.DMatrix(X[tr, , drop = FALSE], label = y[tr])
    m   <- xgb.train(par, dtr, nrounds = XGB_ROUNDS, verbose = 0)
    dte <- xgb.DMatrix(X[te, , drop = FALSE])
    score[te] <- predict(m, dte)
    if (want_shap) {
      sc <- predict(m, dte, predcontrib = TRUE)      # exact TreeSHAP
      sc <- sc[, seq_len(ncol(X)), drop = FALSE]     # drop the bias column
      a  <- colMeans(abs(sc))
      shap_acc <- if (is.null(shap_acc)) a else shap_acc + a
      shap_n   <- shap_n + 1
    }
  }
  keep <- !is.na(score)
  list(auc = r8_auc(score[keep], y[keep] == 1),
       n_scored = sum(keep), n_pos_scored = sum(y[keep] == 1),
       shap = if (want_shap && shap_n > 0) shap_acc / shap_n else NULL)
}

# Mean AUC over N_REPEAT repeats; the null uses the identical function and repeat count.
run_statistic <- function(X, y, hh_i, hh_j, households, n_repeat, want_shap = FALSE) {
  res <- lapply(seq_len(n_repeat), function(r)
    run_one_repeat(X, y, hh_i, hh_j, households, want_shap = want_shap))
  list(auc = mean(vapply(res, function(z) z$auc, numeric(1)), na.rm = TRUE),
       auc_sd = sd(vapply(res, function(z) z$auc, numeric(1))),
       n_scored = res[[1]]$n_scored, n_pos_scored = res[[1]]$n_pos_scored,
       shap = if (want_shap) do.call(rbind, lapply(res, function(z) z$shap)) else NULL)
}

# ---------------------------------------------------------------------------
# Full run
# ---------------------------------------------------------------------------
logmsg("")
logmsg("== T40d: classifier, %d repeats, %d permutations, %d specifications ==",
       N_REPEAT, N_PERM, length(SPECS))
logmsg("   layer '%s'; %d parallel workers, %d xgboost thread(s) each",
       PRIMARY_LAYER, N_CORE, N_THREAD)

rows <- list(); null_draws <- list(); shap_store <- list(); shap_nulls <- list()
for (spec in SPECS) {
  logmsg("")
  logmsg("  -- specification '%s': %s", spec$id, spec$label)
  for (s in STRATA) {
    p_s  <- pairs[pairs$stratum == s, ]
    X    <- r8_encode(LAY$layers[[PRIMARY_LAYER]], p_s, spec$encoding, spec$summaries)
    y    <- as.integer(p_s$same_hh)
    hh_i <- key$household[p_s$i]; hh_j <- key$household[p_s$j]
    want_shap <- spec$id == "primary"

    t0  <- Sys.time()
    obs <- run_statistic(X, y, hh_i, hh_j, key$household, N_REPEAT,
                         want_shap = want_shap)

    # null: permute household labels within compartment, rebuild the pair labels, and
    # re-run the whole pipeline including fold construction
    nb <- mclapply(seq_len(N_PERM), function(b) {
      hp <- r8_permute_households(key)
      r  <- run_statistic(X, as.integer(hp[p_s$i] == hp[p_s$j]),
                          hp[p_s$i], hp[p_s$j], hp, N_REPEAT,
                          want_shap = want_shap && b <= 50)
      list(auc = r$auc,
           shap = if (!is.null(r$shap)) colMeans(r$shap, na.rm = TRUE) else NULL)
    }, mc.cores = N_CORE, mc.set.seed = TRUE)

    failed <- vapply(nb, function(z) inherits(z, "try-error") || is.null(z$auc), logical(1))
    if (any(failed)) stop(sprintf("%d permutation workers failed for %s / %s",
                                  sum(failed), s, spec$id))
    nul <- vapply(nb, function(z) z$auc, numeric(1))
    el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

    z   <- (obs$auc - mean(nul, na.rm = TRUE)) / sd(nul, na.rm = TRUE)
    pv  <- (1 + sum(nul >= obs$auc, na.rm = TRUE)) / (sum(!is.na(nul)) + 1)
    # smallest AUC separable from the null at 80% power, two-sided alpha 0.05
    mde <- mean(nul, na.rm = TRUE) + 2.802 * sd(nul, na.rm = TRUE)

    tag <- paste(spec$id, s, sep = "|")
    rows[[tag]] <- data.frame(
      spec = spec$id, spec_label = spec$label,
      stratum = R8_STRATUM_LABEL[s], role = R8_STRATUM_ROLE[s],
      layer = PRIMARY_LAYER, encoding = spec$encoding, summaries = spec$summaries,
      n_pairs = nrow(p_s), n_features = ncol(X), n_pos = sum(y),
      n_scored = obs$n_scored, n_pos_scored = obs$n_pos_scored,
      AUC = obs$auc, AUC_sd_across_repeats = obs$auc_sd,
      null_mean = mean(nul, na.rm = TRUE), null_sd = sd(nul, na.rm = TRUE),
      z = z, p = pv, detectable_AUC_80pct = mde,
      n_perm = N_PERM, n_repeat = N_REPEAT, minutes = el,
      stringsAsFactors = FALSE)
    null_draws[[tag]] <- nul
    if (want_shap) {
      shap_store[[s]] <- obs$shap
      sn <- do.call(rbind, lapply(nb, function(z) z$shap))
      shap_nulls[[s]] <- sn
    }
    logmsg("     %-13s (%-7s) AUC=%.3f null=%.3f+/-%.4f z=%+.2f p=%.4f detect>=%.3f [%.1f min]",
           R8_STRATUM_LABEL[s], R8_STRATUM_ROLE[s], obs$auc, mean(nul, na.rm = TRUE),
           sd(nul, na.rm = TRUE), z, pv, mde, el)
  }
}
t40d <- do.call(rbind, rows)
# BH within the primary family of the primary specification (two tests: stool to soil, stool
# to dung). Controls are reported uncorrected and labelled; sensitivities as a block.
t40d$q <- NA_real_
prim <- t40d$role == "primary" & t40d$spec == "primary"
t40d$q[prim] <- p.adjust(t40d$p[prim], "BH")
save_tab(t40d, "T40d_pair_classifier.tsv")

# ---------------------------------------------------------------------------
# Precondition: does a control fire?
# ---------------------------------------------------------------------------
ctrl <- t40d[t40d$role == "control" & t40d$spec == "primary", ]
gate_dungsoil <- ctrl$p[ctrl$stratum == "dung to soil"] < 0.05
gate_soilsoil <- ctrl$p[ctrl$stratum == "soil to soil"] < 0.05
GATE_MET <- isTRUE(gate_dungsoil) && isTRUE(gate_soilsoil)
# GATE_MET governs whether a null in the primary strata is informative; it needs the
# dung-to-soil control. RUN_SHAP governs whether anything is worth attributing, which needs
# one stratum to have detected something. Every SHAP row carries its own stratum's
# detection status.
RUN_SHAP <- isTRUE(gate_soilsoil) || isTRUE(gate_dungsoil)
detected <- setNames(ctrl_all <- t40d$p[t40d$spec == "primary"] < 0.05,
                     t40d$stratum[t40d$spec == "primary"])
logmsg("")
logmsg("== precondition ==")
logmsg("  soil to soil control fires : %s", gate_soilsoil)
logmsg("  dung to soil control fires : %s   <- the operative one", gate_dungsoil)
if (GATE_MET) {
  logmsg("  MET. A null in the primary strata is informative: the model has been")
  logmsg("  shown to detect a cross-compartment signal that the fixed statistic misses.")
} else if (isTRUE(gate_soilsoil)) {
  logmsg("  PARTIALLY MET. Soil to soil recovered, but the cross-compartment control did not")
  logmsg("  fire, so a null in the primary strata is uninterpretable and not evidence of absence.")
} else {
  logmsg("  NOT MET. Neither control fired; nothing in the primary strata can be read.")
}

# ===========================================================================
# T40e: SHAP attribution, conditional on the precondition
# ===========================================================================
if (RUN_SHAP) {
  logmsg("")
  logmsg("== T40e: SHAP attribution ==")
  # feature names come from the primary encoding, which carries six summary columns
  # beside the per-variant indicators
  feat <- colnames(r8_encode(LAY$layers[[PRIMARY_LAYER]], pairs[1:2, ],
                             PRIMARY_ENCODING, summaries = TRUE))
  stopifnot(length(feat) == ncol(shap_store[[STRATA[1]]]))
  t40e <- do.call(rbind, lapply(STRATA, function(s) {
    S <- shap_store[[s]]; if (is.null(S)) return(NULL)
    colnames(S) <- feat
    rk <- t(apply(-S, 1, rank, ties.method = "average"))
    nullmat <- shap_nulls[[s]]
    data.frame(
      stratum = R8_STRATUM_LABEL[s], role = R8_STRATUM_ROLE[s],
      stratum_detected_signal = isTRUE(detected[[R8_STRATUM_LABEL[s]]]),
      feature = feat,
      mean_abs_shap = colMeans(S, na.rm = TRUE),
      rank_median = apply(rk, 2, median), rank_iqr = apply(rk, 2, IQR),
      null_shap_mean = colMeans(nullmat, na.rm = TRUE),
      null_shap_p95  = apply(nullmat, 2, quantile, 0.95, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
  t40e$exceeds_null <- t40e$mean_abs_shap > t40e$null_shap_p95
  # univariate cross-check with no model; the six summary columns have no per-variant
  # analogue and are left NA
  Pm <- LAY$layers[[PRIMARY_LAYER]] > 0
  vfeat <- colnames(Pm)
  t40e$univariate_z <- NA_real_
  for (s in STRATA) {
    p_s <- pairs[pairs$stratum == s, ]
    co  <- Pm[p_s$i, , drop = FALSE] & Pm[p_s$j, , drop = FALSE]
    obs <- colMeans(co[p_s$same_hh, , drop = FALSE]) -
           colMeans(co[!p_s$same_hh, , drop = FALSE])
    nbu <- vapply(seq_len(min(199, N_PERM)), function(b) {
      hp <- r8_permute_households(key); L <- hp[p_s$i] == hp[p_s$j]
      colMeans(co[L, , drop = FALSE]) - colMeans(co[!L, , drop = FALSE])
    }, numeric(ncol(co)))
    zz <- (obs - rowMeans(nbu, na.rm = TRUE)) /
          pmax(apply(nbu, 1, sd, na.rm = TRUE), 1e-12)
    m <- t40e$stratum == R8_STRATUM_LABEL[s]
    t40e$univariate_z[m] <- unname(zz[match(t40e$feature[m], vfeat)])
  }
  save_tab(t40e, "T40e_shap_attribution.tsv")
  for (s in unique(t40e$stratum)) {
    top <- t40e[t40e$stratum == s, ]
    top <- top[order(-top$mean_abs_shap), ][1:5, ]
    logmsg("  %-13s top5: %s", s, paste(sprintf("%s(rank IQR %.0f)", top$feature, top$rank_iqr),
                                        collapse = ", "))
  }
} else {
  t40e <- NULL
}

# ===========================================================================
# Figures
# ===========================================================================
r8_figure_classifier(t40d)
if (RUN_SHAP) r8_figure_shap(t40e)

writeLines(log_lines, file.path(LOG_DIR, "22_pair_classifier.log"))
cat("\nPair classifier complete.\n")
