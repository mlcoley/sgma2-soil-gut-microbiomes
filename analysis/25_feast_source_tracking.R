#!/usr/bin/env Rscript
# FEAST source tracking in a household contrast.
# Is the own-household source proportion higher than an unrelated household's? 
# Outputs: T42a (per-sink proportions), T42b (contrast and test per stratum), T42c
# (rarefaction-depth sensitivity), T42d (sign-flip null draws), Figure S6.
# Overrides: R5_NDECOY, R5_NREP, R5_NPERM, R5_NCORE, R5_COVERAGE, R5_LAYER.
# Usage: Rscript analysis/25_feast_source_tracking.R   (from the repository root)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(parallel); library(FEAST)
})
source("analysis/00_load_features.R")
source("analysis/util_pair_features.R")
source("analysis/util_fig_feast.R")

PRIMARY_LAYER <- Sys.getenv("R5_LAYER", "ARG variant (clinical)")
COVERAGE <- as.integer(Sys.getenv("R5_COVERAGE", "848"))
N_DECOY  <- as.integer(Sys.getenv("R5_NDECOY",   "5"))
N_REP    <- as.integer(Sys.getenv("R5_NREP",     "3"))
N_PERM   <- as.integer(Sys.getenv("R5_NPERM",    "9999"))
N_CORE   <- as.integer(Sys.getenv("R5_NCORE",    "6"))
EM_ITER  <- 1000L

set.seed(1)
log_lines <- c(sprintf("FEAST source tracking - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
L   <- r8_layers()
A   <- L$layers[[PRIMARY_LAYER]]
key <- L$key

# assert the sample alignment and the integer counts FEAST requires
stopifnot(identical(rownames(A), key$sample_id))
if (any(A != floor(A))) stop("FEAST requires integer counts; layer is not integral")
if (min(rowSums(A)) < COVERAGE)
  stop(sprintf("COVERAGE %d exceeds the shallowest sample (%d). FEAST_rarefy would silently leave %d samples unrarefied.",
               COVERAGE, min(rowSums(A)), sum(rowSums(A) < COVERAGE)))

key$compartment <- as.character(key$compartment)
HH <- sort(unique(key$household))
logmsg("cohort %d samples, %d households; layer '%s' with %d features",
       nrow(A), length(HH), PRIMARY_LAYER, ncol(A))
logmsg("COVERAGE %d (fixed); %d decoy draws x %d rarefaction replicates; EM %d",
       COVERAGE, N_DECOY, N_REP, EM_ITER)

samples_of <- function(h, cp) key$sample_id[key$household == h & key$compartment == cp]

# ---------------------------------------------------------------------------
# Strata
# ---------------------------------------------------------------------------
# `n_arm` holds the two arms to equal sample counts
# `loo` marks the control where the sink comes from the same compartment as its source arm and is excluded
STRATA <- list(
  `gut to soil`  = list(sink = "gut",  src = "soil", n_arm = 3, loo = FALSE, role = "primary"),
  `gut to dung`  = list(sink = "gut",  src = "dung", n_arm = 1, loo = FALSE, role = "secondary"),
  `dung to soil` = list(sink = "dung", src = "soil", n_arm = 3, loo = FALSE, role = "control"),
  `soil to soil` = list(sink = "soil", src = "soil", n_arm = 2, loo = TRUE,  role = "control")
)
R5_STRAT_ORDER <- names(STRATA)

# ---------------------------------------------------------------------------
# One FEAST run: one sink, own arm vs decoy arm, both as sources
# ---------------------------------------------------------------------------
# Returns the summed proportion of each arm plus the unknown
feast_once <- function(sink_id, own_ids, decoy_ids, seed, M = A, coverage = COVERAGE) {
  ids <- c(sink_id, own_ids, decoy_ids)
  C   <- M[ids, , drop = FALSE]
  md  <- data.frame(
    Env        = c("sink", rep("own", length(own_ids)), rep("decoy", length(decoy_ids))),
    SourceSink = c("Sink", rep("Source", length(own_ids) + length(decoy_ids))),
    # different_sources_flag = 1 pairs sink to sources by matching `id`
    id         = rep(1L, length(ids)),
    row.names  = ids, stringsAsFactors = FALSE)

  set.seed(seed)
  of <- tempfile("feast_")
  res <- suppressWarnings(invisible(utils::capture.output(
    out <- FEAST(C = C, metadata = md, EM_iterations = EM_ITER, COVERAGE = coverage,
                 different_sources_flag = 1, dir_path = tempdir(), outfile = of))))
  unlink(paste0(of, "_source_contributions_matrix.txt"))

  p  <- unlist(out)
  nm <- names(p)
  own_p   <- sum(p[nm %in% paste0(own_ids,   "_own")])
  decoy_p <- sum(p[nm %in% paste0(decoy_ids, "_decoy")])
  unk     <- unname(p[nm == "Unknown"][1])
  # the three must close on 1
  if (!isTRUE(all.equal(own_p + decoy_p + unk, 1, tolerance = 1e-6)))
    stop(sprintf("proportions do not sum to 1 for sink %s (own %.4f decoy %.4f unknown %.4f)",
                 sink_id, own_p, decoy_p, unk))
  c(own = own_p, decoy = decoy_p, unknown = unk)
}

# ---------------------------------------------------------------------------
# One stratum
# ---------------------------------------------------------------------------
run_stratum <- function(sname, M = A, coverage = COVERAGE,
                        ndecoy = N_DECOY, nrep = N_REP) {
  s   <- STRATA[[sname]]
  eligible <- HH[sapply(HH, function(h) length(samples_of(h, s$src)) >= s$n_arm + as.integer(s$loo))]
  sinks <- do.call(rbind, lapply(HH, function(h) {
    sk <- samples_of(h, s$sink)
    if (!length(sk)) return(NULL)
    data.frame(household = h, sink_id = sk, stringsAsFactors = FALSE)
  }))

  rows <- mclapply(seq_len(nrow(sinks)), mc.cores = N_CORE, FUN = function(i) {
    h <- sinks$household[i]; sk <- sinks$sink_id[i]
    own_pool <- setdiff(samples_of(h, s$src), sk)          # sink drops out
    if (length(own_pool) < s$n_arm) return(NULL)
    # deterministic per sink
    set.seed(1000L + i)
    own_ids  <- sort(sample(own_pool, s$n_arm))
    decoys   <- sample(setdiff(eligible, h), min(ndecoy, length(setdiff(eligible, h))))
    do.call(rbind, lapply(seq_along(decoys), function(d) {
      dh   <- decoys[d]
      dpool <- samples_of(dh, s$src)
      dec_ids <- sort(sample(dpool, s$n_arm))              # arms held equal
      do.call(rbind, lapply(seq_len(nrep), function(r) {
        v <- feast_once(sk, own_ids, dec_ids, seed = 7919L * i + 31L * d + r,
                        M = M, coverage = coverage)
        data.frame(stratum = sname, role = s$role, household = h, sink_id = sk,
                   decoy_household = dh, rep = r,
                   own = v[["own"]], decoy = v[["decoy"]], unknown = v[["unknown"]],
                   contrast = v[["own"]] - v[["decoy"]], stringsAsFactors = FALSE)
      }))
    }))
  })
  bad <- vapply(rows, function(x) inherits(x, "try-error") || is.character(x), logical(1))
  if (any(bad)) stop("a worker failed: ", rows[[which(bad)[1]]])
  do.call(rbind, rows)
}

t0 <- Sys.time()
T42a <- do.call(rbind, lapply(R5_STRAT_ORDER, function(sn) {
  ts <- Sys.time()
  r  <- run_stratum(sn)
  logmsg("  %-13s %-9s %3d sinks x %2d decoys x %d reps = %4d runs  unknown %.3f  [%.0fs]",
         sn, STRATA[[sn]]$role, length(unique(r$sink_id)), length(unique(r$decoy_household)),
         N_REP, nrow(r), mean(r$unknown), as.numeric(difftime(Sys.time(), ts, units = "secs")))
  r
}))
logmsg("FEAST sweep: %d runs in %.1f min", nrow(T42a),
       as.numeric(difftime(Sys.time(), t0, units = "mins")))
save_tab(T42a, "T42a_feast_source_proportions.tsv")

# ---------------------------------------------------------------------------
# Household-level contrast and sign-flip test
# ---------------------------------------------------------------------------
# Averaged to household before the test
hh_contrast <- T42a %>%
  group_by(stratum, role, household) %>%
  summarise(contrast = mean(contrast), own = mean(own), decoy = mean(decoy),
            unknown = mean(unknown), n_sink = n_distinct(sink_id), .groups = "drop")


signflip <- function(x, nperm) {
  obs <- mean(x)
  nul <- replicate(nperm, mean(x * sample(c(-1, 1), length(x), replace = TRUE)))
  list(obs = obs, nul = nul,
       p = (1 + sum(abs(nul) >= abs(obs))) / (nperm + 1))
}

null_store <- list()
T42b <- do.call(rbind, lapply(R5_STRAT_ORDER, function(sn) {
  d <- hh_contrast[hh_contrast$stratum == sn, ]
  if (!nrow(d)) return(NULL)
  sf <- signflip(d$contrast, N_PERM)
  null_store[[sn]] <<- sf$nul
  data.frame(stratum = sn, role = STRATA[[sn]]$role, n_household = nrow(d),
             mean_own = mean(d$own), mean_decoy = mean(d$decoy),
             mean_unknown = mean(d$unknown),
             contrast = sf$obs, contrast_sd = sd(d$contrast),
             null_mean = mean(sf$nul), null_sd = sd(sf$nul),
             n_positive = sum(d$contrast > 0),
             p = sf$p, n_perm = N_PERM, stringsAsFactors = FALSE)
}))
# BH over primary family only
T42b$q <- NA_real_
pf <- T42b$role == "primary"
T42b$q[pf] <- p.adjust(T42b$p[pf], method = "BH")
save_tab(T42b, "T42b_source_contrast_by_stratum.tsv")
save_tab(as.data.frame(hh_contrast), "T42b2_source_contrast_by_household.tsv")

T42d <- do.call(rbind, lapply(names(null_store), function(sn)
  data.frame(stratum = sn, draw = seq_along(null_store[[sn]]),
             contrast = null_store[[sn]], stringsAsFactors = FALSE)))
save_tab(T42d, "T42d_signflip_null_draws.tsv")

logmsg("")
for (i in seq_len(nrow(T42b))) with(T42b[i, ], logmsg(
  "  %-13s %-9s contrast %+.5f (own %.4f vs decoy %.4f, unknown %.3f) %2d/%d hh positive  p=%.4f",
  stratum, role, contrast, mean_own, mean_decoy, mean_unknown, n_positive, n_household, p))

# ---------------------------------------------------------------------------
# Precondition: control
# ---------------------------------------------------------------------------
fires <- function(sn) { r <- T42b[T42b$stratum == sn, ]; nrow(r) == 1 && r$p < 0.05 && r$contrast > 0 }
g_ss <- fires("soil to soil"); g_ds <- fires("dung to soil")
logmsg("")
logmsg("== precondition ==")
logmsg("  soil to soil control fires : %s", g_ss)
logmsg("  dung to soil control fires : %s   <- the operative one", g_ds)
if (g_ss && g_ds) {
  logmsg("  MET. FEAST recovers a cross-compartment household signal, so a null")
  logmsg("  in the primary stratum is informative.")
} else if (g_ss) {
  logmsg("  PARTIALLY MET. Machinery works; the cross-compartment control did not fire, so the")
  logmsg("  primary result is uninterpretable and not evidence of absence.")
} else {
  logmsg("  NOT MET.")
}

# ---------------------------------------------------------------------------
# T42c: rarefaction-depth sensitivity
# depth varied downward (every sample still exceeds the target)
# ---------------------------------------------------------------------------
sens_rows <- list()
for (cov in c(500L, 700L)) {
  for (sn in c("gut to soil", "dung to soil")) {
    r  <- run_stratum(sn, coverage = cov, ndecoy = 3, nrep = 2)
    hc <- r %>% group_by(household) %>% summarise(contrast = mean(contrast),
                                                  unknown = mean(unknown), .groups = "drop")
    sf <- signflip(hc$contrast, 1999)
    sens_rows[[length(sens_rows) + 1]] <- data.frame(
      analysis = "rarefaction depth", setting = sprintf("COVERAGE = %d", cov),
      layer = PRIMARY_LAYER, stratum = sn, n_household = nrow(hc),
      mean_unknown = mean(hc$unknown), contrast = sf$obs, p = sf$p,
      note = "downward only; every sample still exceeds the target depth",
      stringsAsFactors = FALSE)
  }
}

T42c <- do.call(rbind, sens_rows)
save_tab(T42c, "T42c_sensitivity.tsv")
logmsg("")
for (i in seq_len(nrow(T42c))) with(T42c[i, ], logmsg(
  "  %-22s %-18s %-14s %-13s unknown %.3f  contrast %+.5f  p=%.4f",
  analysis, setting, substr(layer, 1, 14), stratum, mean_unknown, contrast, p))

writeLines(log_lines, file.path(LOG_DIR, "25_feast_source_tracking.log"))

# Figure S6
r5_figure_contrast(T42b, T42d)

cat("\nFEAST source tracking complete.\n")
