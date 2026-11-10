#!/usr/bin/env Rscript
# Dietary intake and child gut.
# n = 39 the minimum detectable correlation is |r| = 0.44
# null outcome excludes large effects only.
# Outputs: T45a (recall, child), T45b (FFQ, household), T45c (feasibility screen),
# T45d (compositional R2 against the floor).
# Household metadata is read from $SGMA2_HOUSEHOLD.
# Usage: Rscript analysis/17_diet.R   (from the repository root)

suppressPackageStartupMessages({
  library(dplyr); library(vegan)
})
source("analysis/00_load_features.R")
source("analysis/util_exposure_panel.R")

set.seed(1)
N_PERM <- as.integer(Sys.getenv("R6_NPERM", "999"))

log_lines <- c(sprintf("Dietary intake - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines)+1]] <<- m }

# ---------------------------------------------------------------------------
# Diet frame
# ---------------------------------------------------------------------------
key  <- load_key()
diet <- diet_panel()
logmsg("diet frame: %d households, verified against source files", nrow(diet))

# ---------------------------------------------------------------------------
# Feasibility screen
# ---------------------------------------------------------------------------
MIN_NONMISS  <- 30
MIN_MINORITY <- 8
MIN_LEVELS   <- 5

CANDIDATES <- list(
  A = c(child_dairy_any = "any dairy products yesterday (child)",
        child_fish_any  = "any fish yesterday (child)",
        child_dds       = "dietary diversity, food groups (child, 24h)"),
  B = c(ffq_asf_score     = "animal-source food score (household, habitual)",
        ffq_dds           = "dietary diversity, food groups (household, habitual)",
        ffq_dairy_regular = "dairy products at least 'Sometimes' (household)"))

screen_one <- function(nm) {
  x <- diet[[nm]]; x <- x[!is.na(x)]
  nn <- length(x); u <- length(unique(x))
  is_bin <- all(x %in% c(0, 1))
  minority <- if (is_bin) min(sum(x == 1), sum(x == 0)) else NA_integer_
  ok <- if (is_bin) (nn >= MIN_NONMISS && minority >= MIN_MINORITY)
        else        (nn >= MIN_NONMISS && u >= MIN_LEVELS)
  reason <- if (ok) "usable"
    else if (nn < MIN_NONMISS) sprintf("only %d non-missing (<%d)", nn, MIN_NONMISS)
    else if (is_bin) sprintf("minority class %d (<%d)", minority, MIN_MINORITY)
    else sprintf("only %d distinct values (<%d)", u, MIN_LEVELS)
  data.frame(variable = nm, type = if (is_bin) "binary" else "continuous",
             n_nonmissing = nn, n_distinct = u, minority_n = minority,
             usable = ok, reason = reason, stringsAsFactors = FALSE)
}

t45c <- do.call(rbind, lapply(names(CANDIDATES), function(fam)
  cbind(family = fam, do.call(rbind, lapply(names(CANDIDATES[[fam]]), screen_one)))))
t45c$label <- unlist(CANDIDATES)[t45c$variable]
save_tab(t45c, "T45c_diet_feasibility_screen.tsv")

logmsg("")
logmsg("== feasibility screen ==")
for (i in seq_len(nrow(t45c)))
  logmsg("  [%s] %-18s %-11s %s", t45c$family[i], t45c$variable[i], t45c$type[i], t45c$reason[i])

# there is no habitual dietary-diversity contrast in ffq
logmsg("")
logmsg("children reaching 5 food groups in the 24-hour recall: %d of 39", sum(diet$child_dds >= 5))

usable <- split(t45c$variable[t45c$usable], t45c$family[t45c$usable])
logmsg("")
for (fam in names(usable))
  logmsg("family %s usable predictors: %s", fam, paste(usable[[fam]], collapse = ", "))

# ---------------------------------------------------------------------------
# Responses, three scalars, stool only
# ---------------------------------------------------------------------------
arg <- gate_layer(get_layer("arg", key), load_gating())
ko  <- get_layer("ko_all", key)
gf  <- genome_factor(key)
clin <- clinical_args()

gut_ids <- key$sample_id[key$compartment == "gut"]
gut_ids <- intersect(gut_ids, rownames(arg$M))
gut_ids <- intersect(gut_ids, rownames(ko$M))
argg <- arg$M[gut_ids, , drop = FALSE]
kog  <- ko$M[gut_ids, , drop = FALSE]
argp <- per_genome(argg, gf)[gut_ids, , drop = FALSE]
cl   <- intersect(clin, colnames(argp))

resp <- data.frame(
  sample_id     = gut_ids,
  hh_id         = key$household[match(gut_ids, key$sample_id)],
  arg_abundance = rowSums(argp[, cl, drop = FALSE]),
  arg_richness  = rowSums(argg[, cl, drop = FALSE] > 0),
  ko_shannon    = vegan::diversity(kog, index = "shannon"),
  stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(resp$hh_id))    
logmsg("")
logmsg("responses: %d gut samples, one per household; normalization %s",
       nrow(resp), gf$method)

RESPONSES <- c(arg_abundance = "clinical ARG abundance (per genome)",
               arg_richness  = "clinical ARG richness",
               ko_shannon    = "KO Shannon diversity")

# ---------------------------------------------------------------------------
# Scalar arm: standardized OLS, beta in SD units
# ---------------------------------------------------------------------------
z <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)

scalar_test <- function(xv, rv) {
  d <- data.frame(y = z(resp[[rv]]), x = diet[[xv]][match(resp$hh_id, diet$hh_id)])
  d <- d[complete.cases(d), ]
  d$x <- z(d$x)                           
  fit <- lm(y ~ x, data = d)
  ci <- confint(fit)["x", ]
  data.frame(response = rv, variable = xv, n = nrow(d),
             beta = unname(coef(fit)["x"]), ci_lo = ci[1], ci_hi = ci[2],
             p = summary(fit)$coefficients["x", 4], stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Compositional arm
# ---------------------------------------------------------------------------
clr_m <- function(M) { L <- log(M + 1); sweep(L, 1, rowMeans(L), "-") }
gower <- function(d) {
  D2 <- as.matrix(d)^2; n <- nrow(D2)
  J <- diag(n) - matrix(1 / n, n, n)
  -0.5 * J %*% D2 %*% J
}
r2_1df <- function(G, x, trG) {
  xc <- x - mean(x); ss <- sum(xc^2)
  if (ss <= 0) return(NA_real_)
  as.numeric(crossprod(xc, G %*% xc)) / (ss * trG)
}

permanova_1df <- function(G, trG, x, nperm = N_PERM) {
  obs  <- r2_1df(G, x, trG)
  perm <- vapply(seq_len(nperm), function(i) r2_1df(G, x[sample.int(length(x))], trG),
                 numeric(1))
  list(R2 = obs, R2_floor = mean(perm, na.rm = TRUE), R2_floor_sd = sd(perm, na.rm = TRUE),
       R2_excess = obs - mean(perm, na.rm = TRUE),
       p = (sum(perm >= obs - 1e-12, na.rm = TRUE) + 1) / (nperm + 1))
}

comp_test <- function(xv, dname) {
  M <- kog
  # predictor is aligned to resp$hh_id, so M's row order must equal resp$sample_id
  stopifnot(identical(rownames(M), resp$sample_id))
  x <- diet[[xv]][match(resp$hh_id, diet$hh_id)]
  keep <- !is.na(x)
  Mk <- M[keep, , drop = FALSE]; Mk <- Mk[, colSums(Mk) > 0, drop = FALSE]
  d <- if (dname == "Bray-Curtis") vegdist(tss(Mk), method = "bray") else dist(clr_m(Mk))
  G <- gower(d); trG <- sum(diag(G))
  r <- permanova_1df(G, trG, x[keep])
  #dispersion
  xk <- x[keep]
  grp <- if (all(xk %in% c(0, 1))) factor(xk) else factor(xk > median(xk, na.rm = TRUE))
  disp_p <- if (nlevels(grp) == 2)
    vegan::permutest(vegan::betadisper(d, grp), permutations = N_PERM)$tab$`Pr(>F)`[1] else NA_real_
  data.frame(response = sprintf("gut composition (%s)", dname), variable = xv,
             n = sum(keep), R2 = r$R2, R2_floor = r$R2_floor, R2_excess = r$R2_excess,
             p = r$p, dispersion_p = disp_p, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Run both families
# ---------------------------------------------------------------------------
run_family <- function(fam) {
  xs <- usable[[fam]]
  if (!length(xs)) return(NULL)
  sc <- do.call(rbind, lapply(xs, function(xv)
    do.call(rbind, lapply(names(RESPONSES), function(rv) scalar_test(xv, rv)))))
  cp <- do.call(rbind, lapply(xs, function(xv) comp_test(xv, "Bray-Curtis")))
  sc$arm <- "scalar"; cp$arm <- "compositional"
  sc$R2 <- NA_real_; sc$R2_floor <- NA_real_; sc$R2_excess <- NA_real_; sc$dispersion_p <- NA_real_
  cp$beta <- NA_real_; cp$ci_lo <- NA_real_; cp$ci_hi <- NA_real_
  out <- rbind(sc[, c("arm","response","variable","n","beta","ci_lo","ci_hi",
                      "R2","R2_floor","R2_excess","dispersion_p","p")],
               cp[, c("arm","response","variable","n","beta","ci_lo","ci_hi",
                      "R2","R2_floor","R2_excess","dispersion_p","p")])
  # BH within family, across both arms
  out$q_BH <- p.adjust(out$p, method = "BH")
  out$family <- fam
  out[order(out$q_BH, out$p), ]
}

t45a <- run_family("A")
t45b <- run_family("B")
save_tab(t45a, "T45a_diet_recall_child.tsv")
save_tab(t45b, "T45b_diet_ffq_household.tsv")

for (fam in c("A", "B")) {
  tt <- if (fam == "A") t45a else t45b
  logmsg("")
  logmsg("== family %s: %d tests, BH within family; minimum BH-adjusted p = %.3f ==",
         fam, nrow(tt), min(tt$q_BH))
  logmsg("  reaching BH-adjusted p < 0.05: %d", sum(tt$q_BH < 0.05))
  for (i in seq_len(min(3, nrow(tt))))
    logmsg("   %-14s %-28s %-18s p=%.3f BH-adjusted p=%.3f", tt$arm[i], tt$response[i], tt$variable[i],
           tt$p[i], tt$q_BH[i])
}

# ---------------------------------------------------------------------------
# T45d: compositional
# ---------------------------------------------------------------------------
all_usable <- unlist(usable, use.names = FALSE)
t45d <- do.call(rbind, lapply(c("Bray-Curtis", "Aitchison"), function(dn)
  do.call(rbind, lapply(all_usable, function(xv) comp_test(xv, dn)))))
t45d$R2_closed_form_check <- 1 / (nrow(resp) - 1)
save_tab(t45d, "T45d_diet_compositional_vs_floor.tsv")

cf   <- 1 / (nrow(resp) - 1)
worst <- max(abs(t45d$R2_floor - cf))
logmsg("")
logmsg("R2 floor check: empirical floor vs closed form 1/(n-1) = %.5f, max |diff| = %.5f",
       cf, worst)
if (worst > 0.01)
  stop(sprintf("R2 floor estimator FAILED its known-answer check (max |diff| = %.5f)", worst))

logmsg("")
logmsg("== compositional arm: R2 against its floor ==")
for (i in seq_len(nrow(t45d)))
  logmsg("  %-28s %-18s R2=%.4f floor=%.4f excess=%+.4f p=%.3f disp_p=%.3f",
         t45d$response[i], t45d$variable[i], t45d$R2[i], t45d$R2_floor[i],
         t45d$R2_excess[i], t45d$p[i], t45d$dispersion_p[i])
n_below <- sum(t45d$R2_excess <= 0)
logmsg("")
logmsg("%d of %d compositional tests sit at or below their permutation floor", n_below, nrow(t45d))

writeLines(unlist(log_lines), file.path(LOG_DIR, "17_diet.log"))
cat("\nDietary intake complete.\n")
