#!/usr/bin/env Rscript
# Supplementary tables as Word documents (Tables S1-S15) and CSV files (Data S1-S6).
#   Rscript analysis/27_manuscript_tables.R          # all items
#   Rscript analysis/27_manuscript_tables.R T3       # one item by id (see SPECS)
# Styling and the reader-facing vocabulary checks are in util_manuscript_tables.R. to be run after

suppressPackageStartupMessages({
  library(dplyr)
})

source("analysis/00_load_features.R")
source("analysis/util_manuscript_tables.R")

set.seed(1)

N_PERM_STD <- 999                              # permutation min = 1/(N_PERM_STD + 1)
P_FLOOR    <- 1 / (N_PERM_STD + 1)

N_PERM_FEAST     <- 9999                       # 25_feast_source_tracking.R
N_PERM_CLASSIFY  <- 199                        # 24_pair_classifier.R
PRIMARY_COVERAGE <- 848L                       # fixed FEAST rarefaction depth

read_tab <- function(f) {
  path <- file.path(TAB_DIR, f)
  if (!file.exists(path)) stop("missing source table: ", path, call. = FALSE)
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# T3 -- PERMANOVA

spec_T3 <- function() {
  d <- read_tab("T3_permanova.tsv")
  u <- read_tab("T37a_unweighted_permanova.tsv")

  stopifnot(all(c("layer", "distance", "term", "n", "groups", "R2",
                  "R2_null_expected", "R2_excess", "F", "p",
                  "note", "contamination_filter") %in% names(d)))
  stopifnot(all(c("layer", "presence_rule", "distance", "term", "n", "groups",
                  "R2", "R2_null_expected", "R2_excess", "F", "p") %in% names(u)),
            nrow(u) == 12)

  dup <- u$presence_rule == "count >= 1" & u$distance == "Bray-Curtis (abundance)"
  stopifnot(sum(dup) == 2)
  ref <- d[d$term == "compartment" & d$distance == "Bray-Curtis", ]
  chk <- merge(u[dup, c("layer", "R2")], ref[, c("layer", "R2")], by = "layer",
               suffixes = c("_unweighted", "_main"))
  stopifnot(nrow(chk) == 2,
            all(abs(chk$R2_unweighted - chk$R2_main) < 1e-4))
  u <- u[!dup, ]

  rng <- function(lay) range(u$R2[u$layer == lay & grepl("presence", u$distance)])
  stopifnot(all(abs(rng("KO (function)")   - c(0.6261, 0.6942)) < 5e-4),
            all(abs(rng("ARG (resistome)") - c(0.4792, 0.6463)) < 5e-4))

  main <- data.frame(
    `Feature layer`  = layer_display(d$layer, wrap = FALSE),
    `Presence rule`  = rule_display("count >= 1"),
    `Distance`       = d$distance,
    `Term`           = term_display(d$term),
    `n (groups)`     = sprintf("%s (%s)", fmt_int(d$n), fmt_int(d$groups)),
    `R2`             = fmt_num(d$R2, 3),
    `Expected R2`    = fmt_num(d$R2_null_expected, 3),
    `Excess R2`      = fmt_num(d$R2_excess, 3),
    `F`              = fmt_num(d$F, 2),
    `p`              = fmt_p(d$p, floor_at = P_FLOOR),
    check.names = FALSE, stringsAsFactors = FALSE)

  sens <- data.frame(
    `Feature layer`  = layer_display(u$layer, wrap = FALSE),
    `Presence rule`  = rule_display(u$presence_rule),
    `Distance`       = sub(" \\(abundance\\)$", "", u$distance),
    `Term`           = term_display(u$term),
    `n (groups)`     = sprintf("%s (%s)", fmt_int(u$n), fmt_int(u$groups)),
    `R2`             = fmt_num(u$R2, 3),
    `Expected R2`    = fmt_num(u$R2_null_expected, 3),
    `Excess R2`      = fmt_num(u$R2_excess, 3),
    `F`              = fmt_num(u$F, 2),
    `p`              = fmt_p(u$p, floor_at = P_FLOOR),
    check.names = FALSE, stringsAsFactors = FALSE)

  df <- rbind(main, sens)
  stopifnot(nrow(df) == 20)

  filt <- unique(d$contamination_filter[grepl("^ARG", d$layer)])
  notes <- unique(d$note)

  tbl_spec(
    title = paste0("Table S9. Variation in functional and resistance-gene composition ",
                   "explained by compartment, and by household among soil samples (PERMANOVA)."),
    df = df,
    caption = c(
      paste0("PERMANOVA (vegan::adonis2), 999 permutations. R2 is the proportion of variation in the ",
             "distance matrix attributed to the term. Expected R2 is the value a random grouping of the ",
             "same number and size produces by construction, and excess R2 is the difference between ",
             "them. F is the pseudo-F ratio."),
      paste0("Compartment contrasts child stool, cattle dung and soil across all 194 samples. Household ",
             "(soils only) runs on the 117 soil samples alone and contrasts the three microsites of one ",
             "homestead against soils of other homesteads; those microsites lie metres apart, so ",
             "household and physical distance are confounded in this test. Compartment (soil mean ",
             "composition) repeats the compartment test with each household's three soil profiles ",
             "averaged to one."),
      paste0("The lower ten rows test the compartment term on Jaccard and Sorensen distances, which use ",
             "detection alone, under both presence rules. Neither presence rule is primary in this ",
             "study."),
      paste0("Bray-Curtis distances are computed on relative abundance. Aitchison distances are ",
             "Euclidean distances on the centred log-ratio of counts + 1."),
      paste0("p = 0.001 is the lower bound attainable with 999 permutations. No test in this table falls ",
             "below it."),
      paste0("Resistance-gene features are ", filt, " against the negative controls; KO features are\n             not filtered.")))
}

# Tclin -- clinical ARG abundance per average genome, parts A and B

spec_Tclin <- function() {
  a <- read_tab("T_clinical_arg_abundance_summary.tsv")
  b <- read_tab("T_clinical_arg_abundance_tests.tsv")

  stopifnot(all(c("compartment", "n", "median", "q1", "q3", "min", "max",
                  "quantity", "denominator", "contamination_filter") %in% names(a)),
            all(c("contrast", "test", "statistic", "p", "p_adj_BH", "median_ratio",
                  "p_superiority", "quantity", "denominator",
                  "contamination_filter") %in% names(b)))


  const <- function(d) unique(paste(d$quantity, d$denominator, d$contamination_filter))
  stopifnot(length(const(a)) == 1, identical(const(a), const(b)))

  quantity <- unique(a$quantity)
  denom    <- unique(a$denominator)
  filt     <- unique(a$contamination_filter)

  a <- a[match(COMP_LEVELS, a$compartment), ]
  stopifnot(!any(is.na(a$compartment)))

  df_a <- data.frame(
    `Compartment` = compartment_display(a$compartment),
    `n`           = fmt_int(a$n),
    `Median`      = fmt_num(a$median, 3),
    `Q1`          = fmt_num(a$q1, 3),
    `Q3`          = fmt_num(a$q3, 3),
    `Minimum`     = fmt_num(a$min, 3),
    `Maximum`     = fmt_num(a$max, 3),
    check.names = FALSE, stringsAsFactors = FALSE)

  is_kw  <- grepl("Kruskal", b$test)
  st_txt <- ifelse(is_kw, fmt_num(b$statistic, 1), fmt_int(b$statistic))

  df_b <- data.frame(
    `Contrast`                    = contrast_display(b$contrast),
    `Test`                        = b$test,
    `Statistic`                   = st_txt,
    `p`                           = fmt_p_tail(b$p),
    `BH-adjusted p`               = fmt_p_tail(b$p_adj_BH),
    `Median ratio`                = fmt_num(b$median_ratio, 2),
    `Probability of superiority`  = fmt_num(b$p_superiority, 2),
    check.names = FALSE, stringsAsFactors = FALSE)

  part_a <- tbl_spec(
    title = "(A) Distribution by compartment.",
    df = df_a,
    spanner = c("", "", "", "Interquartile range", "Interquartile range", "Range", "Range"),
    caption = c(
      paste0("Q1 and Q3 are the first and third quartiles, so the interquartile range is the ",
             "interval between them. Minimum and maximum are the lowest and highest sample values ",
             "in the compartment."),
      paste0("Values use the ribosomal-protein panel denominator. The same medians under the ",
             "MicrobeCensus genome-equivalent denominator are in Table S5; the two give different ",
             "absolute levels and the same compartment ordering.")))

  part_b <- tbl_spec(
    title = "(B) Omnibus and pairwise tests.",
    df = df_b,
    caption = c(
      paste0("The omnibus test is across all three compartments. The three pairwise tests are two-sided ",
             "Wilcoxon rank-sum tests, adjusted together by the Benjamini-Hochberg procedure. All ",
             "p-values are asymptotic and have no lower bound; values below 0.0001 are given as an ",
             "inequality."),
      paste0("The median ratio is the first compartment's median divided by the second's. The ",
             "probability of superiority is the proportion of all sample pairs, one drawn from each ",
             "compartment, in which the first compartment's value is the higher. A value of 0.50 means ",
             "the two are indistinguishable and 1.00 means every sample of the first exceeds every ",
             "sample of the second."),
      paste0("Child stool and the two environmental compartments separate completely: all 39 child stool ",
             "samples exceed all 155 cattle dung and soil samples, in 6,045 of 6,045 pairs. This was ",
             "verified under the primary clinical-subset definition; the three alternative definitions ",
             "were not checked."),
      paste0("Cattle dung and soil overlap almost entirely, at a probability of superiority of 0.609. ",
             "The BH-adjusted p of 0.045 for that contrast describes two distributions that overlap ",
             "across most of their range.")))

  tbl_parts(
    title = paste0("Table S12. Clinically relevant resistance-gene abundance per average ",
                   "genome in child stool, cattle dung and soil."),
    parts = list(part_a, part_b),
    caption = c(
      paste0("Values are ", quantity, ". Abundance is summed over the clinically relevant resistance\n             determinants and divided by the ", denom, " in the same sample, so 1.00 is one copy in\n             every genome present."),
      paste0("Resistance-gene features are ", filt, " against the negative controls."),
      paste0("n is the analysed cohort: 39 child stool, 38 cattle dung and 117 soil samples, three soil ",
             "microsites per household.")))
}

# Talpha -- within-compartment spread of alpha diversity

spec_Talpha <- function() {
  loc <- read_tab("T_alpha_location_tests.tsv")
  d   <- read_tab("T_alpha_dispersion.tsv")

  stopifnot(nrow(loc) == 4, nrow(d) == 4,
            all(c("layer", "metric", "chi2", "p", "p_gut_dung", "p_gut_soil",
                  "p_dung_soil") %in% names(loc)),
            all(c("layer", "metric", "chi2", "p", "mad_gut", "mad_dung", "mad_soil",
                  "iqr_gut", "iqr_dung", "iqr_soil") %in% names(d)))

  key <- function(x) paste(x$layer, x$metric)
  stopifnot(identical(key(loc), key(d)))

  METRIC <- c(richness = "Observed richness", shannon = "Shannon")

  LAYER  <- c(KO = layer_display("KO"), ARG = layer_display("ARG"))

  df <- data.frame(
    `Feature layer`  = unname(LAYER[d$layer]),
    `Metric`         = unname(METRIC[d$metric]),
    `Location χ²`    = fmt_num(loc$chi2, 1),
    `Location p`     = fmt_p_tail(loc$p),
    `Child stool vs cattle dung` = fmt_p_tail(loc$p_gut_dung),
    `Child stool vs soil`        = fmt_p_tail(loc$p_gut_soil),
    `Cattle dung vs soil`        = fmt_p_tail(loc$p_dung_soil),
    `Spread χ²`      = fmt_num(d$chi2, 1),
    `Spread p`       = fmt_p_tail(d$p),
    `MAD child stool` = fmt_num(d$mad_gut, 3),
    `MAD cattle dung` = fmt_num(d$mad_dung, 3),
    `MAD soil`        = fmt_num(d$mad_soil, 3),
    check.names = FALSE, stringsAsFactors = FALSE)

  tbl_spec(
    title = paste0("Table S11. Compartment differences in functional and resistance-gene alpha ",
                   "diversity, in location and in spread."),
    df = df,
    spanner = c("", "", "Kruskal-Wallis", "Kruskal-Wallis",
                "Pairwise Wilcoxon, BH-adjusted p", "Pairwise Wilcoxon, BH-adjusted p",
                "Pairwise Wilcoxon, BH-adjusted p",
                "Fligner-Killeen", "Fligner-Killeen",
                "Median absolute deviation", "Median absolute deviation",
                "Median absolute deviation"),
    caption = c(
      paste0("Location was tested by Kruskal-Wallis across the three compartments, followed by ",
             "pairwise two-sided Wilcoxon rank-sum tests with Benjamini-Hochberg adjustment ",
             "across the three pairs. Spread was tested by median-centred Fligner-Killeen on the ",
             "same values."),
      paste0("A location test compares whole distributions, so a rejection where spread also ",
             "differs is a statement about the distributions together. Both tests are reported ",
             "for every layer and metric for that reason."),
      paste0("Median absolute deviation is given per compartment, in the units of the metric: ",
             "gene counts for observed richness and index units for Shannon diversity."),
      paste0("Child stool is the most variable compartment on every metric in both layers. At ",
             "ARG-family resolution, Shannon diversity separates the compartments more weakly ",
             "than richness does, and cattle dung and soil are indistinguishable."),
      paste0("The ARG-family layer is the extended AMRFinderPlus annotation, which includes metal, ",
             "biocide and virulence-associated features. KO features are not filtered against the ",
             "negative controls."),
      paste0("All p-values are asymptotic and have no lower bound; values below 0.0001 are given ",
             "as an inequality."),
      "Analysed cohort: 39 child stool, 38 cattle dung and 117 soil samples."))
}

# T8g -- the feature partition of each compartment pair at equal effort

spec_T8g <- function() {
  d <- read_tab("T8g_pairwise_partition.tsv")

  stopifnot(all(c("layer", "threshold", "pair", "a_total", "a_total_sd",
                  "b_total", "b_total_sd", "shared", "shared_sd", "a_specific",
                  "a_specific_sd", "b_specific", "b_specific_sd") %in% names(d)))


  stopifnot(all(abs(d$shared + d$a_specific - d$a_total) < 0.15),
            all(abs(d$shared + d$b_specific - d$b_total) < 0.15))

  PAIR_DISPLAY <- c("gut-dung"  = "Child stool and cattle dung",
                    "gut-soil"  = "Child stool and soil",
                    "dung-soil" = "Cattle dung and soil")
  stopifnot(all(d$pair %in% names(PAIR_DISPLAY)))


  d$layer <- factor(d$layer, levels = unique(read_tab("T8f_private_repertoire_partition.tsv")$layer))
  d$pair  <- factor(d$pair, levels = names(PAIR_DISPLAY))
  d <- d[order(d$layer, d$threshold, d$pair), ]
  stopifnot(nlevels(d$layer) == 5, nrow(d) == 30)

  THRESH_DISPLAY <- c("1" = "≥ 1 read", "5" = "≥ 5 reads")

  df <- data.frame(
    `Feature layer`           = layer_display(d$layer, wrap = FALSE),
    `Detection threshold`     = unname(THRESH_DISPLAY[as.character(d$threshold)]),
    `Compartment pair`        = unname(PAIR_DISPLAY[as.character(d$pair)]),
    `First compartment`       = fmt_pm(d$a_total, d$a_total_sd),
    `Second compartment`      = fmt_pm(d$b_total, d$b_total_sd),
    `Co-detected`             = fmt_pm(d$shared, d$shared_sd),
    `Specific to first`       = fmt_pm(d$a_specific, d$a_specific_sd),
    `Specific to second`      = fmt_pm(d$b_specific, d$b_specific_sd),
    check.names = FALSE, stringsAsFactors = FALSE)

  stopifnot(!any(grepl("private", names(df), ignore.case = TRUE)))

  tbl_spec(
    title = paste0("Table S13. Co-detected and compartment-specific features of each pair of compartments, ",
                   "at equal sampling effort."),
    df = df,
    caption = c(
      paste0("Each pair of compartments is compared at equal sampling effort of 38 samples per ",
             "compartment over 100 subsampling draws. The first and second compartment are the two named ",
             "in the pair, in the order given. Values are the mean across draws with the standard ",
             "deviation across draws."),
      paste0("Co-detected counts features detected in both compartments. The two specific columns count ",
             "features detected in one compartment and absent from the other. Co-detected plus the specific ",
             "count equals the size of each compartment's set."),
      "Detection thresholds of 1 and 5 reads are both shown.",
      paste0("Resistance-gene features are contamination-filtered against the negative controls; KO ",
             "features are not filtered.")),
    spanner = c("", "", "", "Features detected", "Features detected",
                "Partition", "Partition", "Partition"))
}


# T13 --  120 household-sharing tests, spreadsheet supplement

spec_T13 <- function() {
  d <- read_tab("T13_rare_feature_sharing.tsv")

  stopifnot(all(c("layer", "prev_cap", "n_features", "pct_features_kept",
                  "n_samples_emptied", "metric", "pair", "n_households",
                  "within", "null_mean", "z", "ratio", "mde_pct", "p", "q_BH",
                  "contamination_filter") %in% names(d)),
            nrow(d) == 120)


  CAP_DISPLAY <- c("1" = "\u2264 100% (unrestricted)", "0.5" = "\u2264 50%",
                   "0.25" = "\u2264 25%", "0.1" = "\u2264 10%")
  cap <- unname(CAP_DISPLAY[as.character(d$prev_cap)])
  stopifnot(!any(is.na(cap)))

  stopifnot(all(d$p >= P_FLOOR - 1e-12))

  df <- data.frame(
    `Feature layer`                            = layer_display(d$layer, wrap = FALSE),
    `Household prevalence cap`                 = cap,
    `Features retained`                        = as.integer(d$n_features),
    `Features retained (%)`                    = round(d$pct_features_kept, 1),
    `Samples with no feature detected`         = as.integer(d$n_samples_emptied),
    `Similarity metric`                        = d$metric,
    `Compartment pair`                         = pair_display(d$pair),
    `Households`                               = as.integer(d$n_households),
    `Within-household mean similarity`         = round(d$within, 4),
    `Permuted mean`                            = round(d$null_mean, 4),
    `z`                                        = round(d$z, 2),
    `Ratio to permuted mean`                   = round(d$ratio, 3),
    `Smallest effect reaching significance (%)`= round(d$mde_pct, 1),
    `p`                                        = round(d$p, 3),
    `BH-adjusted p`                            = round(d$q_BH, 4),
    `Contamination filter`                     = d$contamination_filter,
    check.names = FALSE, stringsAsFactors = FALSE)

  stopifnot(!any(grepl("\\bq\\b|q_BH", names(df))))

  defs <- c(
    `Feature layer` = paste0("Feature set the test was run on. ARG family layers aggregate ",
      "genes to their resistance-gene symbol; ARG variant layers use the 95%-identity ",
      "catalog gene. Clinical panels are restricted to acquired resistance determinants; ",
      "complete panels add metal or biocide tolerance and virulence-associated features."),
    `Household prevalence cap` = paste0("Upper limit on the proportion of households in ",
      "which a feature may be detected for it to be retained. The unrestricted rows apply ",
      "no limit and are the full feature set; the three capped levels progressively remove ",
      "the features present in the most households."),
    `Features retained` = "Number of features passing the cap.",
    `Features retained (%)` = "Features retained as a percentage of the unrestricted set.",
    `Samples with no feature detected` = paste0("Number of the 194 samples left with zero ",
      "detected features after the cap is applied. These samples contribute nothing to the ",
      "similarities, which is the mechanical reason the tightest caps carry little power."),
    `Similarity metric` = paste0("Jaccard uses presence and absence only; Bray-Curtis is ",
      "weighted by abundance."),
    `Compartment pair` = "The two compartments compared.",
    `Households` = "Households contributing a sample of both compartments.",
    `Within-household mean similarity` = paste0("Mean similarity across pairs of samples ",
      "drawn from the same household."),
    `Permuted mean` = paste0("Mean of the same statistic across permutations in which ",
      "household labels are shuffled, so that pairs are matched at random rather than by ",
      "household. This is the value the observation is tested against."),
    `z` = paste0("Within-household mean minus the permuted mean, in units of the standard ",
      "deviation of the permuted distribution."),
    `Ratio to permuted mean` = paste0("Within-household mean divided by the permuted mean. ",
      "A value of 1 means household membership makes no difference."),
    `Smallest effect reaching significance (%)` = paste0("1.645 x the standard deviation of ",
      "the permuted distribution, as a percentage of its mean: the smallest effect that ",
      "would reach significance in this test. It corresponds to roughly 50% power, so it is ",
      "smaller than an 80%-power detection threshold, which would need about 2.49 standard ",
      "deviations. It rises steeply at the tighter caps, where the feature set is small and ",
      "the permuted distribution is wide."),
    `p` = paste0("One-sided permutation p-value from 999 permutations of the household ",
      "labels."),
    `BH-adjusted p` = paste0("The p-value after Benjamini-Hochberg correction across all 120 ",
      "tests in this file."),
    `Contamination filter` = "Whether the feature set was filtered against the negative controls.")

  tbl_spec(
    title = paste0("Within-household similarity against a household-permuted null, ",
                   "across feature layers, prevalence caps, similarity metrics and ",
                   "compartment pairs: 120 tests."),
    df = df,
    definitions = defs,
    caption = c(
      paste0("Each test compares the mean similarity between samples of two compartments drawn from the ",
             "same household against the distribution obtained when household labels are permuted. A ",
             "positive z means samples from one household resemble each other more than the study ",
             "population at large."),
      "All 120 tests were corrected together by the Benjamini-Hochberg procedure.",
      paste0("p is a permutation p-value with a lower bound of 1/(999 + 1) = 0.001. Five of the 120 rows ",
             "sit on that bound, which is the resolution limit of the test."),
      paste0("The unrestricted rows apply no prevalence cap and cover the full feature set. The three ",
             "capped levels remove features detected in every household."),
      paste0("Analysed cohort: 39 child stool, 38 cattle dung and 117 soil samples, three soil ",
             "microsites per household. A pair contributes a household where both of its compartments ",
             "were sampled.")))
}

# T7b / T7c -- the univariate exposure tests - two supplementary data files

.t7_shared <- function() {
  list(
    eff = paste0("Standardized regression coefficient: change in the response, in ",
                 "standard deviations of the response, per one-unit change in the ",
                 "predictor for binary predictors, or per standard deviation for ",
                 "continuous ones."),
    response = paste0("Clinical ARG abundance is copies per average genome; clinical ARG ",
                      "richness is the number of clinical resistance determinants detected; ",
                      "KO Shannon diversity is the Shannon index of the functional gene ",
                      "profile."),
    se = "Standard error of the effect.",
    ci_lo = "Lower bound of the 95% confidence interval on the effect.",
    ci_hi = "Upper bound of the 95% confidence interval on the effect.",
    p = "Two-sided p-value for the effect.",
    note_univariate = paste0("Each row is one univariate model fitted to one response; ",
                             "predictors were tested one at a time rather than jointly. ",
                             "Multivariate compositional versions of the same question are ",
                             "reported separately."),
    note_cohort = paste0("Analysed cohort: 39 households, with 39 child stool, 38 cattle ",
                         "dung and 117 soil samples."))
}

spec_T7b <- function() {
  sh <- .t7_shared()
  b <- read_tab("T7b_exposure_associations.tsv")

  stopifnot(nrow(b) == 135,
            all(c("compartment", "response", "exposure", "family", "n_samples",
                  "n_households", "estimate_SD", "se", "p", "ci_lo", "ci_hi",
                  "q_BH", "mde_r_at_80pct") %in% names(b)))

  filt <- unique(b$contamination_filter)
  norm <- unique(b$normalization)
  stopifnot(length(filt) == 1, length(norm) == 1)

  stopifnot(length(unique(b$mde_r_at_80pct)) == 1)

  b <- b[order(match(b$family, names(FAMILY_DISPLAY)), b$exposure,
               b$response, b$compartment), ]

  df <- data.frame(
    `Analysis`            = family_display(b$family),
    `Exposure`            = exposure_display(b$exposure),
    `Compartment`         = compartment_display(b$compartment),
    `Response`            = response_display(b$response),
    `Households`          = as.integer(b$n_households),
    `Samples`             = as.integer(b$n_samples),
    `Effect (SD)`         = round(b$estimate_SD, 3),
    `SE`                  = round(b$se, 3),
    `95% CI lower`        = round(b$ci_lo, 3),
    `95% CI upper`        = round(b$ci_hi, 3),
    `p`                   = round(b$p, 4),
    `BH-adjusted p`       = round(b$q_BH, 4),
    `Minimum detectable |r| at 80% power` = round(b$mde_r_at_80pct, 3),
    check.names = FALSE, stringsAsFactors = FALSE)

  defs <- c(
    `Analysis` = paste0("Primary tests cover the animal-exposure and antibiotic axes the study ",
      "was designed around. Exploratory tests cover water, sanitation and hygiene. Post hoc ",
      "socioeconomic tests were added after the primary results had been seen. Benjamini-Hochberg ",
      "adjustment was applied separately within each of these three analyses."),
    `Exposure` = paste0("Household exposure variable. Two names in the analysis code ",
      "misdescribe their contents and are corrected here: the cattle-and-chickens count is ",
      "the sum of cattle and chicken head only, with no goats, sheep or donkeys; and whether ",
      "cattle are kept in a pen or stable is a single binary with no day or night split."),
    `Compartment` = "Compartment the response was measured in.",
    `Response` = sh$response,
    `Households` = "Households contributing to the test.",
    `Samples` = "Samples contributing to the test.",
    `Effect (SD)` = sh$eff,
    `SE` = sh$se,
    `95% CI lower` = sh$ci_lo,
    `95% CI upper` = sh$ci_hi,
    `p` = sh$p,
    `BH-adjusted p` = paste0("The p-value after Benjamini-Hochberg correction within its own ",
      "analysis. The three were adjusted separately, so a value here is comparable within one ",
      "analysis and across the three it is not."),
    `Minimum detectable |r| at 80% power` = paste0("The smallest correlation this design ",
      "could detect at 80% power given 39 households: 0.44, or about 19% of variance. It is a ",
      "design constant and identical in every row. A null result here excludes associations ",
      "at least this large and says nothing about smaller ones."))

  tbl_spec(
    title = paste0("Associations between household and participant exposures and compartment ",
                   "gene content: 135 univariate tests."),
    df = df, definitions = defs,
    caption = c(
      paste0("Each row is one univariate model fitted to one response. Predictors were tested one at a ",
             "time."),
      paste0("Resistance-gene features are ", filt, " against the negative controls, and abundances are\n             normalized to the ", norm, "."),
      paste0("Figure S4 draws the 54 primary tests. The 45 exploratory and 36 post hoc socioeconomic ",
             "tests appear in no figure."),
      paste0("Soil physico-chemical properties are tested separately, in ",
             "exposure_associations_soil_properties.csv."),
      "Analysed cohort: 39 households, with 39 child stool, 38 cattle dung and 117 soil samples."))
}

spec_T7c <- function() {
  sh <- .t7_shared()
  d <- read_tab("T7c_soil_drivers.tsv")

  stopifnot(nrow(d) == 24,
            all(c("driver_label", "response", "n_soils", "n_households",
                  "estimate_SD", "se", "p", "q_BH", "ci_lo", "ci_hi",
                  "chemistry_source") %in% names(d)))

  filt <- unique(d$contamination_filter)
  chem <- unique(d$chemistry_source)
  stopifnot(length(filt) == 1, length(chem) == 1)

  d <- d[order(d$q_BH), ]

  df <- data.frame(
    `Soil property` = d$driver_label,
    `Response`      = response_display(d$response),
    `Soils`         = as.integer(d$n_soils),
    `Households`    = as.integer(d$n_households),
    `Effect (SD)`   = round(d$estimate_SD, 3),
    `SE`            = round(d$se, 3),
    `95% CI lower`  = round(d$ci_lo, 3),
    `95% CI upper`  = round(d$ci_hi, 3),
    `p`             = round(d$p, 5),
    `BH-adjusted p` = round(d$q_BH, 4),
    check.names = FALSE, stringsAsFactors = FALSE)

  defs <- c(
    `Soil property` = "Soil chemical property tested as a predictor.",
    `Response` = sh$response,
    `Soils` = "Soil samples contributing to the test.",
    `Households` = "Households contributing soils to the test.",
    `Effect (SD)` = sh$eff,
    `SE` = sh$se,
    `95% CI lower` = sh$ci_lo,
    `95% CI upper` = sh$ci_hi,
    `p` = sh$p,
    `BH-adjusted p` = paste0("The p-value after Benjamini-Hochberg correction across the 24 ",
      "tests in this file."))

  tbl_spec(
    title = paste0("Associations between soil physico-chemical properties and soil gene ",
                   "content: 24 univariate tests."),
    df = df, definitions = defs,
    caption = c(
      paste0("Each row is one univariate model fitted to one response. Predictors were tested one at a ",
             "time."),
      paste0("Soil chemistry values are an ", chem, ". Per-analyte reliability is reported in Table S6."),
      paste0("Resistance-gene features are ", filt, " against the negative controls."),
      paste0("Household and participant exposures are tested separately, in ",
             "exposure_associations_household.csv."),
      "Analysed cohort: 39 households and 117 soil samples."))
}

# T42c -- FEAST rarefaction-depth sensitivity

spec_T42c <- function() {
  prim <- read_tab("T42b_source_contrast_by_stratum.tsv")
  sens <- read_tab("T42c_sensitivity.tsv")

  stopifnot(all(c("stratum", "n_household", "contrast", "p") %in% names(prim)),
            all(c("analysis", "setting", "layer", "stratum", "n_household",
                  "contrast", "p", "note") %in% names(sens)))

  sens <- sens[sens$analysis == "rarefaction depth", , drop = FALSE]
  if (!nrow(sens)) stop("no rarefaction-depth rows in T42c", call. = FALSE)

  stopifnot(length(unique(sens$layer)) == 1,
            !any(grepl("taxonomy|species", sens$analysis, ignore.case = TRUE)))

  strata <- c("gut to soil", "dung to soil")
  stopifnot(all(strata %in% sens$stratum), all(strata %in% prim$stratum))

  depth <- as.integer(sub(".*?([0-9]+).*", "\\1", sens$setting))
  stopifnot(!any(is.na(depth)))

  pieces <- lapply(strata, function(st) {
    pr <- prim[prim$stratum == st, , drop = FALSE]
    se <- sens[sens$stratum == st, , drop = FALSE]
    d  <- depth[sens$stratum == st]
    o  <- order(-d)
    data.frame(
      stratum     = st,
      depth       = c(PRIMARY_COVERAGE, d[o]),
      is_primary  = c(TRUE, rep(FALSE, nrow(se))),
      n_household = c(pr$n_household, se$n_household[o]),
      contrast    = c(pr$contrast, se$contrast[o]),
      p           = c(pr$p, se$p[o]),
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, pieces)

  df <- data.frame(
    `Comparison`         = pair_display(sub(" to ", "-", d$stratum)),
    `Rarefaction depth`  = ifelse(d$is_primary,
                                  paste0(fmt_int(d$depth), " (primary)"),
                                  fmt_int(d$depth)),
    `Households`         = fmt_int(d$n_household),
    `Contrast`           = fmt_signed(d$contrast, 4),
    `p`                  = fmt_p(d$p, floor_at = 1 / (N_PERM_FEAST + 1)),
    check.names = FALSE, stringsAsFactors = FALSE)

  tbl_spec(
    title = paste0("Table S14. FEAST own-household minus decoy-household source ",
                   "contrast at three rarefaction depths."),
    df = df,
    caption = c(
      paste0("For each sink sample, own-household and decoy-household source samples were supplied to ",
             "the same FEAST run at equal sample counts. The contrast is the estimated own-household ",
             "source proportion minus the decoy-household proportion, averaged to the household before ",
             "testing."),
      paste0("Households is the number contributing to the contrast. p is a sign-flip permutation ",
             "p-value over 9,999 permutations, with a lower bound of 0.0001."),
      paste0("Rarefaction depth is a per-sample count-sum target. 848 is the primary depth, and is the ",
             "total in the shallowest sample of the ARG-variant table. Depth was varied downward only: a ",
             "sample already below the target is passed through unrarefied."),
      "Five decoy households were drawn per sink, with three rarefaction replicates per draw."))
}

# T40 -- fixed similarity statistics against the trained classifier, dung-soil

spec_T40 <- function() {
  fx <- read_tab("T40c_fixed_statistic_calibration.tsv")
  cl <- read_tab("T40d_pair_classifier.tsv")

  stopifnot(all(c("layer", "statistic", "stratum", "AUC", "null_sd", "p",
                  "n_pos", "n_pairs") %in% names(fx)),
            all(c("spec", "stratum", "layer", "AUC", "null_sd", "p",
                  "detectable_AUC_80pct", "n_pos", "n_pairs") %in% names(cl)))

  LAYER   <- "ARG variant (clinical)"
  STRATUM <- "dung to soil"

  fx <- fx[fx$layer == LAYER & fx$stratum == STRATUM, , drop = FALSE]
  cl <- cl[cl$layer == LAYER & cl$stratum == STRATUM & cl$spec == "primary", , drop = FALSE]
  if (!nrow(fx) || nrow(cl) != 1L)
    stop("expected the variant-layer dung-soil rows and exactly one primary ",
         "classifier row; got ", nrow(fx), " and ", nrow(cl), call. = FALSE)

  stopifnot(all(fx$layer == LAYER), all(cl$layer == LAYER))

  stopifnot(length(unique(c(fx$n_pairs, cl$n_pairs))) == 1L,
            length(unique(c(fx$n_pos,   cl$n_pos)))   == 1L)

  d <- data.frame(
    statistic = c(fx$statistic, "Trained classifier (household-grouped cross-validation)"),
    AUC       = c(fx$AUC, cl$AUC),
    null_sd   = c(fx$null_sd, cl$null_sd),
    p         = c(fx$p, cl$p),
    mde       = c(rep(NA_real_, nrow(fx)), cl$detectable_AUC_80pct),
    learned   = c(rep(FALSE, nrow(fx)), TRUE),
    stringsAsFactors = FALSE)
  d <- d[order(d$learned, -d$AUC), , drop = FALSE]

  pp <- character(nrow(d))
  pp[!d$learned] <- fmt_p(d$p[!d$learned], floor_at = 1 / (N_PERM_FEAST + 1))
  pp[ d$learned] <- fmt_p(d$p[ d$learned], floor_at = 1 / (N_PERM_CLASSIFY + 1))

  df <- data.frame(
    `Statistic`                = statistic_display(d$statistic),
    `AUC`                      = fmt_num(d$AUC, 3),
    `Null SD`                  = fmt_num(d$null_sd, 3),
    `p`                        = pp,
    `Detectable AUC (80% power)` = fmt_num(d$mde, 3),
    check.names = FALSE, stringsAsFactors = FALSE)

  tbl_spec(
    title = paste0("Table S15. Direct measures of pairwise similarity and a trained classifier, ",
                   "applied to the same cattle dung and soil sample pairs."),
    df = df,
    caption = c(
      paste0("All rows are computed on the same 4,446 cattle dung and soil sample pairs, 114 of which ",
             "are same-household. AUC is the probability that a same-household pair scores above a ",
             "different-household pair."),
      paste0("The upper seven rows are computed from the sample pairs directly, without a model: ",
             "co-detected-variant richness, then Jaccard and co-presence over variants detected in a ",
             "quarter of households or fewer and in a tenth or fewer, then Bray-Curtis similarity. ",
             "Co-presence counts the variants both samples carry; Jaccard is that count as a ",
             "proportion of the variants either carries. Each was tested against 995 permutations of ",
             "the household labels."),
      paste0("The trained classifier is a gradient-boosted model with household-level ",
             "cross-validation, tested against 199 permutations, each re-running the complete fitting ",
             "procedure. It was supplied the same summary quantities as features alongside the ",
             "per-variant indicators, so it has access to everything the direct measures use."),
      paste0("The joint rare-variant richness row is a confound probe: it counts rare variants ",
             "present in either sample without requiring both, so a household signal there would ",
             "indicate a shared detectability artifact."),
      "Null SD is the standard deviation of the permutation distribution for that statistic.",
      paste0("Detectable AUC at 80% power is given for the classifier alone, at the null mean plus ",
             "2.802 null standard deviations. Six of the direct measures record values above chance ",
             "that fall below that threshold.")))
}

# T1a -- cohort inventory and sequencing depth  [supplement: Table S8]

spec_T1a <- function() {
  d <- read_tab("T1a_compartment_inventory.tsv")

  stopifnot(all(c("compartment", "n_key", "n_analyzed", "n_households",
                  "depth_min", "depth_mean", "depth_median", "depth_max") %in% names(d)))
  d <- d[match(COMP_LEVELS, d$compartment), ]
  stopifnot(!any(is.na(d$compartment)))

  M <- function(x) x / 1e6

  df <- data.frame(
    `Compartment`               = compartment_display(d$compartment),
    `Samples collected`         = fmt_int(d$n_key),
    `Samples analysed`          = fmt_int(d$n_analyzed),
    `Households`                = fmt_int(d$n_households),
    `Mean depth (M reads)`      = fmt_num(M(d$depth_mean), 2),
    `Median depth (M reads)`    = fmt_num(M(d$depth_median), 2),
    `Range (M reads)`           = sprintf("%.2f–%.2f", M(d$depth_min), M(d$depth_max)),
    check.names = FALSE, stringsAsFactors = FALSE)

  tbl_spec(
    title = paste0("Table S8. Analysed cohort and post-host-filtering sequencing depth ",
                   "by compartment."),
    df = df,
    spanner = c("", "Samples", "Samples", "Samples",
                "Sequencing depth", "Sequencing depth", "Sequencing depth"),
    caption = c(
      paste0("Samples collected is the biological cohort. Samples analysed is the subset present in the ",
             "gene-catalog feature tables, which every analysis in this study uses. The two differ for ",
             "cattle dung alone: one library returned 25,052 reads and was excluded, leaving complete ",
             "sample sets for 38 of the 39 households. A further 35 soil samples, collected in five ",
             "households for a separate pilot analysis, sit outside this cohort. The analysed cohort is ",
             "194 samples."),
      paste0("Depth is total reads after adapter trimming, quality filtering and host depletion, ",
             "counting read 1 and read 2 separately. Range is the minimum and maximum single library of ",
             "that compartment."),
      paste0("Soil is sampled at three microsites per household. Child stool and cattle dung are one ",
             "sample per household.")))
}

# T7d -- per-analyte reliability of the ICRAF predictions  [supplement: Table S6]

spec_T7d <- function() {
  d <- read_tab("T7e_soil_analyte_reliability.tsv")

  stopifnot(all(c("icraf_property", "analyte", "model", "icraf_r2_val",
                  "icraf_rmse_val", "sd_this_cohort", "reliability",
                  "screen_outcome") %in% names(d)))
  stopifnot(sum(d$screen_outcome == "retained") == 8,
            sum(d$screen_outcome == "excluded") == 5)

  stopifnot(sum(d$screen_outcome == "not measured in this cohort") == 2)
  d <- d[d$screen_outcome != "not measured in this cohort", , drop = FALSE]

  fmt_reliability <- function(x)
    ifelse(is.na(x), EMDASH,
           ifelse(x <= 0, "≤ 0", sprintf("%.0f%%", 100 * x)))

  SCREEN_DISPLAY <- c(
    "retained"                    = "Retained",
    "excluded"                    = "Excluded",
    "not tested (soil texture)"   = "Not a candidate (texture)",
    "not measured in this cohort" = "Not in our data")

  df <- data.frame(
    `Analyte`                     = ifelse(is.na(d$analyte), EMDASH, d$analyte),
    `Property as reported`        = d$icraf_property,
    `Model`                       = d$model,
    `Validation R²`          = fmt_num(d$icraf_r2_val, 3),
    `Validation RMSE`             = fmt_num(d$icraf_rmse_val, 3),
    `SD in these soils`           = ifelse(is.na(d$sd_this_cohort), EMDASH,
                                           formatC(d$sd_this_cohort, format = "f",
                                                   digits = 3, big.mark = ",")),
    `Variance resolved`           = fmt_reliability(d$reliability),
    `Screen`                      = unname(SCREEN_DISPLAY[d$screen_outcome]),
    check.names = FALSE, stringsAsFactors = FALSE)

  stopifnot(!any(is.na(df$Screen)))

  tbl_spec(
    title = paste0("Table S6. Reliability of the spectroscopic soil-property predictions, ",
                   "and the analytes retained for analysis."),
    df = df,
    spanner = c("", "", "", "ICRAF global model", "ICRAF global model",
                "These soils", "These soils", ""),
    caption = c(
      paste0("Soil chemistry values in this study are mid-infrared spectroscopy predictions from ",
             "CIFOR-ICRAF's global soil spectral library. Rows are ordered by the variance each ",
             "prediction resolves in these soils."),
      paste0("Validation R2 and RMSE are ICRAF's own, measured on a 30% hold-out of ICRAF's calibration ",
             "library, and are reproduced from Table 2 of their report. Their stated criterion is R2 of ",
             "0.6 or above. These statistics describe the models and say nothing about these soils."),
      paste0("Variance resolved is 1 - (RMSE / SD)^2, computed here from ICRAF's validation RMSE and the ",
             "standard deviation of the 152 soils in this study. It appears nowhere in ICRAF's report. ",
             "The two statistics order the analytes differently: potassium is 0.698 by R2 and 95% here, ",
             "and aluminium is 0.855 by R2 and 13% here."),
      paste0("Where prediction error exceeds the spread of these soils the expression falls below zero ",
             "and the prediction resolves none of the local variation. Such values are shown as 0 or ",
             "below."),
      paste0("Eight analytes resolved 69% or more and were carried into the association analyses; five ",
             "were excluded. Clay and sand are soil texture and were never candidates. Sulfur and copper ",
             "were returned by ICRAF and are not columns in our per-sample data.")))
}

# T31a -- compositional exposure tests  [supplementary data]

spec_T31a <- function() {
  d <- read_tab("T31a_exposure_compositional_associations.tsv")

  stopifnot(nrow(d) == 180,
            all(c("compartment", "layer", "distance", "exposure", "exposure_label",
                  "family", "type", "n_samples", "n_households", "n_exposed", "R2",
                  "R2_null_expected", "R2_excess", "permutation", "p",
                  "dispersion_p", "q_BH") %in% names(d)))

  stopifnot(all(is.na(d$q_BH[d$distance == "Aitchison"])),
            !any(is.na(d$q_BH[d$distance == "Bray-Curtis"])),
            sum(d$distance == "Bray-Curtis") == 90)
  stopifnot(setequal(unique(d$layer), c("function (KO)", "resistome (ARG)")))

  filt <- unique(d$contamination_filter)
  stopifnot(length(filt) == 1)

  stopifnot(min(d$q_BH, na.rm = TRUE) >= 0.05)

  d <- d[order(match(d$family, names(FAMILY_DISPLAY)), d$distance, d$exposure,
               d$compartment, d$layer), ]

  df <- data.frame(
    `Analysis`          = family_display(d$family),
    `Exposure`          = d$exposure_label,
    `Predictor type`    = d$type,
    `Compartment`       = compartment_display(d$compartment),
    `Feature layer`     = d$layer,
    `Distance`          = d$distance,
    `Samples`           = as.integer(d$n_samples),
    `Households`        = as.integer(d$n_households),
    `Exposed`           = as.integer(d$n_exposed),
    `R2`                = round(d$R2, 4),
    `Expected R2`       = round(d$R2_null_expected, 4),
    `Excess R2`         = round(d$R2_excess, 4),
    `Permutation`       = d$permutation,
    `p`                 = round(d$p, 4),
    `Dispersion p`      = round(d$dispersion_p, 4),
    `BH-adjusted p`     = round(d$q_BH, 4),
    check.names = FALSE, stringsAsFactors = FALSE)

  defs <- c(
    `Analysis` = paste0("Primary tests cover the animal-exposure and antibiotic axes the study ",
      "was designed around. Exploratory tests cover water, sanitation and hygiene. Post hoc ",
      "socioeconomic tests were added after the primary results had been seen. Benjamini-Hochberg ",
      "adjustment was applied separately within each of these three analyses."),
    `Exposure` = "Household exposure variable, tested one at a time rather than jointly.",
    `Predictor type` = "Binary or continuous, which determines the permutation null used.",
    `Compartment` = "Compartment whose composition was tested.",
    `Feature layer` = paste0("Gene content layer the distance matrix was built on: functional ",
      "orthologs, resistance genes, or genus-level taxonomy."),
    `Distance` = paste0("Bray-Curtis is the primary distance and carries the multiple-testing ",
      "adjustment. Aitchison is Euclidean distance on the centred log-ratio and is reported ",
      "beside it as an unadjusted companion on the compositional choice."),
    `Samples` = "Samples contributing to the test.",
    `Households` = "Households contributing to the test.",
    `Exposed` = paste0("Samples in the exposed group, for binary predictors. NA for continuous ",
      "predictors, where the term does not apply."),
    `R2` = "Proportion of variation in the distance matrix attributed to the exposure.",
    `Expected R2` = paste0("The value a random grouping of the same number and size produces by ",
      "construction. An R2 must be read against it, never alone."),
    `Excess R2` = "R2 minus expected R2. This is the effect size.",
    `Permutation` = paste0("Soil has three microsites per household, so permutation is ",
      "restricted within households to preserve those triples. Child stool and cattle dung have ",
      "one sample per household and permute freely."),
    `p` = paste0("Permutation p over 999 permutations. It cannot fall below 1/(999 + 1) = ",
      "0.001; a value at that bound is the resolution limit of the test, not a tail probability."),
    `Dispersion p` = paste0("Test of whether the groups differ in spread rather than in location ",
      "(vegan::betadisper). A diagnostic, not a result: a small value means the main p should ",
      "not be read as a difference in composition alone."),
    `BH-adjusted p` = paste0("Benjamini-Hochberg adjusted p, applied within each of the three analyses named in the ",
      "first column and across the Bray-Curtis rows only. NA on every Aitchison row, which is an ",
      "unadjusted companion."))

  tbl_spec(
    title = paste0("Associations between household exposures and whole compartment composition: ",
                   "180 permutation tests."),
    df = df, definitions = defs,
    caption = c(
      paste0("Each row is one exposure permuted against one distance matrix. These are the multivariate ",
             "counterparts of the univariate tests in exposure_associations_household.csv."),
      paste0("90 tests were run under each distance: 45 on functional orthologs and 45 on resistance ",
             "genes, in each case across the three compartments. No test reaches an adjusted p below ",
             "0.05, and the smallest raw p is 0.022."),
      paste0("Resistance-gene features are ", filt, " against the negative controls; KO features are not\n             filtered."),
      "Analysed cohort: 39 households, with 39 child stool, 38 cattle dung and 117 soil samples."))
}

# T39b -- child covariate tests  [supplementary data]

spec_T39b <- function() {
  d <- read_tab("T39b_child_covariate_associations.tsv")

  stopifnot(nrow(d) == 18,
            all(c("analysis", "response", "covariate", "n_households", "estimate_SD",
                  "se", "p", "ci_lo", "ci_hi", "q_BH", "mde_r_at_80pct",
                  "dob_source") %in% names(d)))
  stopifnot(length(unique(d$mde_r_at_80pct)) == 1,
            length(unique(d$dob_source)) == 1,
            all(d$analysis == "primary"))

  filt <- unique(d$contamination_filter)
  norm <- unique(d$normalization)
  dob  <- unique(d$dob_source)
  stopifnot(length(filt) == 1, length(norm) == 1)

  d <- d[order(d$q_BH, d$p), ]

  df <- data.frame(
    `Covariate`     = covariate_display(d$covariate),
    `Response`      = response_display(d$response),
    `Households`    = as.integer(d$n_households),
    `Effect (SD)`   = round(d$estimate_SD, 3),
    `SE`            = round(d$se, 3),
    `95% CI lower`  = round(d$ci_lo, 3),
    `95% CI upper`  = round(d$ci_hi, 3),
    `p`             = round(d$p, 4),
    `BH-adjusted p` = round(d$q_BH, 4),
    `Minimum detectable |r| at 80% power` = round(d$mde_r_at_80pct, 3),
    check.names = FALSE, stringsAsFactors = FALSE)

  sh <- .t7_shared()
  defs <- c(
    `Covariate` = paste0("Child-level covariate, measured once per household. Child sex is ",
      "coded 1 = female, so the sign of its effect is relative to male children."),
    `Response` = sh$response,
    `Households` = "Households contributing to the test. All 39 contribute to every test.",
    `Effect (SD)` = sh$eff,
    `SE` = sh$se,
    `95% CI lower` = sh$ci_lo,
    `95% CI upper` = sh$ci_hi,
    `p` = sh$p,
    `BH-adjusted p` = "Benjamini-Hochberg adjusted p across the 18 tests in this file.",
    `Minimum detectable |r| at 80% power` = paste0("The smallest correlation this design could ",
      "detect at 80% power given 39 households: 0.44, or about 19% of variance. It is a design ",
      "constant and identical in every row. A null result here excludes associations at least ",
      "this large and says nothing about smaller ones."))

  tbl_spec(
    title = "Associations between child covariates and child stool gene content: 18 tests.",
    df = df, definitions = defs,
    caption = c(
      paste0("Each row is one univariate model fitted to one response. Covariates were tested one at a ",
             "time."),
      paste0("Child age is derived from the date of birth recorded at ", dob, ". A second date of birth\n             was recorded in the demographic module and the two disagree in some households; neither\n             can be shown correct from the data. The screening date is used throughout, and the\n             alternative is run as a full sensitivity analysis."),
      paste0("Three further covariates were constructed and failed the feasibility screen on sample size ",
             "or group balance: current breastfeeding, daytime location and length of residence. They ",
             "are reported in the screen."),
      paste0("Resistance-gene features are ", filt, " against the negative controls, and abundances are\n             normalized to the ", norm, "."),
      "Responses are measured in child stool only. All 39 households contribute."))
}

# T45 -- dietary tests, both families  [supplementary data]

spec_T45 <- function() {
  a <- read_tab("T45a_diet_recall_child.tsv")
  b <- read_tab("T45b_diet_ffq_household.tsv")

  stopifnot(nrow(a) == 12, nrow(b) == 8,
            identical(sort(names(a)), sort(names(b))))
  d <- rbind(a, b[, names(a)])
  stopifnot(all(c("arm", "response", "variable", "n", "beta", "ci_lo", "ci_hi",
                  "R2", "R2_floor", "R2_excess", "dispersion_p", "p", "q_BH",
                  "family") %in% names(d)),
            setequal(unique(d$family), c("A", "B")))

  sc <- d$arm == "scalar"
  stopifnot(all(!is.na(d$beta[sc])), all(is.na(d$R2[sc])),
            all(is.na(d$beta[!sc])), all(!is.na(d$R2[!sc])))


  stopifnot(min(d$q_BH, na.rm = TRUE) >= 0.05)

  d <- d[order(d$family, d$arm, d$q_BH, d$p), ]

  df <- data.frame(
    `Dietary analysis` = diet_family_display(d$family),
    `Measure`         = diet_display(d$variable),
    `Analysis`        = diet_arm_display(d$arm),
    `Response`        = response_display(d$response),
    `Households`      = as.integer(d$n),
    `Effect (SD)`     = round(d$beta, 3),
    `95% CI lower`    = round(d$ci_lo, 3),
    `95% CI upper`    = round(d$ci_hi, 3),
    `R2`              = round(d$R2, 4),
    `R2 floor`        = round(d$R2_floor, 4),
    `Excess R2`       = round(d$R2_excess, 4),
    `Dispersion p`    = round(d$dispersion_p, 4),
    `p`               = round(d$p, 4),
    `BH-adjusted p`   = round(d$q_BH, 4),
    check.names = FALSE, stringsAsFactors = FALSE)

  sh <- .t7_shared()
  defs <- c(
    `Dietary analysis` = paste0("Analysis A covers what the child ate in the 24 hours before the ",
      "visit. Analysis B covers what the household reports eating habitually. Benjamini-Hochberg ",
      "adjustment was applied separately within each, so an adjusted p is comparable within one ",
      "analysis and across the two it is not."),
    `Measure` = paste0("Dietary variable tested. The two milk-product items cover cheese, ",
      "yogurt, kefir, milk and other milk products together, and neither identifies where the ",
      "food came from."),
    `Analysis` = paste0("A scalar analysis regresses one summary of the gene content on the ",
      "dietary measure. A gut composition analysis permutes the measure against the whole ",
      "Bray-Curtis distance matrix. The two carry different columns."),
    `Response` = sh$response,
    `Households` = "Households contributing to the test.",
    `Effect (SD)` = paste0(sh$eff, " NA on gut composition rows, which have no single ",
      "coefficient."),
    `95% CI lower` = paste0(sh$ci_lo, " NA on gut composition rows."),
    `95% CI upper` = paste0(sh$ci_hi, " NA on gut composition rows."),
    `R2` = paste0("Proportion of variation in the distance matrix attributed to the dietary ",
      "measure. NA on scalar rows."),
    `R2 floor` = paste0("The R2 a single predictor of this kind reaches by construction at this ",
      "sample size. An R2 must be read against it, never alone. NA on scalar rows."),
    `Excess R2` = "R2 minus the floor. This is the effect size. NA on scalar rows.",
    `Dispersion p` = paste0("Test of whether the groups differ in spread rather than in location. ",
      "A diagnostic, not a result. NA on scalar rows."),
    `p` = sh$p,
    `BH-adjusted p` = paste0("Benjamini-Hochberg adjusted p within the dietary analysis named in the first column."))

  tbl_spec(
    title = "Associations between recent and habitual dietary intake and child stool gene content: 20 tests.",
    df = df, definitions = defs,
    caption = c(
      paste0("Twelve tests cover what the child ate in the 24 hours before the visit and eight cover the ",
             "household's habitual intake. Each measure was run against three scalar responses and ",
             "against whole gut composition. No test reaches an adjusted p below 0.05."),
      paste0("The milk-product items ask about cheese, yogurt, kefir, milk or other milk products as one ",
             "category. Food sourcing is recorded nowhere in this study, so neither item distinguishes ",
             "the household's own cattle from market or neighbour supply."),
      paste0("NA marks a column that does not apply to that row: scalar analyses have no R2 and ",
             "composition analyses have no single coefficient."),
      "Responses are measured in child stool only. All 39 households contribute."))
}

# T20 -- clinical subset definition and sensitivity  [supplement: Table S1]

spec_T20 <- function() {
  d <- read_tab("T20_clinical_subset_sensitivity.tsv")

  stopifnot(nrow(d) == 4,
            all(c("def", "definition", "n_symbols", "n_symbols_detected", "n_alleles",
                  "rich_gut", "rich_dung", "rich_soil", "cpg_gut", "cpg_dung",
                  "cpg_soil", "gut_over_soil") %in% names(d)))


  stopifnot(d$def[1] == "A", d$n_symbols_detected[1] == 170, d$n_alleles[1] == 1681)

  rr <- range(d$gut_over_soil)
  stopifnot(abs(rr[1] - 5.76) < 0.01, abs(rr[2] - 9.39) < 0.01)

  df <- data.frame(
    `Definition`                = sprintf("%s. %s", d$def,
                                          sub(" \\(current\\)$", "", d$definition)),
    `ARG families`              = fmt_int(d$n_symbols_detected),
    `ARG variants`              = fmt_int(d$n_alleles),
    `Child stool`               = fmt_num(d$cpg_gut, 3),
    `Cattle dung`               = fmt_num(d$cpg_dung, 3),
    `Soil`                      = fmt_num(d$cpg_soil, 3),
    `Child stool / soil`        = sprintf("%.2f×", d$gut_over_soil),
    `Child stool richness`      = fmt_int(d$rich_gut),
    `Cattle dung richness`      = fmt_int(d$rich_dung),
    `Soil richness`             = fmt_int(d$rich_soil),
    check.names = FALSE, stringsAsFactors = FALSE)

  tbl_spec(
    title = paste0("Table S1. Clinically relevant resistance-gene abundance and richness under ",
                   "four definitions of the clinical subset."),
    df = df,
    spanner = c("", "Feature set", "Feature set",
                "Median copies per average genome", "Median copies per average genome",
                "Median copies per average genome", "Median copies per average genome",
                "Median richness", "Median richness", "Median richness"),
    caption = c(
      paste0("Definition A is reported throughout the manuscript: AMRFinderPlus elements of type AMR at ",
             "core scope. B adds nine chromosomal efflux and beta-lactamase determinants called at ",
             "--plus scope. C adds the three core quaternary-ammonium biocide determinants. D adds both. ",
             "Virulence determinants are excluded under every definition."),
      paste0("The child stool to soil ratio spans ", sprintf("%.2f", rr[1]), " to ", sprintf("%.2f", rr[2]),
             " across the four definitions and should be read as a range. The ordering of child stool\n             above both environmental compartments, and the complete separation of child stool from\n             them, hold under all four."),
      paste0("ARG families are AMRFinderPlus element symbols detected in the 194 analysed samples after ",
             "contamination filtering. ARG variants are the catalog gene representatives, clustered at ",
             "95% identity, belonging to those families."),
      paste0("One construction is held fixed here, so the four rows are comparable to each other. ",
             "Equal-effort subsampling is not applied, and these values are not comparable to the levels ",
             "reported in the main text.")))
}

# T44a / T44b -- what contamination filtering removed  [supplement: Table S2]

spec_T44 <- function() {
  a <- read_tab("T44a_contamination_filter_features.tsv")
  b <- read_tab("T44b_contamination_filter_contrast.tsv")
  f <- read_tab("T_arg_abundance_filter_sensitivity.tsv")

  stopifnot(nrow(a) == 4, nrow(b) == 2, nrow(f) == 3,
            all(c("resolution", "feature_set", "exclusion_list", "n_before",
                  "n_excluded", "n_after", "pct_removed") %in% names(a)),
            all(c("contrast", "ratio_filtered", "ratio_unfiltered") %in% names(b)))

  cap1 <- function(x) sub("^(.)", "\\U\\1", x, perl = TRUE)
  df_a <- data.frame(
    `Resolution`      = a$resolution,
    `Feature set`     = cap1(a$feature_set),
    `Exclusion list`  = cap1(a$exclusion_list),
    `Detected`        = fmt_int(a$n_before),
    `Excluded`        = fmt_int(a$n_excluded),
    `Retained`        = fmt_int(a$n_after),
    `Removed (%)`     = fmt_num(a$pct_removed, 1),
    check.names = FALSE, stringsAsFactors = FALSE)

  f <- f[match(COMP_LEVELS, f$compartment), ]
  df_b <- data.frame(
    `Comparison`  = c(compartment_display(f$compartment), b$contrast),
    `Filtered`    = c(fmt_num(f$clinical_filtered, 3),
                      sprintf("%.2f×", b$ratio_filtered)),
    `Not filtered` = c(fmt_num(f$clinical_unfiltered, 3),
                       sprintf("%.2f×", b$ratio_unfiltered)),
    check.names = FALSE, stringsAsFactors = FALSE)

  part_a <- tbl_spec(
    title = "(A) Features excluded, at each resolution.",
    df = df_a,
    caption = c(
      paste0("A feature was excluded where it reached 5 or more reads in any negative extraction or PCR ",
             "control library. The threshold was set on the ARG detection curve, which collapses to zero ",
             "by 50 reads."),
      paste0("Two exclusion lists were used and were never crossed. The family-level list keys on the ",
             "element symbol and the gene-level list keys on the catalog gene identifier. The 15 ",
             "excluded families span 1,173 catalog genes. Both lists hold 15 entries and the entries are ",
             "different objects: 15 element symbols in one, 15 catalog genes in the other."),
      paste0("KO features are not filtered. A 1-read rule would strip about 37% of them, and they are ",
             "largely universal housekeeping functions present in any blank carrying reads.")))

  part_b <- tbl_spec(
    title = "(B) Clinical resistance-gene abundance per average genome, with and without filtering.",
    df = df_b,
    caption = c(
      paste0("The upper three rows are median copies per average genome. The lower two are the ratios ",
             "between them."),
      paste0("Filtering moves the contrast between compartments as well as each level, and it moves the ",
             "two contrasts in opposite directions: child stool to soil widens from 6.54 to 9.39, and ",
             "child stool to cattle dung narrows from 10.16 to 8.18. Proportionally more is removed from ",
             "soil than from child stool.")))

  tbl_parts(
    title = paste0("Table S2. Effect of contamination filtering on the resistance-gene feature ",
                   "space and on clinical resistance-gene abundance."),
    parts = list(part_a, part_b))
}

# T44c -- reference-based taxonomic classification  [supplement: Table S3]

spec_T44c <- function() {
  d <- read_tab("T44c_taxonomic_classification.tsv")

  stopifnot(nrow(d) == 4,
            all(c("sample_type", "compartment", "n", "pct_classified_median",
                  "pct_classified_min", "pct_classified_max",
                  "n_species_median") %in% names(d)))
  stopifnot(sum(d$n[1:3]) == 193)

  df <- data.frame(
    `Sample type`          = d$sample_type,
    `Compartment`          = c(compartment_display(d$compartment[1:3]), d$compartment[4]),
    `n`                    = fmt_int(d$n),
    `Median reads classified (%)` = fmt_num(d$pct_classified_median, 1),
    `Range (%)`            = sprintf("%.1f–%.1f", d$pct_classified_min, d$pct_classified_max),
    `Median species detected` = fmt_int(d$n_species_median),
    check.names = FALSE, stringsAsFactors = FALSE)

  tbl_spec(
    title = paste0("Table S3. Reference-based taxonomic classification rates by compartment, ",
                   "against the mock communities."),
    df = df,
    caption = c(
      paste0("MetaPhlAn 4.1.2 against the mpa_vJun23_CHOCOPhlAnSGB_202307 database. Values are the ",
             "percentage of reads in a library assigned to a species-level genome bin, summarised across ",
             "the libraries of one compartment."),
      paste0("The ZymoBIOMICS mock libraries are a defined eight-bacterium community sequenced alongside ",
             "the study samples. They establish that classification is near-complete where the community ",
             "is represented in the reference."),
      paste0("All cross-compartment analyses in this study are conducted in functional and ",
             "resistance-gene feature space, which is assembled from these libraries."),
      paste0("n is 193 here and 194 elsewhere: one soil library is present in the gene catalog and ",
             "absent from the taxonomic summary.")))
}

# markers -- 57 ribosomal-protein marker KOs  [supplement: Table S4]

spec_markers <- function() {
  d <- read.delim(file.path(FT_DIR, "single_copy_marker_kos.tsv"),
                  stringsAsFactors = FALSE)

  stopifnot(nrow(d) == 57,
            all(c("ko", "subunit", "gene_symbols", "protein", "source", "accessed") %in% names(d)))
  src <- unique(d$source); acc <- unique(d$accessed)
  stopifnot(length(src) == 1, length(acc) == 1)

  d <- d[order(d$subunit, d$ko), ]

  df <- data.frame(
    `KO identifier` = d$ko,
    `Subunit`       = d$subunit,
    `Gene symbols`  = d$gene_symbols,
    `Protein`       = d$protein,
    check.names = FALSE, stringsAsFactors = FALSE)

  n_sub <- table(d$subunit)

  tbl_spec(
    title = paste0("Table S4. The 57 bacterial ribosomal-protein orthologs used as the ",
                   "per-genome normalization panel."),
    df = df,
    caption = c(
      paste0("Retrieved from ", src, ", accessed ", acc, "."),
      paste0("Bacterial ribosomal proteins are near-universal and present at one copy per genome. The ",
             "median abundance of this panel within a sample estimates that sample's average genome copy ",
             "number, and dividing a feature's abundance by it gives copies per average genome."),
      paste0("The panel is quantified from the same gene catalog as the resistance genes it normalizes, ",
             "so the ratio cancels the catalog mapping rate, which differs by compartment. A read-based ",
             "genome-equivalent estimate does not cancel it."),
      paste0("The panel spans both ribosomal subunits: ",
             paste(sprintf("%d %s", as.integer(n_sub), tolower(names(n_sub))), collapse = " and "),
             " proteins.")))
}

# T11c -- the two per-genome denominators compared  [supplement: Table S5]

spec_T11c <- function() {
  d <- read_tab("T11c_normalization_consequence.tsv")

  stopifnot(nrow(d) == 3,
            all(c("compartment", "clinical_per_genome_singlecopy",
                  "clinical_per_genome_GE") %in% names(d)))
  d <- d[match(COMP_LEVELS, d$compartment), ]
  stopifnot(!any(is.na(d$compartment)))

  stopifnot(identical(order(d$clinical_per_genome_singlecopy, decreasing = TRUE),
                      order(d$clinical_per_genome_GE, decreasing = TRUE)))

  df <- data.frame(
    `Compartment`                        = compartment_display(d$compartment),
    `Ribosomal-protein panel`            = fmt_num(d$clinical_per_genome_singlecopy, 3),
    `MicrobeCensus genome equivalents`   = fmt_num(d$clinical_per_genome_GE, 3),
    check.names = FALSE, stringsAsFactors = FALSE)

  tbl_spec(
    title = paste0("Table S5. Clinical resistance-gene abundance per average genome under two ",
                   "normalization denominators."),
    df = df,
    spanner = c("", "Median copies per average genome", "Median copies per average genome"),
    caption = c(
      paste0("The ribosomal-protein panel is the denominator used throughout the manuscript: the median ",
             "count of 57 bacterial ribosomal-protein orthologs in the same sample, quantified from the ",
             "same gene catalog. MicrobeCensus estimates genome equivalents from the reads directly, ",
             "without reference to the catalog."),
      paste0("The two denominators give different absolute levels and the same compartment ordering. ",
             "Both are given so these values can be compared against studies using either convention."),
      paste0("The catalog-based denominator is used because the resistance-gene counts come from the ",
             "same mapping, so the ratio cancels the catalog mapping rate, which differs by compartment.")))
}

# T3b / T37b -- multivariate dispersion  [supplement: Table S10]

spec_T3b <- function() {
  d <- read_tab("T3b_dispersion.tsv")
  u <- read_tab("T37b_unweighted_dispersion.tsv")

  stopifnot(all(c("layer", "distance", "term", "F", "p", "gut", "dung", "soil") %in% names(d)),
            all(c("layer", "presence_rule", "distance", "F", "p",
                  "gut", "dung", "soil") %in% names(u)),
            nrow(d) == 4, nrow(u) == 12)

  dup <- u$presence_rule == "count >= 1" & u$distance == "Bray-Curtis (abundance)"
  ref <- d[d$distance == "Bray-Curtis", ]
  chk <- merge(u[dup, c("layer", "F")], ref[, c("layer", "F")], by = "layer",
               suffixes = c("_unweighted", "_main"))
  stopifnot(sum(dup) == 2, nrow(chk) == 2,
            all(abs(chk$F_unweighted - chk$F_main) < 1e-2))
  u <- u[!dup, ]

  bc <- d[d$distance == "Bray-Curtis", ]
  stopifnot(all(bc$gut > bc$soil), all(bc$soil > bc$dung))

  mk <- function(x, rule, dist) data.frame(
    `Feature layer` = layer_display(x$layer, wrap = FALSE),
    `Presence rule` = rule,
    `Distance`      = dist,
    `Child stool`   = fmt_num(x$gut, 3),
    `Cattle dung`   = fmt_num(x$dung, 3),
    `Soil`          = fmt_num(x$soil, 3),
    `F`             = fmt_num(x$F, 2),
    `p`             = fmt_p(x$p, floor_at = P_FLOOR),
    check.names = FALSE, stringsAsFactors = FALSE)

  df <- rbind(mk(d, rule_display("count >= 1"), d$distance),
              mk(u, rule_display(u$presence_rule),
                 sub(" \\(abundance\\)$", "", u$distance)))
  stopifnot(nrow(df) == 14)

  filt <- unique(d$contamination_filter[grepl("^ARG", d$layer)])

  tbl_spec(
    title = paste0("Table S10. Multivariate dispersion of the three compartments, ",
                   "measured as distance to the group centroid."),
    df = df,
    spanner = c("", "", "", "Mean distance to group centroid",
                "Mean distance to group centroid", "Mean distance to group centroid", "", ""),
    caption = c(
      paste0("Dispersion is the mean distance from each sample of a compartment to that compartment's ",
             "centroid, in the same distance space the PERMANOVA used. F and p test whether the three ",
             "dispersions differ, over 999 permutations."),
      paste0("PERMANOVA responds to a difference in group spread as well as group location. The ",
             "compartments differ in dispersion on almost every distance tested, so the compartment R2 ",
             "values reported alongside describe location and spread together."),
      paste0("Under Bray-Curtis, child stool is the most dispersed compartment and cattle dung the ",
             "least, in both feature layers."),
      paste0("Dispersion is assessed for the compartment term only. No household-level dispersion test ",
             "is reported."),
      "p = 0.001 is the lower bound attainable with 999 permutations.",
      paste0("Resistance-gene features are ", filt, " against the negative controls; KO features are not\n             filtered.")))
}

# versions -- software and databases  [supplement: Table S7]

spec_versions <- function() {
  d <- read.delim(file.path(DATA_DIR, "metadata", "software_versions.tsv"),
                  stringsAsFactors = FALSE, na.strings = "NA")

  stopifnot(nrow(d) >= 20,
            all(c("stage", "tool", "version", "database", "role") %in% names(d)))

  df <- data.frame(
    `Stage`    = d$stage,
    `Software` = d$tool,
    `Version`  = d$version,
    `Reference database` = ifelse(is.na(d$database), EMDASH, d$database),
    `Role`     = d$role,
    check.names = FALSE, stringsAsFactors = FALSE)

  tbl_spec(
    title = "Table S7. Software and reference databases.",
    df = df,
    caption = c(
      paste0("Versions of the external tools are those that produced the results reported here, taken ",
             "from the job scripts that ran them. R package versions were queried from the analysis ",
             "library on 2026-09-10."),
      paste0("An em dash in the reference-database column means the tool carries no external reference ",
             "database."),
      "Packages used only to render tables and figures are not listed."))
}

# assemble

SPECS <- list(T20     = spec_T20,      # Table S1
              T44     = spec_T44,      # Table S2
              T44c    = spec_T44c,     # Table S3
              markers = spec_markers,  # Table S4
              T11c    = spec_T11c,     # Table S5
              T7d     = spec_T7d,      # Table S6
              versions = spec_versions, # Table S7
              T1a     = spec_T1a,      # Table S8
              T3      = spec_T3,       # Table S9
              T3b     = spec_T3b,      # Table S10
              Talpha  = spec_Talpha,   # Table S11
              Tclin   = spec_Tclin,    # Table S12
              T8g     = spec_T8g,      # Table S13
              T42c    = spec_T42c,     # Table S14
              T40     = spec_T40,      # Table S15
              T13     = spec_T13,      # Data S1
              T7b     = spec_T7b,      # Data S2
              T39b    = spec_T39b,     # Data S3
              T45     = spec_T45,      # Data S4
              T31a    = spec_T31a,     # Data S5
              T7c     = spec_T7c)      # Data S6

FILES <- c(T20     = "TS1_clinical_subset.docx",
           T44     = "TS2_contamination_filter.docx",
           T44c    = "TS3_taxonomic_classification.docx",
           markers = "TS4_marker_kos.docx",
           T11c    = "TS5_normalization.docx",
           T7d     = "TS6_soil_analyte_reliability.docx",
           versions = "TS7_software_versions.docx",
           T1a     = "TS8_cohort_depth.docx",
           T3      = "TS9_permanova.docx",
           T3b     = "TS10_dispersion.docx",
           Talpha  = "TS11_alpha_dispersion.docx",
           Tclin   = "TS12_clinical_arg_abundance.docx",
           T8g     = "TS13_shared_and_specific.docx",
           T42c    = "TS14_feast_rarefaction_sensitivity.docx",
           T40     = "TS15_fixed_statistics_vs_classifier.docx",
           T13     = "DataS1_household_sharing_tests",
           T7b     = "DataS2_exposure_associations_household",
           T39b    = "DataS3_covariate_associations_child",
           T45     = "DataS4_dietary_associations",
           T31a    = "DataS5_exposure_associations_compositional",
           T7c     = "DataS6_exposure_associations_soil_properties")

args <- commandArgs(trailingOnly = TRUE)
args <- args[!grepl("^[0-9]+$", args)]       # run_all.R passes script numbers
ids  <- if (length(args)) args else names(SPECS)
bad  <- setdiff(ids, names(SPECS))
if (length(bad)) stop("unknown table id(s): ", paste(bad, collapse = ", "),
                      "\nknown: ", paste(names(SPECS), collapse = ", "), call. = FALSE)

tables <- lapply(ids, function(id) SPECS[[id]]())
names(tables) <- ids

ext     <- tools::file_ext(FILES[ids])
is_docx <- ext == "docx"
is_xlsx <- ext == "xlsx"
is_csv  <- ext == ""                       # no extension: a CSV pair, see FILES

stopifnot(all(is_docx | is_xlsx | is_csv))

for (i in seq_along(ids)) {
  if (is_csv[i]) {
    write_table_csv(tables[[ids[i]]], FILES[[ids[i]]])
  } else if (is_xlsx[i]) {
    write_table_xlsx(tables[[ids[i]]], FILES[[ids[i]]])
  } else {
    write_table_docx(tables[[ids[i]]], FILES[[ids[i]]])
  }
}

doc_ids <- ids[is_docx]
if (length(doc_ids) > 1) write_tables_docx(tables[doc_ids])

message("\nmanuscript tables -> ", DOCX_DIR)
if (any(is_csv)) message("supplementary data  -> ", CSV_DIR)
