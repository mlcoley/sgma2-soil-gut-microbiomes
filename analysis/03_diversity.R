#!/usr/bin/env Rscript
# Beta-diversity structure: PERMANOVA of compartment and of household within soil,
# multivariate dispersion, and the within-soil household distance figure.
# Usage: Rscript analysis/03_diversity.R   (from the repository root)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})
source("analysis/00_load_features.R")

set.seed(1)
N_PERM <- 999

log_lines <- c(sprintf("Diversity structure - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines) + 1]] <<- m }

theme_sgma <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
key <- load_key()

dist_bray <- function(M) vegdist(tss(M), method = "bray")
dist_ait  <- function(M) {
  L <- log(M + 1)
  dist(L - rowMeans(L))
}
DISTS <- list("Bray-Curtis" = dist_bray, "Aitchison (CLR)" = dist_ait)

# KO and ARG span all analyzed samples
ko  <- get_layer("ko_all", key)
GATING <- load_gating()
arg <- gate_layer(get_layer("arg", key), GATING)
logmsg("ARG filter: %s%s", GATING$note,
       if (GATING$applied) sprintf(" -> %d ARG features retained", ncol(arg$M)) else "")

# Tables with a `layer` column mix KO and ARG rows, so the filter status is per row.
stamp_gating <- function(df) {
  df$contamination_filter <- if ("layer" %in% names(df)) {
    ifelse(grepl("ARG", df$layer), gating_col(GATING), "n/a (KO not filtered)")
  } else gating_col(GATING)
  df
}
LAY <- list("KO (function)" = ko, "ARG (resistome)" = arg)
logmsg("layers: %s", paste(sprintf("%s n=%d x %d features", names(LAY),
       vapply(LAY, function(l) nrow(l$M), 1L), vapply(LAY, function(l) ncol(l$M), 1L)),
       collapse = " | "))

D <- lapply(LAY, function(l) lapply(DISTS, function(f) f(l$M)))

# ===========================================================================
# PERMANOVA
# ===========================================================================
logmsg("== PERMANOVA ==")

# One row per (layer, distance, term). R2_null_expected = (g - 1)/(n - 1) is the
# R2 a random grouping into g levels would reach; R2_excess is the observed R2
# minus that.
permanova_row <- function(d, grp, layer, distance, term, note = "") {
  grp <- droplevels(factor(grp))
  if (nlevels(grp) < 2) return(NULL)
  a  <- vegan::adonis2(d ~ grp, permutations = N_PERM, by = "terms")
  r2_null <- (nlevels(grp) - 1) / (length(grp) - 1)
  data.frame(layer = layer, distance = distance, term = term,
             n = length(grp), groups = nlevels(grp),
             R2 = round(a$R2[1], 4), R2_null_expected = round(r2_null, 4),
             R2_excess = round(a$R2[1] - r2_null, 4), F = round(a$F[1], 2),
             p = a$`Pr(>F)`[1], note = note)
}

t3 <- list()
for (nm in names(LAY)) for (dn in names(DISTS)) {
  l <- LAY[[nm]]; d <- D[[nm]][[dn]]
  t3[[length(t3) + 1]] <- permanova_row(d, l$key$compartment, nm, dn, "compartment",
                                        "soil contributes 3 samples/household")
  ks <- l$key$compartment == "soil"
  if (sum(ks) >= 10)
    t3[[length(t3) + 1]] <- permanova_row(as.dist(as.matrix(d)[ks, ks]),
                                          l$key$household[ks], nm, dn,
                                          "household (within soil)", "3 soils/household")
}
t3 <- bind_rows(t3)

# Soil-collapse sensitivity: soil replicates averaged to one household-mean
# relative-abundance profile.
collapse_soil <- function(l) {
  R <- tss(l$M)
  ks <- l$key$compartment == "soil"
  if (!any(ks)) return(NULL)
  sm <- rowsum(R[ks, , drop = FALSE], l$key$household[ks]) /
        as.vector(table(l$key$household[ks])[unique(sort(l$key$household[ks]))])
  M2 <- rbind(R[!ks, , drop = FALSE], sm)
  cmp <- c(as.character(l$key$compartment[!ks]), rep("soil", nrow(sm)))
  list(M = M2, compartment = factor(cmp, levels = COMP_LEVELS))
}
for (nm in names(LAY)) {
  cs <- collapse_soil(LAY[[nm]])
  if (is.null(cs) || nlevels(droplevels(cs$compartment)) < 2) next
  t3 <- bind_rows(t3, permanova_row(vegdist(cs$M, "bray"), cs$compartment, nm,
                                    "Bray-Curtis", "compartment (soil collapsed)",
                                    "soil averaged to 1 profile/household"))
}
save_tab(stamp_gating(t3), "T3_permanova.tsv")
for (i in seq_len(nrow(t3)))
  logmsg("PERMANOVA [%s | %s] %s: R2=%.3f, F=%.1f, p=%.3f (n=%d, %d groups)",
         t3$layer[i], t3$distance[i], t3$term[i], t3$R2[i], t3$F[i], t3$p[i],
         t3$n[i], t3$groups[i])

# ===========================================================================
# Household effect within soil: same-household vs different-household distances
# ===========================================================================
soil_pairs <- function(nm) {
  l <- LAY[[nm]]; keep <- l$key$compartment == "soil"
  if (sum(keep) < 10) return(NULL)
  m  <- as.matrix(D[[nm]][["Bray-Curtis"]])[keep, keep]
  hh <- as.character(l$key$household[keep])
  iu <- which(upper.tri(m), arr.ind = TRUE)
  data.frame(layer = nm, dist = m[iu],
             pair = ifelse(hh[iu[, 1]] == hh[iu[, 2]], "same household", "different households"))
}
sp <- bind_rows(lapply(names(LAY), soil_pairs))
sp$pair <- factor(sp$pair, levels = c("same household", "different households"))
save_tab(sp %>% group_by(layer, pair) %>%
           summarise(n_pairs = n(), median_dist = round(median(dist), 4),
                     mean_dist = round(mean(dist), 4), .groups = "drop"),
         "T3f_soil_household_distances.tsv")
sp$layer <- factor(sp$layer, levels = names(LAY))

# Half-violins scaled to the same
# maximum width, so the 117 same-household and 6,669 different-household pairs are
# comparable in shape but not in number.
HALF_W_7B <- 0.42
f7dens <- do.call(rbind, lapply(split(sp, list(sp$layer, sp$pair), drop = TRUE), function(d) {
  if (nrow(d) < 3 || length(unique(d$dist)) < 2) return(NULL)
  dn <- density(d$dist, n = 256)
  data.frame(layer = d$layer[1], pair = d$pair[1],
             dist = c(dn$x, rev(dn$x)),
             off  = c(dn$y / max(dn$y) * HALF_W_7B, rep(0, length(dn$x))))
}))
f7dens$ypos <- as.integer(f7dens$pair) + f7dens$off
f7div <- data.frame(layer = factor("ARG (resistome)", levels = names(LAY)), xint = -Inf)
f7b <- ggplot(sp, aes(dist, as.integer(pair))) +
  geom_vline(data = f7div, aes(xintercept = xint), inherit.aes = FALSE,
             colour = "grey70", linewidth = .5) +
  geom_polygon(data = f7dens, aes(x = dist, y = ypos, group = interaction(layer, pair),
                                  fill = pair), colour = NA, alpha = .55) +
  geom_boxplot(aes(group = pair), width = .18, outlier.shape = NA, fill = "white",
               colour = "grey25", linewidth = .35, show.legend = FALSE) +
  facet_wrap(~layer, scales = "free_x", labeller = layer_labeller()) +
  scale_y_continuous(breaks = seq_along(levels(sp$pair)), labels = levels(sp$pair),
                     limits = c(0.6, length(levels(sp$pair)) + HALF_W_7B + 0.25)) +
  scale_fill_manual(values = c("same household" = COMP_COLS[["soil"]],
                               "different households" = "grey70"), guide = "none") +
  labs(x = "Bray-Curtis distance", y = NULL) +
  theme_sgma +
  theme(panel.spacing.x = unit(1.4, "lines"))
save_fig(f7b, "F7b_soil_household_effect.png", width = 9.5, height = 3.4)

# ===========================================================================
# Multivariate dispersion (betadisper) for compartment PERMANOVA
# ===========================================================================
t3b <- list()
for (nm in names(LAY)) for (dn in names(DISTS)) {
  l <- LAY[[nm]]; grp <- droplevels(l$key$compartment)
  if (nlevels(grp) < 2) next
  bd <- vegan::betadisper(D[[nm]][[dn]], grp)
  pt <- vegan::permutest(bd, permutations = N_PERM)
  t3b[[length(t3b) + 1]] <- data.frame(
    layer = nm, distance = dn, term = "compartment",
    F = round(pt$tab$F[1], 2), p = pt$tab$`Pr(>F)`[1],
    t(round(tapply(bd$distances, grp, mean), 3)))
}
t3b <- bind_rows(t3b)
save_tab(stamp_gating(t3b), "T3b_dispersion.tsv")
logmsg("dispersion (betadisper) significant in %d of %d compartment tests",
       sum(t3b$p < 0.05), nrow(t3b))

writeLines(log_lines, file.path(LOG_DIR, "03_diversity.log"))
cat("\nDiversity structure complete.\n")
