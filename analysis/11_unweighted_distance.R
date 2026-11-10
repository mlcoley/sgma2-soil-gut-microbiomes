#!/usr/bin/env Rscript
# Compartment separation on presence-only versus abundance-weighted distances.
# Bray-Curtis (on relative abundance) is compared with Jaccard and Sorensen, under both presence def
# for the KO and ARG-family layers.
# Outputs: T37a (PERMANOVA), T37b (dispersion), T37c (Procrustes agreement between the
# Bray-Curtis and presence-only ordinations), T37d (depth and richness), T37e (depth by compartment, Kruskal-Wallis), and Figure 1.
# Inputs: feature_tables/{ko,arg}_abundance_wide.tsv.gz and the *_cf50 versions.
# Usage: Rscript analysis/11_unweighted_distance.R   (from the repository root)

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })
source("analysis/00_load_features.R")

set.seed(1)
N_PERM <- 999

log_lines <- c(sprintf("Unweighted distance - %s", Sys.time()), "")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); log_lines[[length(log_lines) + 1]] <<- m }

theme_sgma <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

key    <- load_key()
GATING <- load_gating()

## ---- distances -----------------------------------------------------------
# Jaccard and Sorensen are monotone transforms of each other (J = S / (2 - S)); both
# are computed because PERMANOVA partitions the distance values and returns slightly
# different R2. Jaccard is in the figure.
dist_bray  <- function(M) vegdist(tss(M), method = "bray")
dist_jacc  <- function(M) vegdist(M, method = "jaccard", binary = TRUE)
dist_soren <- function(M) vegdist(M, method = "bray",    binary = TRUE)
DISTS <- list("Bray-Curtis (abundance)" = dist_bray,
              "Jaccard (presence)"      = dist_jacc,
              "Sorensen (presence)"     = dist_soren)

pcoa <- function(d, k = 2) {
  cm  <- cmdscale(d, k = k, eig = TRUE)
  pos <- cm$eig[cm$eig > 0]
  list(pts = cm$points, var = 100 * cm$eig[seq_len(k)] / sum(pos))
}

## ---- matrices: 2 layers x 2 presence defs -------------------------------
# The cf50 tables are not registered in LAYERS
align_to_key <- function(M) {
  common <- intersect(rownames(M), key$sample_id)
  M <- M[common, , drop = FALSE]
  M[, colSums(M) > 0, drop = FALSE]
}

ko_std  <- get_layer("ko_all", key)$M
arg_std <- gate_layer(get_layer("arg", key), GATING)$M

cf_files <- c(ko  = file.path(FT_DIR, "ko_abundance_wide_cf50.tsv.gz"),
              arg = file.path(FT_DIR, "arg_abundance_wide_cf50.tsv.gz"))
if (!all(file.exists(cf_files)))
  stop("breadth (cf50) tables missing: ", paste(cf_files[!file.exists(cf_files)], collapse = ", "))
ko_cf  <- align_to_key(load_matrix(cf_files[["ko"]],  "woltka"))
arg_cf <- align_to_key(drop_gated(load_matrix(cf_files[["arg"]], "woltka"), GATING))

MATS <- list(
  list(layer = "KO (function)",   rule = "count >= 1",              M = ko_std),
  list(layer = "KO (function)",   rule = "covered_fraction >= 0.5", M = ko_cf),
  list(layer = "ARG (resistome)", rule = "count >= 1",              M = arg_std),
  list(layer = "ARG (resistome)", rule = "covered_fraction >= 0.5", M = arg_cf))


for (l in c("KO (function)", "ARG (resistome)")) {
  s <- lapply(Filter(function(x) x$layer == l, MATS), function(x) sort(rownames(x$M)))
  if (!identical(s[[1]], s[[2]]))
    stop(sprintf("%s: presence rules cover different samples (%d vs %d)",
                 l, length(s[[1]]), length(s[[2]])))
}
for (m in MATS)
  logmsg("%-16s | %-24s : %d samples x %d features", m$layer, m$rule, nrow(m$M), ncol(m$M))

## ---- depth and richness by compartment -----------------------------------
diag_rows <- bind_rows(lapply(MATS, function(m) {
  kk  <- key[match(rownames(m$M), key$sample_id), ]
  rich <- rowSums(m$M > 0)
  bind_rows(lapply(COMP_LEVELS, function(cc) {
    i <- kk$compartment == cc
    if (sum(i) < 3) return(NULL)
    data.frame(layer = m$layer, presence_rule = m$rule, compartment = cc, n = sum(i),
               median_depth_reads = round(median(kk$posthost_reads[i])),
               median_richness    = round(median(rich[i])),
               rho_richness_depth = round(suppressWarnings(
                 cor(rich[i], kk$posthost_reads[i], method = "spearman")), 3))
  }))
}))
save_tab(diag_rows, "T37d_depth_richness_diagnostic.tsv")

# Analysed samples only: load_key() still carries P181 (dung, 25,052 reads), 
# which was excluded due to low deptj, so restrict to samples present in the matrices above.
key_dep <- key[key$sample_id %in% rownames(ko_std), ]
dep <- key_dep %>% group_by(compartment) %>%
  summarise(n = n(), median_depth = median(posthost_reads),
            min_depth = min(posthost_reads), max_depth = max(posthost_reads),
            .groups = "drop")
kw <- kruskal.test(posthost_reads ~ compartment, data = key_dep)
save_tab(transform(as.data.frame(dep),
                   kruskal_chisq = round(unname(kw$statistic), 3),
                   kruskal_df = unname(kw$parameter),
                   kruskal_p = signif(kw$p.value, 3)),
         "T37e_depth_by_compartment.tsv")
logmsg("")
logmsg("== depth by compartment ==")
for (i in seq_len(nrow(dep)))
  logmsg("  %-5s median depth %.1fM reads", dep$compartment[i], dep$median_depth[i] / 1e6)
logmsg("  Kruskal-Wallis depth ~ compartment: chi2 = %.1f, p = %.3g  [%s]",
       kw$statistic, kw$p.value,
       if (kw$p.value < 0.05) "differs by compartment" else "no compartment difference in depth")
logmsg("  richness-vs-depth within compartment (Spearman), count >= 1 rule:")
for (i in which(diag_rows$presence_rule == "count >= 1"))
  logmsg("    %-16s %-5s rho = %+.3f (median richness %d)", diag_rows$layer[i],
         diag_rows$compartment[i], diag_rows$rho_richness_depth[i], diag_rows$median_richness[i])

## ---- PERMANOVA + dispersion ----------------------------------------------
D <- lapply(MATS, function(m) lapply(DISTS, function(f) f(m$M)))

perm_rows <- list(); disp_rows <- list()
for (i in seq_along(MATS)) {
  m  <- MATS[[i]]
  kk <- key[match(rownames(m$M), key$sample_id), ]
  grp <- droplevels(factor(kk$compartment))
  for (dn in names(DISTS)) {
    d <- D[[i]][[dn]]
    a <- vegan::adonis2(d ~ grp, permutations = N_PERM, by = "terms")
    r2_null <- (nlevels(grp) - 1) / (length(grp) - 1)
    perm_rows[[length(perm_rows) + 1]] <- data.frame(
      layer = m$layer, presence_rule = m$rule, distance = dn, term = "compartment",
      n = length(grp), groups = nlevels(grp), n_features = ncol(m$M),
      R2 = round(a$R2[1], 4), R2_null_expected = round(r2_null, 4),
      R2_excess = round(a$R2[1] - r2_null, 4),
      F = round(a$F[1], 2), p = a$`Pr(>F)`[1])
    
    bd <- vegan::betadisper(d, grp)
    pt <- vegan::permutest(bd, permutations = N_PERM)
    disp_rows[[length(disp_rows) + 1]] <- data.frame(
      layer = m$layer, presence_rule = m$rule, distance = dn,
      F = round(pt$tab$F[1], 2), p = pt$tab$`Pr(>F)`[1],
      t(round(tapply(bd$distances, grp, mean), 3)))
  }
}
T37a <- bind_rows(perm_rows); T37b <- bind_rows(disp_rows)
save_tab(T37a, "T37a_unweighted_permanova.tsv")
save_tab(T37b, "T37b_unweighted_dispersion.tsv")

# Control: Bray-Curtis under the count rule should reproduce R2 from 03_diversity.R.
chk <- T37a %>% filter(distance == "Bray-Curtis (abundance)", presence_rule == "count >= 1")
pub <- c("KO (function)" = 0.856, "ARG (resistome)" = 0.681)
for (i in seq_len(nrow(chk))) {
  got <- chk$R2[i]; want <- pub[[chk$layer[i]]]
  logmsg("control: %s Bray-Curtis R2 %.4f vs T3 %.3f  [%s]", chk$layer[i], got, want,
         if (abs(got - want) < 0.002) "MATCH" else "MISMATCH")
  if (abs(got - want) >= 0.002)
    stop(sprintf("control failed: %s Bray-Curtis R2 %.4f != T3 %.3f", chk$layer[i], got, want))
}

## ---- Procrustes agreement between abundance and presence ordinations -----
proc_rows <- bind_rows(lapply(seq_along(MATS), function(i) {
  m <- MATS[[i]]
  a <- pcoa(D[[i]][["Bray-Curtis (abundance)"]], k = 3)$pts
  bind_rows(lapply(c("Jaccard (presence)", "Sorensen (presence)"), function(dn) {
    b  <- pcoa(D[[i]][[dn]], k = 3)$pts
    pr <- vegan::protest(a, b, permutations = N_PERM)
    data.frame(layer = m$layer, presence_rule = m$rule,
               comparison = sprintf("Bray-Curtis vs %s, PCoA (3 axes)", dn),
               procrustes_r = round(pr$t0, 3), p = pr$signif)
  }))
}))
save_tab(proc_rows, "T37c_unweighted_procrustes.tsv")

## ---- Figure 1: Bray-Curtis and Jaccard ordinations (count >= 1 def) ----------
panel_plot <- function(i, dn) {
  m  <- MATS[[i]]
  kk <- key[match(rownames(m$M), key$sample_id), ]
  p  <- pcoa(D[[i]][[dn]])
  df <- data.frame(compartment = kk$compartment, Axis1 = p$pts[, 1], Axis2 = p$pts[, 2])
  st <- T37a[T37a$layer == m$layer & T37a$presence_rule == m$rule & T37a$distance == dn, ]
  ggplot(df, aes(Axis1, Axis2, color = compartment)) +
    geom_point(size = 1.5, alpha = .75) +
    stat_ellipse(level = .95, linewidth = .4, show.legend = FALSE) +
    annotate("text", x = Inf, y = Inf, hjust = 1.04, vjust = 1.5, size = 3.5, color = "grey25",
             label = sprintf("R\u00b2 = %.3f, %s", st$R2[1], fmt_p(st$p[1]))) +
    scale_color_manual(values = COMP_COLS, name = NULL, limits = COMP_LEVELS,
                       labels = compartment_display) +
    scale_y_continuous(expand = expansion(mult = c(.05, .15))) +  
    # layer_display() maps the stored key to the reader-facing name
    labs(title = sprintf("%s \u00b7 %s", layer_display(m$layer, wrap = FALSE),
                         sub(" \\(.*\\)$", "", dn)), 
         x = sprintf("PCoA 1 (%.1f%%)", p$var[1]),
         y = sprintf("PCoA 2 (%.1f%%)", p$var[2])) +
    theme_sgma +
    theme(plot.title = element_text(face = "plain", size = 10))
}

f39 <- patchwork::wrap_plots(
  unlist(lapply(which(vapply(MATS, function(m) m$rule == "count >= 1", TRUE)), function(i)
    lapply(c("Bray-Curtis (abundance)", "Jaccard (presence)"),
           function(dn) panel_plot(i, dn))), recursive = FALSE),
  nrow = 2) +
  patchwork::plot_layout(guides = "collect")
save_fig(f39, "F39_unweighted_ordination.png", width = 9.5, height = 7.5)

## ---- report --------------------------------------------------------------
logmsg("")
logmsg("== PERMANOVA, compartment (R2_null_expected = %.3f) ==", T37a$R2_null_expected[1])
for (i in seq_len(nrow(T37a))) with(T37a[i, ], logmsg(
  "  %-16s | %-24s | %-24s R2 = %.4f  F = %6.1f  p = %.3f  (%d features)",
  layer, presence_rule, distance, R2, F, p, n_features))
logmsg("")
logmsg("== dispersion (betadisper) ==")
for (i in seq_len(nrow(T37b))) with(T37b[i, ], logmsg(
  "  %-16s | %-24s | %-24s F = %6.2f  p = %.3f", layer, presence_rule, distance, F, p))
logmsg("")
logmsg("== Procrustes ==")
for (i in seq_len(nrow(proc_rows))) with(proc_rows[i, ], logmsg(
  "  %-16s | %-24s | %-45s r = %.3f  p = %.3f", layer, presence_rule, comparison,
  procrustes_r, p))

writeLines(log_lines, file.path(LOG_DIR, "11_unweighted_distance.log"))
cat("\nUnweighted distance complete.\n")
