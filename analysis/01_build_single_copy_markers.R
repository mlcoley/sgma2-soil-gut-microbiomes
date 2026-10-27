#!/usr/bin/env Rscript
# Single-copy marker KOs for per-genome normalization.
# Bacterial ribosomal-protein KOs are near-universal and single copy, so their median
# abundance per sample estimates average genome copy number (cf. MUSiCC, Manor &
# Borenstein 2015). Source: KEGG BRITE ko03011 (Ribosome), Ribosomal proteins >
# Bacteria; the raw file is downloaded to feature_tables/kegg_ko03011_ribosome.keg
# if absent (retrieved 2026-07-22).
# Output: feature_tables/single_copy_marker_kos.tsv
# Usage: Rscript analysis/01_build_single_copy_markers.R

ft  <- file.path(Sys.getenv("SGMA2_DATA", "data"), "feature_tables")
keg <- file.path(ft, "kegg_ko03011_ribosome.keg")
out <- file.path(ft, "single_copy_marker_kos.tsv")
if (!file.exists(keg)) download.file("https://rest.kegg.jp/get/br:ko03011", keg, quiet = TRUE)
lines <- readLines(keg)

in_prot <- FALSE; in_bact <- FALSE; sub_now <- NA_character_
rows <- list()
for (ln in lines) {
  if (grepl("^ARibosomal proteins", ln)) in_prot <- TRUE
  else if (grepl("^ARibosomal RNAs", ln)) in_prot <- FALSE
  if (in_prot && grepl("^B  Bacteria", ln)) { in_bact <- TRUE; next }
  else if (in_prot && grepl("^B  ", ln)) in_bact <- FALSE
  if (in_bact && grepl("^C ", ln)) sub_now <- trimws(sub("^C\\s*", "", ln))
  if (in_bact && grepl("^D", ln)) {
    m <- regmatches(ln, regexpr("K[0-9]{5,}", ln))
    if (length(m) != 1) next
    # "D      K02863  RP-L1, MRPL1, rplA; large subunit ribosomal protein L1"
    rest <- trimws(sub(paste0("^D\\s*", m, "\\s*"), "", ln))
    sym  <- trimws(sub(";.*$", "", rest))
    desc <- if (grepl(";", rest)) trimws(sub("^[^;]*;\\s*", "", rest)) else NA_character_
    rows[[length(rows) + 1]] <- data.frame(
      ko = m, subunit = sub_now, gene_symbols = sym, protein = desc,
      stringsAsFactors = FALSE)
  }
}
mk <- do.call(rbind, rows)
mk <- mk[!duplicated(mk$ko), ]
mk <- mk[order(mk$ko), ]
mk$source   <- "KEGG BRITE ko03011 (Ribosome), Ribosomal proteins > Bacteria"
mk$accessed <- "2026-07-22"

stopifnot(nrow(mk) == 57, !any(is.na(mk$subunit)), !any(is.na(mk$protein)))
write.table(mk, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s : %d bacterial ribosomal-protein marker KOs (%s)\n",
            out, nrow(mk), paste(sprintf("%s %d", names(table(mk$subunit)),
                                         as.integer(table(mk$subunit))), collapse = ", ")))
