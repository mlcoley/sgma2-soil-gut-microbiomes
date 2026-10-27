#!/usr/bin/env Rscript
# Word-table and CSV machinery for 25_manuscript_tables.R: fonts, borders, header
# typography, captions, the .docx and .csv writers, and a vocabulary guard that stops
# before writing if banned wording would reach a reader-facing file. Defines functions
# only.
# Conventions: fmt_p() never prints "<0.001" (every permutation p has a floor at
# 1/(n_perm + 1)); fmt_p_tail() is for asymptotic p-values (Kruskal-Wallis, Wilcoxon);
# adjusted p-values are labelled "BH-adjusted p"; layer and compartment names go through
# layer_display() and compartment_display() so tables and figures print the same words.

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

FONT_DOCX <- "Arial"
DOCX_DIR  <- file.path(OUT_DIR, "tables_docx")
CSV_DIR   <- file.path(OUT_DIR, "supplementary_data")
dir.create(DOCX_DIR, recursive = TRUE, showWarnings = FALSE)

# --- formatting -------------------------------------------------------------

EMDASH <- "—"

fmt_num <- function(x, d = 3) ifelse(is.na(x), EMDASH, sprintf(paste0("%.", d, "f"), x))

fmt_int <- function(x) ifelse(is.na(x), EMDASH, formatC(x, format = "d", big.mark = ","))

fmt_signed <- function(x, d = 4) ifelse(is.na(x), EMDASH,
                                        sprintf(paste0("%+.", d, "f"), x))

fmt_pm <- function(m, s, d = 1) {
  f <- function(x) formatC(x, format = "f", digits = d, big.mark = ",")
  ifelse(is.na(m), EMDASH, sprintf("%s \u00b1 %s", f(m), f(s)))
}

fmt_p <- function(p, floor_at = NULL) {
  if (!is.null(floor_at) && any(!is.na(p) & p < floor_at - 1e-12))
    stop("p below the stated permutation floor of ", floor_at,
         " -- check n_perm before rendering.", call. = FALSE)
  ifelse(is.na(p), EMDASH, sprintf("%.3f", p))
}

fmt_p_tail <- function(p, dp = 4) {
  floor_p <- 10^-dp
  ifelse(is.na(p), EMDASH,
         ifelse(p < floor_p,
                sprintf(paste0("< %.", dp, "f"), floor_p),
                sub("\\.$", "", sub("0+$", "", sprintf(paste0("%.", dp, "f"), p)))))
}

TERM_DISPLAY <- c(
  "household (within soil)"      = "household (soils only)",
  "compartment (soil collapsed)" = "compartment (soil mean composition)")

term_display <- function(x) {
  out <- unname(TERM_DISPLAY[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

# Presence rules are stored as the code expression that defines them. A table
PRESENCE_RULE_DISPLAY <- c(
  "count >= 1"              = "≥ 1 read",
  "count >= 5"              = "≥ 5 reads",
  "covered_fraction >= 0.5" = "gene covered ≥ 50%")

rule_display <- function(x) {
  out <- unname(PRESENCE_RULE_DISPLAY[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

CONTRAST_DISPLAY <- c("compartment (omnibus)" = "All three compartments")

contrast_display <- function(x) {
  x <- as.character(x)
  mapped <- unname(CONTRAST_DISPLAY[x])
  out <- vapply(seq_along(x), function(i) {
    if (!is.na(mapped[i])) return(mapped[i])
    parts <- strsplit(x[i], " vs ", fixed = TRUE)[[1]]
    paste(compartment_display(parts), collapse = " vs ")
  }, character(1))
  out
}

METRIC_DISPLAY <- c(richness = "Observed richness", shannon = "Shannon")

metric_display <- function(x) {
  out <- unname(METRIC_DISPLAY[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

pair_display <- function(x) {
  x <- as.character(x)
  stopifnot(!any(grepl("-", names(COMPARTMENT_DISPLAY), fixed = TRUE)))
  vapply(x, function(p) {
    parts <- strsplit(p, "-", fixed = TRUE)[[1]]
    paste(compartment_display(parts), collapse = "\u2013")   # en dash
  }, character(1), USE.NAMES = FALSE)
}

STATISTIC_DISPLAY <- c(
  "Bray-Curtis similarity"               = "Bray\u2013Curtis similarity",
  "shared richness"                      = "Shared variant richness",
  "rare co-presence (prev<=0.25)"        = "Rare-variant co-presence (prevalence \u2264 25%)",
  "rare Jaccard (prev<=0.25)"            = "Rare-variant Jaccard (prevalence \u2264 25%)",
  "rare co-presence (prev<=0.10)"        = "Rare-variant co-presence (prevalence \u2264 10%)",
  "rare Jaccard (prev<=0.10)"            = "Rare-variant Jaccard (prevalence \u2264 10%)",
  "joint rare richness (confound probe)" = "Joint rare-variant richness (confound probe)")

statistic_display <- function(x) .map_display(x, STATISTIC_DISPLAY)

FAMILY_DISPLAY <- c(primary = "Primary", exploratory = "Exploratory",
                    posthoc_SES = "Post hoc (socioeconomic)")

COVARIATE_DISPLAY <- c(
  child_age_months = "child age (months)",
  child_weight_kg  = "child weight (kg)",
  child_height_cm  = "child height (cm)",
  child_female     = "child sex (1 = female)",
  bsc_sample       = "Bristol Stool Chart, sample collected",
  bsc_typical      = "Bristol Stool Chart, child's usual stool")

DIET_DISPLAY <- c(
  child_dairy_any   = "any milk products yesterday (child)",
  child_fish_any    = "any fish yesterday (child)",
  child_dds         = "dietary diversity, food groups (child, 24 h)",
  ffq_asf_score     = "animal-source food score (household, habitual)",
  ffq_dairy_regular = "milk products at least 'sometimes' (household)")

DIET_FAMILY_DISPLAY <- c(
  A = "A: child 24-hour recall",
  B = "B: household food frequency")

DIET_ARM_DISPLAY <- c(
  scalar        = "Scalar response",
  compositional = "Gut composition")

covariate_display <- function(x) .map_display(x, COVARIATE_DISPLAY)
diet_display <- function(x) .map_display(x, DIET_DISPLAY)
diet_family_display <- function(x) .map_display(x, DIET_FAMILY_DISPLAY)
diet_arm_display <- function(x) .map_display(x, DIET_ARM_DISPLAY)
family_display   <- function(x) .map_display(x, FAMILY_DISPLAY)

BANNED_HARD <- c("gating", "gated", "ungated", "blocklist", "blocked",
                 "private", "privacy", "repertoire")
BANNED_SOFT <- c("drivers", "q-value", "q_BH", "qBH")

assert_reader_safe <- function(spec) {
  # A lettered table holds its cells one level down, and a guard that silently
  # checks nothing is worse than no guard. Recurse before collecting text.
  if (!is.null(spec$parts)) {
    invisible(lapply(spec$parts, assert_reader_safe))
    spec <- list(title = spec$title, df = data.frame(), caption = spec$caption)
  }
  if (!is.null(spec$sheets)) {
    invisible(lapply(spec$sheets, assert_reader_safe))
    spec <- list(title = spec$title, df = data.frame(), caption = spec$caption)
  }
  txt <- c(spec$title, names(spec$df), spec$spanner,
           unlist(lapply(spec$df, as.character)), spec$caption,
           names(spec$definitions), unname(spec$definitions))
  txt <- txt[!is.na(txt)]
  hit <- function(words) {
    re <- paste0("\\b(", paste(words, collapse = "|"), ")\\b")
    unique(txt[grepl(re, txt, ignore.case = TRUE, perl = TRUE)])
  }
  h <- hit(BANNED_HARD)
  if (length(h))
    stop("Banned vocabulary would reach a .docx:\n  ",
         paste(h, collapse = "\n  "), call. = FALSE)
  s <- hit(BANNED_SOFT)
  if (length(s))
    warning("Check reader-facing wording in '", spec$title, "':\n  ",
            paste(s, collapse = "\n  "), call. = FALSE)
  invisible(TRUE)
}

tbl_spec <- function(title, df, caption = character(0), indent = logical(0),
                     spanner = character(0), definitions = NULL) {
  if (length(spanner) && length(spanner) != ncol(df))
    stop("spanner must have one entry per column (", ncol(df), ")", call. = FALSE)
  list(title = title, df = df, caption = caption, indent = indent,
       spanner = spanner, definitions = definitions)
}

tbl_workbook <- function(title, sheets, caption = character(0)) {
  if (is.null(names(sheets)) || any(!nzchar(names(sheets))))
    stop("every sheet needs a name", call. = FALSE)
  if (any(nchar(names(sheets)) > 31))
    stop("Excel sheet names are capped at 31 characters", call. = FALSE)
  list(title = title, sheets = sheets, caption = caption)
}

tbl_parts <- function(title, parts, caption = character(0))
  list(title = title, parts = parts, caption = caption)

# --- text runs --------------------------------------------------------------

.txt   <- function(x, size = 9, bold = FALSE, italic = FALSE)
  ftext(x, fp_text(font.family = FONT_DOCX, font.size = size, bold = bold, italic = italic))
.title_fpar <- function(x) fpar(.txt(x, size = 11, bold = TRUE))
.part_fpar  <- function(x) fpar(.txt(x, size = 10, bold = TRUE))
.fn    <- function(x) .txt(x, size = 9)
.fni   <- function(x) .txt(x, size = 9, italic = TRUE)
.hdr   <- function(x) as_chunk(x, props = fp_text(bold = TRUE, font.family = FONT_DOCX,
                                                  font.size = 10))
.hdri  <- function(x) as_chunk(x, props = fp_text(bold = TRUE, italic = TRUE,
                                                  font.family = FONT_DOCX, font.size = 10))
.hdrsup <- function(x) as_chunk(x, props = fp_text(bold = TRUE, font.family = FONT_DOCX,
                                                   font.size = 10,
                                                   vertical.align = "superscript"))

#' Footnotes are plain strings; *stars* mark italic runs.
fn_fpar <- function(s) {
  parts <- strsplit(s, "\\*")[[1]]
  chunks <- lapply(seq_along(parts), function(i)
    if (i %% 2 == 0) .fni(parts[i]) else .fn(parts[i]))
  do.call(fpar, chunks)
}

hdr_paragraph <- function(x) {
  # exact matches first
  if (x == "p") return(as_paragraph(.hdri("p")))
  if (x == "F") return(as_paragraph(.hdri("F")))
  if (x == "n") return(as_paragraph(.hdri("n")))
  # "BH-adjusted p" and friends: italicise only the trailing symbol
  m <- regmatches(x, regexec("^(.*\\s)(p|F|n)$", x))[[1]]
  if (length(m) == 3) {
    lead <- hdr_chunks_r2(m[2])
    return(do.call(as_paragraph, c(lead, list(.hdri(m[3])))))
  }
  do.call(as_paragraph, hdr_chunks_r2(x))
}

hdr_chunks_r2 <- function(x) {
  x <- gsub("R²", "R2", x)
  if (!grepl("R2", x, fixed = TRUE)) return(list(.hdr(x)))
  parts <- strsplit(x, "R2", fixed = TRUE)[[1]]
  if (length(parts) == 0) parts <- c("", "")
  if (length(parts) == 1) parts <- c(parts, "")
  out <- list()
  for (i in seq_along(parts)) {
    if (nzchar(parts[i])) out <- c(out, list(.hdr(parts[i])))
    if (i < length(parts)) out <- c(out, list(.hdr("R"), .hdrsup("2")))
  }
  out
}

sty <- function(ft) {
  ft |>
    font(fontname = FONT_DOCX, part = "all") |>
    fontsize(size = 10, part = "all") |>
    bold(part = "header") |>
    border_remove() |>
    hline_top(part = "header", border = fp_border(width = 1.5)) |>
    hline_bottom(part = "header", border = fp_border(width = 1.0)) |>
    hline_bottom(part = "body", border = fp_border(width = 1.5)) |>
    padding(padding.top = 2, padding.bottom = 2, part = "all") |>
    line_spacing(space = 1, part = "all")
}

build_ft <- function(t) {
  df <- t$df
  hdr <- names(df)
  names(df) <- paste0("c", seq_along(df))
  df[] <- lapply(df, as.character)
  df[is.na(df)] <- ""

  ft <- flextable(df)
  for (j in seq_along(hdr))
    ft <- compose(ft, j = j, part = "header", value = hdr_paragraph(hdr[j]))

  ft <- sty(ft) |>
    align(part = "header", align = "center") |>
    align(j = 1, part = "header", align = "left") |>
    align(part = "body", align = "center") |>
    align(j = 1, part = "body", align = "left") |>
    valign(valign = "bottom", part = "header")

  if (length(t$spanner)) {
    sp  <- t$spanner
    run <- rle(sp)
    ft  <- add_header_row(ft, values = run$values, colwidths = run$lengths, top = TRUE)
    ft  <- bold(ft, i = 1, part = "header") |>
      fontsize(i = 1, size = 10, part = "header") |>
      align(i = 1, align = "center", part = "header")
    ends <- cumsum(run$lengths)
    for (k in which(nzchar(run$values))) {
      cols <- (ends[k] - run$lengths[k] + 1):ends[k]
      ft <- hline(ft, i = 1, j = cols, part = "header",
                  border = fp_border(width = 0.75))
    }
  }

  if (length(t$indent) && any(t$indent))
    ft <- padding(ft, i = which(t$indent), j = 1, padding.left = 14, part = "body")

  set_table_properties(ft, layout = "autofit", width = 1)
}

# --- writers ----------------------------------------------------------------

#' One part: heading, its own caption, table. `heading` chooses the weight,
#' so the same routine sets a standalone table and a lettered part.
.add_block <- function(doc, t, heading = .title_fpar) {
  doc <- body_add_fpar(doc, heading(t$title))
  for (cp in t$caption) doc <- body_add_fpar(doc, fn_fpar(cp))
  doc |>
    body_add_par("", style = "Normal") |>
    body_add_flextable(build_ft(t)) |>
    body_add_par("", style = "Normal")
}

.add_table <- function(doc, t) {
  if (is.null(t$parts)) return(.add_block(doc, t))

  # Lettered table: one number, one heading, the shared caption, then each part
  # under its own letter.
  doc <- body_add_fpar(doc, .title_fpar(t$title))
  for (cp in t$caption) doc <- body_add_fpar(doc, fn_fpar(cp))
  for (i in seq_along(t$parts)) {
    doc <- body_add_par(doc, "", style = "Normal")
    doc <- .add_block(doc, t$parts[[i]], heading = .part_fpar)
  }
  doc
}

#' One .docx per table. Landscape, because these tables are wider than the
#' COMETS originals -- a 9-column PERMANOVA table does not fit portrait at 10 pt.
write_table_docx <- function(t, filename, landscape = TRUE) {
  assert_reader_safe(t)
  sec <- if (landscape)
    prop_section(page_size = page_size(orient = "landscape"),
                 type = "continuous",
                 page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1))
  else NULL
  doc <- .add_table(read_docx(), t)
  if (!is.null(sec)) doc <- body_end_block_section(doc, block_section(sec))
  print(doc, target = file.path(DOCX_DIR, filename))
  message("wrote ", filename)
  invisible(file.path(DOCX_DIR, filename))
}

XLSX_SHEET_DATA <- "Data"
XLSX_SHEET_DEFS <- "Definitions"

#' One data worksheet, styled the same way whether it is the only sheet or one
#' of several.
.add_xlsx_sheet <- function(wb, name, df) {
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df, headerStyle = openxlsx::createStyle(
    textDecoration = "bold", fontName = FONT_DOCX, fontSize = 10,
    border = "bottom", wrapText = TRUE, valign = "bottom"))
  openxlsx::freezePane(wb, name, firstRow = TRUE)
  openxlsx::addFilter(wb, name, row = 1, cols = seq_along(df))
  openxlsx::setColWidths(wb, name, cols = seq_along(df), widths = "auto")
  invisible(wb)
}

write_table_xlsx <- function(t, filename) {
  assert_reader_safe(t)

  # One spec and a multi-sheet workbook are the same object with one sheet, so
  # normalise here rather than branching through the whole function.
  sheets <- if (!is.null(t$sheets)) t$sheets else
    stats::setNames(list(t), XLSX_SHEET_DATA)

  for (nm in names(sheets)) {
    sh <- sheets[[nm]]
    if (is.null(sh$definitions))
      stop("sheet '", nm, "' has no `definitions`: a spreadsheet has no caption.",
           call. = FALSE)
    miss <- setdiff(names(sh$df), names(sh$definitions))
    if (length(miss))
      stop("undefined column(s) on sheet '", nm, "': ",
           paste(miss, collapse = ", "), call. = FALSE)
  }

  wb <- openxlsx::createWorkbook()
  for (nm in names(sheets)) .add_xlsx_sheet(wb, nm, sheets[[nm]]$df)

  # Glossary: every sheet's columns, under its own heading when there is more
  # than one, so a definition can never be attached to the wrong sheet.
  defs <- do.call(rbind, lapply(names(sheets), function(nm) {
    d <- sheets[[nm]]$definitions
    body <- data.frame(Column = names(d), Definition = unname(d),
                       check.names = FALSE, stringsAsFactors = FALSE)
    if (length(sheets) == 1) return(body)
    rbind(data.frame(Column = paste0("[", nm, "]"),
                     Definition = sheets[[nm]]$title,
                     check.names = FALSE, stringsAsFactors = FALSE),
          body,
          data.frame(Column = "", Definition = "",
                     check.names = FALSE, stringsAsFactors = FALSE))
  }))
  openxlsx::addWorksheet(wb, XLSX_SHEET_DEFS)
  openxlsx::writeData(wb, XLSX_SHEET_DEFS, t$title, startRow = 1)
  openxlsx::addStyle(wb, XLSX_SHEET_DEFS, openxlsx::createStyle(
    textDecoration = "bold", fontName = FONT_DOCX, fontSize = 11), rows = 1, cols = 1)
  openxlsx::writeData(wb, XLSX_SHEET_DEFS, defs, startRow = 3,
                      headerStyle = openxlsx::createStyle(
                        textDecoration = "bold", fontName = FONT_DOCX, fontSize = 10))
  if (length(t$caption)) {
    # *stars* mark italics in the .docx path; a spreadsheet cell has no runs, so
    # they are stripped before writing.
    notes <- gsub("\\*", "", t$caption)
    openxlsx::writeData(wb, XLSX_SHEET_DEFS, data.frame(Notes = notes),
                        startRow = nrow(defs) + 5,
                        headerStyle = openxlsx::createStyle(
                          textDecoration = "bold", fontName = FONT_DOCX, fontSize = 10))
  }
  openxlsx::setColWidths(wb, XLSX_SHEET_DEFS, cols = 1:2, widths = c(46, 110))

  openxlsx::saveWorkbook(wb, file.path(DOCX_DIR, filename), overwrite = TRUE)
  message("wrote ", filename, " (",
          paste(sprintf("%s: %d rows", names(sheets),
                        vapply(sheets, function(x) nrow(x$df), integer(1))),
                collapse = "; "), ")")
  invisible(file.path(DOCX_DIR, filename))
}

write_table_csv <- function(t, basename) {
  assert_reader_safe(t)
  if (!is.null(t$sheets))
    stop("write_table_csv() takes a single tbl_spec, not a workbook: split it ",
         "into one spec per file.", call. = FALSE)
  if (is.null(t$definitions))
    stop("'", basename, "' has no `definitions`: a data file has no caption.",
         call. = FALSE)
  miss <- setdiff(names(t$df), names(t$definitions))
  if (length(miss))
    stop("undefined column(s) in '", basename, "': ",
         paste(miss, collapse = ", "), call. = FALSE)

  dir.create(CSV_DIR, recursive = TRUE, showWarnings = FALSE)
  data_file <- file.path(CSV_DIR, paste0(basename, ".csv"))
  defs_file <- file.path(CSV_DIR, paste0(basename, "_definitions.csv"))

  utils::write.csv(t$df, data_file, row.names = FALSE, na = "NA")

  strip <- function(x) gsub("\\*", "", x)
  defs <- rbind(
    data.frame(Column = "[Title]", Definition = strip(t$title),
               check.names = FALSE, stringsAsFactors = FALSE),
    data.frame(Column = "[Data file]", Definition = basename(data_file),
               check.names = FALSE, stringsAsFactors = FALSE),
    data.frame(Column = "[Missing values]",
               Definition = paste0("Written as NA. A blank cell would be ",
                                   "indistinguishable from a failed write."),
               check.names = FALSE, stringsAsFactors = FALSE),
    data.frame(Column = names(t$definitions), Definition = strip(unname(t$definitions)),
               check.names = FALSE, stringsAsFactors = FALSE))
  if (length(t$caption))
    defs <- rbind(defs, data.frame(
      Column = sprintf("[Caption %d]", seq_along(t$caption)),
      Definition = strip(t$caption), check.names = FALSE, stringsAsFactors = FALSE))
  utils::write.csv(defs, defs_file, row.names = FALSE, na = "NA")

  message("wrote ", basename(data_file), " (", nrow(t$df), " rows, ",
          ncol(t$df), " columns) + ", basename(defs_file))
  invisible(c(data_file, defs_file))
}

#' Combined document, one table per page, same content and styling.
write_tables_docx <- function(tables, filename = "all_tables.docx", landscape = TRUE) {
  invisible(lapply(tables, assert_reader_safe))
  doc <- read_docx()
  for (i in seq_along(tables)) {
    doc <- .add_table(doc, tables[[i]])
    if (i < length(tables)) doc <- body_add_break(doc)
  }
  if (landscape)
    doc <- body_end_block_section(doc, block_section(prop_section(
      page_size = page_size(orient = "landscape"), type = "continuous",
      page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1))))
  print(doc, target = file.path(DOCX_DIR, filename))
  message("wrote ", filename)
  invisible(file.path(DOCX_DIR, filename))
}
