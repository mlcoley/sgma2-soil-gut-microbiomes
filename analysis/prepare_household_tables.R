#!/usr/bin/env Rscript
# Builds household-level tables read by 05_exposure.R, 15-17 and the exposure panel from
# the long-format household questionnaire exports.
# Input directory $SGMA2_HOUSEHOLD (default data/household):
#   sgma2_fieldmetadata.csv  hq_uuid -> hh_id map (hq_uuid, hh_id; other columns, including
#                            coordinates, are not used)
#   sgma2_hq_module1_childhhdemographics.csv
#   long/sgma2_hq_{animalroster,animalconfinementpractices,childanimalcontact,
#                  childfecescontact,hhmembersroster}.csv  (each with hq_uuid)
# Outputs, written to $SGMA2_HOUSEHOLD: long_with_hh/*_with_hh.csv,
#   sgma2_hq_{animalroster,animalconfinementpractices,childanimalcontact,
#             childfecescontact}_wide.csv, sgma2_hq_module1_hhdem_agestructure.csv.
# The input tables are de-identified and archived at OSF; see inputs/README.md.
# Usage: Rscript analysis/prepare_household_tables.R   (from the repository root)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr); library(stringr)
})

meta_dir <- Sys.getenv("SGMA2_HOUSEHOLD",
                       file.path(Sys.getenv("SGMA2_DATA", "data"), "household"))
long_dir <- file.path(meta_dir, "long")
with_dir <- file.path(meta_dir, "long_with_hh")
dir.create(with_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. add household ids -------------------------------------------------
long_files <- c("sgma2_hq_animalroster.csv", "sgma2_hq_animalconfinementpractices.csv",
                "sgma2_hq_childanimalcontact.csv", "sgma2_hq_childfecescontact.csv",
                "sgma2_hq_hhmembersroster.csv")

ref <- read_csv(file.path(meta_dir, "sgma2_fieldmetadata.csv"), show_col_types = FALSE) %>%
  transmute(uuid = str_to_lower(str_trim(hq_uuid)), hh_id = hh_id)
stopifnot(!any(is.na(ref$uuid)))
if (any(duplicated(ref$uuid)))
  stop("field metadata has duplicated hq_uuid, e.g. ", ref$uuid[duplicated(ref$uuid)][1])

# harmonise Swahili animal names to English labels.
SWAHILI <- c("Ng'ombe" = "Cows", "Kuku" = "Chickens")
label_cols <- c("animal", "anim_confined", "anim_contact", "feces_contact")

add_household_id <- function(f) {
  dat <- read_csv(file.path(long_dir, f), show_col_types = FALSE)
  if (!"hq_uuid" %in% names(dat)) stop("hq_uuid not found in ", f)
  dat2 <- dat %>%
    mutate(uuid_join = str_to_lower(str_trim(hq_uuid))) %>%
    left_join(ref, by = c("uuid_join" = "uuid")) %>%
    mutate(hh_id = hh_id, .after = hq_uuid) %>%
    select(-uuid_join, -any_of(c("hq_index", "hq_uuid")))
  for (lc in intersect(label_cols, names(dat2)))
    dat2[[lc]] <- recode(dat2[[lc]], !!!SWAHILI)
  n_unmatched <- sum(is.na(dat2$hh_id))
  if (n_unmatched > 0) warning(f, ": ", n_unmatched, " row(s) with no household match")
  out <- file.path(with_dir, str_replace(f, "\\.csv$", "_with_hh.csv"))
  write_csv(dat2, out)
  out
}
invisible(map_chr(long_files, add_household_id))

# ---- 2. long to wide --------------------------------------------------------
rd_with_hh <- function(f) read_csv(file.path(with_dir, str_replace(f, "\\.csv$", "_with_hh.csv")),
                                   show_col_types = FALSE)

pivot_roster <- function(dat_long, name_col, out_file, rename_to = NULL) {
  dups <- dat_long %>% count(hh_id, .data[[name_col]]) %>% filter(n > 1)
  if (nrow(dups) > 0)
    stop("duplicate rows per household x ", name_col, "; aggregate before pivoting:\n",
         paste0(capture.output(print(head(dups, 10))), collapse = "\n"))
  wide <- dat_long %>%
    pivot_wider(id_cols = hh_id, names_from = all_of(name_col),
                values_from = -all_of(c("hh_id", name_col)), names_sep = "_") %>%
    arrange(hh_id)
  if (!is.null(rename_to)) names(wide) <- rename_to
  write_csv(wide, file.path(meta_dir, out_file))
}

pivot_roster(rd_with_hh("sgma2_hq_animalroster.csv"), "animal",
             "sgma2_hq_animalroster_wide.csv")
pivot_roster(rd_with_hh("sgma2_hq_animalconfinementpractices.csv"), "anim_confined",
             "sgma2_hq_animalconfinementpractices_wide.csv")
pivot_roster(rd_with_hh("sgma2_hq_childanimalcontact.csv"), "anim_contact",
             "sgma2_hq_childanimalcontact_wide.csv",
             rename_to = c("hh_id", "child_contact_cows", "child_contact_chickens",
                           "child_contact_sheep", "child_contact_goats"))
pivot_roster(rd_with_hh("sgma2_hq_childfecescontact.csv"), "feces_contact",
             "sgma2_hq_childfecescontact_wide.csv")

# ---- 3. household age structure -------------------------------------------
roster <- rd_with_hh("sgma2_hq_hhmembersroster.csv")
hq_dat <- read_csv(file.path(meta_dir, "sgma2_hq_module1_childhhdemographics.csv"),
                   show_col_types = FALSE)

hh_summary <- roster %>%
  mutate(
    age_yr = suppressWarnings(as.numeric(hhmember_age)),
    age_bin = case_when(
      age_yr < 5 ~ "u5",
      age_yr >= 5  & age_yr < 12 ~ "5_12",
      age_yr >= 12 & age_yr < 18 ~ "12_18",
      age_yr >= 18 & age_yr < 35 ~ "18_35",
      age_yr >= 35 & age_yr <= 70 ~ "35_70",
      age_yr > 70 ~ "o70",
      TRUE ~ NA_character_)
  ) %>%
  group_by(hh_id) %>%
  summarise(
    hh_members_total = n(),
    males_u5      = sum(hhmember_gender == "Male"   & age_bin == "u5",    na.rm = TRUE),
    females_u5    = sum(hhmember_gender == "Female" & age_bin == "u5",    na.rm = TRUE),
    males_5_12    = sum(hhmember_gender == "Male"   & age_bin == "5_12",  na.rm = TRUE),
    females_5_12  = sum(hhmember_gender == "Female" & age_bin == "5_12",  na.rm = TRUE),
    males_12_18   = sum(hhmember_gender == "Male"   & age_bin == "12_18", na.rm = TRUE),
    females_12_18 = sum(hhmember_gender == "Female" & age_bin == "12_18", na.rm = TRUE),
    males_18_35   = sum(hhmember_gender == "Male"   & age_bin == "18_35", na.rm = TRUE),
    females_18_35 = sum(hhmember_gender == "Female" & age_bin == "18_35", na.rm = TRUE),
    males_35_70   = sum(hhmember_gender == "Male"   & age_bin == "35_70", na.rm = TRUE),
    females_35_70 = sum(hhmember_gender == "Female" & age_bin == "35_70", na.rm = TRUE),
    males_o70     = sum(hhmember_gender == "Male"   & age_bin == "o70",   na.rm = TRUE),
    females_o70   = sum(hhmember_gender == "Female" & age_bin == "o70",   na.rm = TRUE),
    .groups = "drop")

if (nrow(anti_join(hh_summary, hq_dat, by = "hh_id")) > 0)
  warning("some households in the roster are not in the demographics table and are not merged")

hq_out <- hq_dat %>% left_join(hh_summary, by = "hh_id")
names(hq_out)[names(hq_out) == "household_total_members"] <- "members_reported"
write_csv(hq_out, file.path(meta_dir, "sgma2_hq_module1_hhdem_agestructure.csv"))
cat("household tables written to", meta_dir, "\n")
