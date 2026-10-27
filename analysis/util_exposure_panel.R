#!/usr/bin/env Rscript
# Household exposure and diet frames built from the household metadata files in
# $SGMA2_HOUSEHOLD (META_DIR). Source after 00_load_features.R.
#   exposure_panel()  one row per household, hh_id + 20 exposure variables
#   diet_panel()      one row per household, 24-hour recall (child) and FFQ (household)
# No coordinate is read.

exposure_rd <- function(f) {
  p <- file.path(META_DIR, f)
  if (!file.exists(p)) stop("metadata file missing: ", p)
  read.csv(p, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
           na.strings = c("NA", ""))
}

EXPOSURE_PRIMARY <- c("cows_head", "total_livestock_head", "cows_in_compound",
                      "child_touch_feces", "dung_to_crops",
                      "child_antibiotics_2mo", "vet_meds_cows")
EXPOSURE_EXPLORATORY <- c("chickens_head", "child_eats_soil", "dung_to_floor",
                          "dung_to_walls", "dung_disposed", "soap_available",
                          "water_untreated", "latrine_no_slab", "hh_leftover_abx")
# Added after the primary results were seen; own family, own correction.
EXPOSURE_POSTHOC_SES <- c("ses_income_crops", "ses_income_ownbusiness",
                          "ses_land_over_1acre", "ses_kitchen_garden")

exposure_family_of <- function(nm) {
  ifelse(nm %in% EXPOSURE_PRIMARY, "primary",
  ifelse(nm %in% EXPOSURE_POSTHOC_SES, "posthoc_SES", "exploratory"))
}

exposure_panel <- function() {
  suppressPackageStartupMessages({ library(dplyr) })

  husb  <- exposure_rd("sgma2_hq_module4_hhanimalhusbandry.csv")
  rost  <- exposure_rd("sgma2_hq_animalroster_wide.csv")
  conf  <- exposure_rd("sgma2_hq_animalconfinementpractices_wide.csv")
  wash  <- exposure_rd("sgma2_hq_module3_hhwash.csv")
  chlth <- exposure_rd("sgma2_hq_module5_childhealth.csv")
  ses   <- exposure_rd("sgma2_hq_module2_hhses.csv")

  yn   <- function(x) ifelse(is.na(x), NA, ifelse(x == "Yes", 1L, ifelse(x == "No", 0L, NA)))
  num  <- function(x) suppressWarnings(as.numeric(x))
  # Structural zeros: animal_head_Chickens is NA for the 10 households that keep no
  # chickens (a head count of 0), and the anim_confined_dung* block is NA for the 14
  # households whose cattle are not confined. dung_to_crops is constructed and screened
  # but excluded from modelling as not validly measured (see 05_exposure.R).
  zero_if_na <- function(x) ifelse(is.na(x), 0, x)

  exp_df <- data.frame(hh_id = husb$hh_id, stringsAsFactors = FALSE) %>%
    left_join(rost %>% transmute(hh_id,
                cows_head     = num(animal_head_Cows),
                chickens_head = zero_if_na(num(animal_head_Chickens)),
                vet_meds_cows = yn(animal_medications_Cows)), by = "hh_id") %>%
    left_join(husb %>% transmute(hh_id,
                child_touch_feces = yn(child_touchlivestockfeces),
                child_eats_soil   = yn(child_eatingsoil),
                # Q4.7, "Which of your livestock are kept in stables or pens inside the
                # compound?": cattle penned or stabled, not cattle presence (all 39
                # households keep cattle). It is also the survey routing variable for
                # the confinement module (25 asked / 14 skipped).
                cows_in_compound  = as.integer(grepl("Cows", livestock_compound))),
              by = "hh_id") %>%
    left_join(conf %>% transmute(hh_id,
                dung_to_crops = zero_if_na(num(anim_confined_dungcrops_Cows)),
                dung_to_floor = zero_if_na(num(anim_confined_dungfloor_Cows)),
                dung_to_walls = zero_if_na(num(anim_confined_dungwalls_Cows)),
                dung_disposed = zero_if_na(num(anim_confined_dungdisposed_Cows))),
              by = "hh_id") %>%
    left_join(wash %>% transmute(hh_id,
                soap_available  = yn(handwash_soapavail),
                water_untreated = as.integer(drinkingwater_treatment == "None"),
                latrine_no_slab = as.integer(grepl("without slab", toilet_type))),
              by = "hh_id") %>%
    left_join(chlth %>% transmute(hh_id,
                child_antibiotics_2mo = yn(child_2mo_antibiotics),
                hh_leftover_abx = ifelse(is.na(hh_antibiotics_leftoverusage), NA,
                                  as.integer(hh_antibiotics_leftoverusage == "Yes"))),
              by = "hh_id") %>%
    left_join(ses %>% transmute(hh_id,
                ses_income_crops       = num(hh_incomesources_crops),
                ses_income_ownbusiness = num(hh_incomesources_ownbusiness),
                # owning no land is 0 acres ("under 1 acre"); "No answer" stays NA
                ses_land_over_1acre    = ifelse(hh_landamount_own == "Yes - 1-5 acres", 1L,
                                         ifelse(hh_landamount_own %in%
                                                  c("Yes - less than 1 acre", "No"), 0L, NA)),
                ses_kitchen_garden     = num(hh_kitchengarden_incompound)),
              by = "hh_id")
  exp_df$total_livestock_head <- rowSums(
    cbind(exp_df$cows_head, exp_df$chickens_head), na.rm = TRUE)
  exp_df
}

## ===========================================================================
## Diet frame
## ===========================================================================
# These are dietary-composition covariates only. Food sourcing and raw/boiled
# preparation are not recorded, and the dairy item is "Cheese, yogurt, kefir, milk or
# other milk products" (sgma2_variableindex.csv), so they say nothing about
# animal-to-human gene transfer.

DIET_FFQ_FILE    <- "sgma2_hq_module6_hhffq.csv"      # household, habitual, ordinal
DIET_RECALL_FILE <- "sgma2_hq_module7_24hrrecall.csv" # caregiver + child, single day, counts

# Ordinal coding of the four FFQ levels; any other level is an error.
DIET_FFQ_LEVELS <- c("Never" = 0L, "Rarely" = 1L, "Sometimes" = 2L, "Every day" = 3L)

# FAO/WHO MDD-W food groups mapped onto the instrument's 28 items. The instrument
# merges MDD-W groups 2 and 3 into one item ("Beans, peas, lentils, nuts or seeds"),
# so 9 groups are achievable. `cowpea` (cowpea leaves) is assigned to dark green leafy
# vegetables. Oils, sweets, soda, tea and coffee are not MDD-W groups and are excluded.
DIET_FOOD_GROUPS <- list(
  starchy_staples  = c("ugali", "tubersgourds_white"),
  pulses_nuts      = c("legumes"),
  dairy            = c("milk"),
  meat_fish        = c("organmeats", "meats", "fish", "snails"),
  eggs             = c("eggs"),
  dark_green_leafy = c("spinach", "amaranth", "nightshade", "spiderplant",
                       "cowpea", "jutemallow", "kale", "sweetpotatoleaves"),
  vitA_rich        = c("tubersgourds_orange", "ripemango"),
  other_veg        = c("okra", "africaneggplant"),
  other_fruit      = c("banana", "orange"))

# Animal-source set for the FFQ composite (`snails` fails the degeneracy screen).
DIET_ASF_ITEMS <- c("organmeats", "meats", "fish", "eggs", "milk")

# Recall column names are child_<item> except for organmeats.
DIET_RECALL_COL <- function(item) {
  ifelse(item == "organmeats", "child_organ_meats", paste0("child_", item))
}

# One row per household, both instruments. Caregiver recall columns are not returned.
diet_panel <- function(verify = TRUE) {
  suppressPackageStartupMessages({ library(dplyr) })

  ffq <- exposure_rd(DIET_FFQ_FILE)
  rec <- exposure_rd(DIET_RECALL_FILE)
  all_items <- unlist(DIET_FOOD_GROUPS, use.names = FALSE)

  ffq_ord <- sapply(all_items, function(it) {
    v <- ffq[[paste0("ffq_", it)]]
    if (is.null(v)) stop("FFQ item missing: ffq_", it)
    bad <- setdiff(unique(v[!is.na(v)]), names(DIET_FFQ_LEVELS))
    if (length(bad)) stop("unexpected FFQ level(s) in ffq_", it, ": ",
                          paste(bad, collapse = ", "))
    unname(DIET_FFQ_LEVELS[v])
  })
  ffq_ord <- as.data.frame(ffq_ord)

  rec_cnt <- sapply(all_items, function(it) {
    cn <- DIET_RECALL_COL(it)
    v <- rec[[cn]]
    if (is.null(v)) stop("recall column missing: ", cn)
    suppressWarnings(as.numeric(v))
  })
  rec_cnt <- as.data.frame(rec_cnt)

  # diversity score: number of MDD-W groups with any consumption
  dds <- function(df) {
    rowSums(sapply(DIET_FOOD_GROUPS, function(items) {
      as.integer(rowSums(df[, items, drop = FALSE] > 0, na.rm = TRUE) > 0)
    }))
  }

  out <- data.frame(hh_id = ffq$hh_id, stringsAsFactors = FALSE)

  # FFQ, household, habitual
  out$ffq_asf_score <- rowSums(ffq_ord[, DIET_ASF_ITEMS, drop = FALSE])
  out$ffq_dds       <- dds(ffq_ord)
  out$ffq_dairy     <- ffq_ord$milk
  # The dairy item has three realised levels (no household answers "Never"); it is
  # dichotomised at "Sometimes or more often", fixed on the instrument's wording.
  out$ffq_dairy_regular <- as.integer(ffq_ord$milk >= 2)

  # 24-hour recall, child. Any intake yesterday is the primary form for the two
  # items; the count is carried alongside.
  stopifnot(identical(as.character(rec$hh_id), as.character(ffq$hh_id)))
  out$child_dairy_any   <- as.integer(rec_cnt$milk > 0)
  out$child_fish_any    <- as.integer(rec_cnt$fish > 0)
  out$child_dairy_count <- rec_cnt$milk
  out$child_dds         <- dds(rec_cnt)

  if (verify) diet_panel_verify(out)
  out
}

# Assert the diet frame against counts measured from the raw files.
diet_panel_verify <- function(out) {
  bad <- character(0)
  if (nrow(out) != 39) bad <- c(bad, sprintf("expected 39 households, got %d", nrow(out)))
  if (anyDuplicated(out$hh_id)) bad <- c(bad, "duplicate hh_id")
  exp_counts <- c(child_dairy_any = 17L, child_fish_any = 13L)
  for (nm in names(exp_counts)) {
    got <- sum(out[[nm]], na.rm = TRUE)
    if (got != exp_counts[[nm]])
      bad <- c(bad, sprintf("%s: %d households with intake, expected %d",
                            nm, got, exp_counts[[nm]]))
  }
  if (max(out$ffq_dds, na.rm = TRUE) > length(DIET_FOOD_GROUPS) ||
      max(out$child_dds, na.rm = TRUE) > length(DIET_FOOD_GROUPS))
    bad <- c(bad, "a diversity score exceeds the number of food groups")
  if (length(bad))
    stop("diet frame failed verification against its sources:\n  ",
         paste(bad, collapse = "\n  "))
  invisible(TRUE)
}
