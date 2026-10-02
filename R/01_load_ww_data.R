# Loads the CDC NWSS wastewater subset (data/raw/nwss_county_subset.csv),
# aggregates multiple sites/samples per county-week down to one weekly value
# per county, fills in a complete weekly date grid, interpolates short gaps,
# and writes the cleaned series to data/processed/ww_cleaned.csv.
#
# This replaces the original wastewater_EWS repo's Biobot/MWRA-specific
# loaders (get_county_data / get_MWRA_data), which depended on private data
# files and hardcoded schemas that don't apply here.

source("R/00_config.R")
source("R/utils_data_cleaning.R")
library(readr)
library(dplyr)
library(lubridate)

raw <- read_csv(NWSS_RAW_FILE, show_col_types = FALSE) %>%
  mutate(
    fips = county_fips,
    sample_collect_date = as.Date(sample_collect_date)
  ) %>%
  filter(fips %in% FIPS_SUBSET,
         sample_collect_date >= START_DATE,
         sample_collect_date <= END_DATE)

stopifnot(nrow(raw) > 0)
missing_fips <- setdiff(FIPS_SUBSET, unique(raw$fips))
if (length(missing_fips) > 0) {
  stop(sprintf("No NWSS data found for FIPS: %s", paste(missing_fips, collapse = ", ")))
}

# A county can have multiple sites and multiple samples within the same week;
# collapse to a single weekly value per county (mean concentration).
weekly <- raw %>%
  mutate(date = floor_date(sample_collect_date, unit = "week")) %>%
  group_by(fips, date) %>%
  summarize(Y = mean(pcr_target_avg_conc, na.rm = TRUE), .groups = "drop")

cat("Weekly NWSS rows before interpolation, by county:\n")
print(weekly %>% count(fips))

# Fill in a complete weekly grid and interpolate short gaps, per county.
ww_cleaned <- weekly %>%
  group_by(fips) %>%
  group_modify(~ interp_wrapper(.x, interp_threshold = INTERP_THRESHOLD,
                                 outcome = "Y", cols_to_keep = character(0))) %>%
  ungroup() %>%
  select(fips, date, Y, Y_interp)

# QC checks
stopifnot(!any(duplicated(ww_cleaned[, c("fips", "date")])))
stopifnot(inherits(ww_cleaned$date, "Date"))

missingness <- ww_cleaned %>%
  group_by(fips) %>%
  summarize(n_weeks = n(), pct_missing = mean(is.na(Y)), .groups = "drop")
cat("WW missingness by county (fraction of weeks with no raw sample):\n")
print(missingness)

write_csv(ww_cleaned, WW_CLEANED_FILE)
cat(sprintf("Wrote %s (%d rows, %d counties)\n", WW_CLEANED_FILE, nrow(ww_cleaned), n_distinct(ww_cleaned$fips)))
