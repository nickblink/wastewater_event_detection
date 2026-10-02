# Loads the trimmed JHU cases subset (data/raw/jhu_cases_subset.csv, wide
# format: one row per county, one column per date), converts cumulative
# counts to adjusted weekly new-case counts, and writes the cleaned series to
# data/processed/cases_cleaned.csv.
#
# Cleaning steps (ported from the original wastewater_EWS repo's
# data_preprocessing.R):
#   1) cumulative -> daily diff
#   2) adjust_negative_cases_backwards: absorb negative diffs (data
#      corrections) into preceding days instead of reporting negative counts
#   3) aggregate to weekly totals
#   4) detect_and_average_zeros_balanced: redistribute post-zero-run
#      reporting spikes evenly across the zero run
#
# Difference from the original pipeline: weekly boundaries here use a
# self-contained lubridate::floor_date("week") convention (and must match the
# convention used in R/01_load_ww_data.R), rather than the original's
# non-portable dependency on a separately-loaded deaths file's week
# boundaries.

source("R/00_config.R")
source("R/utils_data_cleaning.R")
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

raw <- read_csv(JHU_RAW_FILE, show_col_types = FALSE)
raw$fips <- sprintf("%05d", as.integer(raw$FIPS))
stopifnot(setequal(raw$fips, FIPS_SUBSET))

date_cols <- names(raw)[!(names(raw) %in% c(
  "UID", "iso2", "iso3", "code3", "FIPS", "Admin2", "Province_State",
  "Country_Region", "Lat", "Long_", "Combined_Key", "fips"
))]

long <- raw %>%
  select(fips, all_of(date_cols)) %>%
  pivot_longer(cols = all_of(date_cols), names_to = "date_str", values_to = "cumulative_cases") %>%
  mutate(
    date = as.Date(date_str, format = "%m/%d/%y"),
    cumulative_cases = pmax(cumulative_cases, 0)
  ) %>%
  select(fips, date, cumulative_cases) %>%
  arrange(fips, date)

daily <- long %>%
  group_by(fips) %>%
  mutate(
    new_cases = cumulative_cases - lag(cumulative_cases, default = 0),
    adjusted_new_cases = adjust_negative_cases_backwards(new_cases)
  ) %>%
  ungroup()

weekly <- daily %>%
  mutate(date = floor_date(date, unit = "week")) %>%
  group_by(fips, date) %>%
  summarize(new_cases = sum(adjusted_new_cases, na.rm = TRUE), n_days = n(), .groups = "drop")

# Drop the first/last week for each county if it's a partial week (fewer than
# 7 days of daily data rolled into it), rather than relying on hardcoded
# boundary dates tied to this dataset's specific vintage.
weekly_complete <- weekly %>%
  group_by(fips) %>%
  filter(!(date == min(date) & n_days < 7),
         !(date == max(date) & n_days < 7)) %>%
  ungroup() %>%
  select(fips, date, new_cases)

cases_cleaned <- weekly_complete %>%
  filter(date >= START_DATE, date <= END_DATE) %>%
  group_by(fips) %>%
  arrange(date) %>%
  mutate(
    new_cases_unsmoothed = new_cases,
    new_cases_smoothed = detect_and_average_zeros_balanced(new_cases)
  ) %>%
  ungroup() %>%
  select(fips, date, new_cases_unsmoothed, new_cases_smoothed)

# QC checks
stopifnot(!any(duplicated(cases_cleaned[, c("fips", "date")])))
stopifnot(inherits(cases_cleaned$date, "Date"))
stopifnot(all(cases_cleaned$new_cases_smoothed >= 0))
missing_fips <- setdiff(FIPS_SUBSET, unique(cases_cleaned$fips))
if (length(missing_fips) > 0) {
  stop(sprintf("No cases data found for FIPS: %s", paste(missing_fips, collapse = ", ")))
}

cat("Weekly cases rows by county:\n")
print(cases_cleaned %>% count(fips))

write_csv(cases_cleaned, CASES_CLEANED_FILE)
cat(sprintf("Wrote %s (%d rows, %d counties)\n", CASES_CLEANED_FILE, nrow(cases_cleaned), n_distinct(cases_cleaned$fips)))
