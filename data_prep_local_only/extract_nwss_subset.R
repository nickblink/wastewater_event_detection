# One-time extraction script (NOT run by repo users / not part of the pipeline).
#
# Filters the full public CDC NWSS "Wastewater Data for SARS-CoV-2" download
# (https://data.cdc.gov/Public-Health-Surveillance/CDC-Wastewater-Data-for-SARS-CoV-2/j9g8-acpt)
# down to a small subset that is committed to this repo at data/raw/nwss_county_subset.csv.
#
# The full source file (~300MB) is never committed. Download it yourself from the
# link above and point NWSS_FULL_FILE at it before running this script.

library(readr)
library(dplyr)

NWSS_FULL_FILE <- "~/Downloads/CDC_Wastewater_Data_for_SARS-CoV-2_20261002.csv"
OUT_FILE <- "data/raw/nwss_county_subset.csv"
CUTOFF_DATE <- as.Date("2023-04-01") # exclusive upper bound; keeps overlap with JHU cases data

raw <- read_csv(
  path.expand(NWSS_FULL_FILE),
  col_select = c(county_fips, sample_collect_date, site, population_served,
                 pcr_target, pcr_target_avg_conc, pcr_target_units),
  show_col_types = FALSE
)

subset_df <- raw %>%
  filter(
    pcr_target == "sars-cov-2",
    pcr_target_units == "copies/l wastewater",   # drop dry-sludge and log10 variants for a consistent scale
    !grepl(",", county_fips),                      # drop shared/multi-county sewersheds
    sample_collect_date < CUTOFF_DATE
  ) %>%
  select(county_fips, sample_collect_date, site, population_served, pcr_target_avg_conc) %>%
  arrange(county_fips, sample_collect_date, site)

cat(sprintf("Rows: %d, distinct counties: %d\n", nrow(subset_df), n_distinct(subset_df$county_fips)))

write_csv(subset_df, OUT_FILE)
cat(sprintf("Wrote %s (%.1f MB)\n", OUT_FILE, file.size(OUT_FILE) / 1e6))
