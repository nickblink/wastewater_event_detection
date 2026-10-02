# One-time extraction script (NOT run by repo users / not part of the pipeline).
#
# Trims the full public JHU CSSE time series
# (time_series_covid19_confirmed_US.csv, wide format: one row per county, one
# column per date) down to just the 10 demo counties, and writes the result to
# data/raw/jhu_cases_subset.csv, which IS committed to this repo.
#
# Source for the full file: https://github.com/CSSEGISandData/COVID-19

source("R/00_config.R")
library(readr)
library(dplyr)

JHU_FULL_FILE <- "~/dev/personal/wastewater_EWS/data/raw/time_series_covid19_confirmed_US.csv"

raw <- read_csv(path.expand(JHU_FULL_FILE), show_col_types = FALSE)

raw$FIPS_padded <- sprintf("%05d", as.integer(raw$FIPS))

subset_df <- raw %>% filter(FIPS_padded %in% FIPS_SUBSET) %>% select(-FIPS_padded)

cat(sprintf("Rows (counties): %d\n", nrow(subset_df)))
stopifnot(nrow(subset_df) == length(FIPS_SUBSET))

write_csv(subset_df, JHU_RAW_FILE)
cat(sprintf("Wrote %s (%.2f MB)\n", JHU_RAW_FILE, file.size(JHU_RAW_FILE) / 1e6))
