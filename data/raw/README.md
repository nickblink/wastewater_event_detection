# `data/raw/`

## `nwss_county_subset.csv`

Extracted from the public CDC NWSS dataset ["CDC Wastewater Data for
SARS-CoV-2"](https://data.cdc.gov/Public-Health-Surveillance/CDC-Wastewater-Data-for-SARS-CoV-2/j9g8-acpt)
(dataset `j9g8-acpt`), which is site-level sample data (one row per
sample/site/date).

The full public download is ~300MB and is **not** committed to this repo. This
subset was produced by `data_prep_local_only/extract_nwss_subset.R`, which:

- keeps only `pcr_target == "sars-cov-2"`
- keeps only `pcr_target_units == "copies/l wastewater"` (drops dry-sludge and
  log10-transformed rows, which use different scales)
- drops rows with comma-separated multi-county `county_fips` (shared
  sewersheds serving more than one county)
- keeps samples collected before 2023-04-01, to overlap with the JHU cases
  data below (JHU stopped updating in March 2023)

Columns: `county_fips`, `sample_collect_date`, `site`, `population_served`,
`pcr_target_avg_conc` (SARS-CoV-2 concentration, copies/L wastewater).

This subset covers 793 counties (anyone can point `R/00_config.R`'s
`FIPS_SUBSET` at a different selection from this pool), not just the 10 used
in the demo pipeline.

## `jhu_cases_subset.csv`

Trimmed from Johns Hopkins CSSE's public
[`time_series_covid19_confirmed_US.csv`](https://github.com/CSSEGISandData/COVID-19)
(cumulative confirmed cases, wide format: one row per county, one column per
date), down to just the 10 demo counties in `R/00_config.R`. Produced by
`data_prep_local_only/extract_jhu_subset.R`.
