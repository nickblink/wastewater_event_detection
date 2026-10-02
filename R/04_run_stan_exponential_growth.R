# Bayesian exponential-growth outbreak detection: for each of the 10
# FIPS_SUBSET counties, fits the windowed exponential-growth Stan model
# (stan/modular_outbreak_detection.stan) separately to the WW series and the
# cases series, then detects outbreak start/end windows from the posterior
# growth-rate (eta) estimates.
#
# REQUIRES rstan (or cmdstanr) with a working C++ toolchain -- this is NOT
# installed by R/00_install_packages.R. Install with:
#   install.packages("rstan")
# and ensure you have a working compiler (see https://mc-stan.org/users/interfaces/).
#
# nsample/burnin are reduced from values used in the original paper so that
# fitting 10 counties x 2 series finishes in a reasonable time on a laptop;
# increase them for more reliable posterior estimates.

source("R/00_config.R")
source("R/stan_exponential_growth.R")
library(readr)
library(dplyr)

if (!requireNamespace("rstan", quietly = TRUE) && !requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("This script requires rstan or cmdstanr. See the header comment in this file for setup instructions.")
}

ww_cleaned <- read_csv(WW_CLEANED_FILE, show_col_types = FALSE)
cases_cleaned <- read_csv(CASES_CLEANED_FILE, show_col_types = FALSE)

NSAMPLE <- 500
BURNIN <- 250
WINDOW <- 5

fit_one_series <- function(df, out_col, family) {
  stan_data <- prep_stan_data(df, time_unit = "days", out_col = out_col, family = family, window = WINDOW)
  stan_lst <- run_stan_EWS(stan_data, nsample = NSAMPLE, burnin = BURNIN)
  organized <- organize_results_data(stan_summary = stan_lst$stan_summary, stan_data = stan_data, out_col = out_col)
  outbreaks <- detect_outbreaks(organized$df, organized$eta_est, outbreak_end = "downtrend")
  outbreaks$matched_outbreaks
}

all_outbreaks <- list()

for (fip in FIPS_SUBSET) {
  county_name <- FIPS_COUNTY_NAMES[[fip]]
  cat(sprintf("Fitting Stan model for %s (%s)...\n", fip, county_name))

  df_WW <- ww_cleaned %>% filter(fips == fip) %>% select(date, Y_interp)
  df_cases <- cases_cleaned %>% filter(fips == fip) %>% select(date, new_cases_smoothed)

  ww_outbreaks <- fit_one_series(df_WW, out_col = "Y_interp", family = "normal")
  cases_outbreaks <- fit_one_series(df_cases, out_col = "new_cases_smoothed", family = "negbin")

  all_outbreaks[[paste0(fip, "_WW")]] <- if (nrow(ww_outbreaks) > 0) {
    ww_outbreaks %>% mutate(fips = fip, county = county_name, series = "WW")
  } else NULL
  all_outbreaks[[paste0(fip, "_cases")]] <- if (nrow(cases_outbreaks) > 0) {
    cases_outbreaks %>% mutate(fips = fip, county = county_name, series = "cases")
  } else NULL
}

stan_results <- bind_rows(all_outbreaks) %>%
  select(fips, county, series, outbreak_start, outbreak_end)

cat("Detected outbreak windows:\n")
print(stan_results)

write_csv(stan_results, "results/stan_results.csv")
cat("Wrote results/stan_results.csv\n")
