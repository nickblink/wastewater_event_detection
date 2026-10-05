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
# nsample/burnin are reduced from the original paper's defaults (nsample =
# 10000, burnin = 1000) so that fitting 10 counties x 2 series finishes in a
# reasonable time on a laptop; increase them for more reliable posterior
# estimates (higher nsample noticeably reduces the divergence rate here).
#
# time_unit = "week": our data (like the original paper's) is weekly-cadence,
# so T is indexed in integer weeks and window = 5 means a 5-*week* lookback
# (~5 data points per fitting window). Passing time_unit = "days" here was an
# earlier bug in this port -- it made T count individual days while samples
# are 7 days apart, so most fitting windows contained only one data point and
# the model was essentially unidentifiable (100% divergent transitions).
#
# family = "negbin" for BOTH series (not "normal" for WW, an earlier bug in
# this port): the paper's model is Yit ~ NB(mu_it, theta_i) for both
# wastewater and outcome data, confirmed by every historical job command in
# the original repo's code/bash_commands/*.txt (all use family=negbin).
#
# min_start = 10: the paper states "Model fit only after outcome values
# exceed 10," matching min_start=10 in every historical job command. An
# earlier version of this port left this at prep_stan_data's own default of 0.
#
# upper = 'X50.' in organize_results_data/detect_outbreaks below: the paper's
# outbreak-end criterion is "eta's posterior median < 0," and every call in
# the original repo's post_processing_WW.R explicitly passes upper = 'X50.'.
# An earlier version of this port relied on organize_results_data's own
# default of upper = 'X95.' (a much stricter, rarely-triggered criterion).

source("R/00_config.R")
source("R/stan_exponential_growth.R")
library(readr)
library(dplyr)

if (requireNamespace("rstan", quietly = TRUE)) {
  library(rstan) # must be attached (not just namespace-qualified) for summary() to S4-dispatch on stanfit objects
  rstan_options(auto_write = TRUE) # cache the compiled model so we don't recompile it 20 times
} else if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("This script requires rstan or cmdstanr. See the header comment in this file for setup instructions.")
}

ww_cleaned <- read_csv(WW_CLEANED_FILE, show_col_types = FALSE)
cases_cleaned <- read_csv(CASES_CLEANED_FILE, show_col_types = FALSE)

NSAMPLE <- 1000
BURNIN <- 500
WINDOW <- 5

message(sprintf(paste(
  "NOTE: this demo uses nsample = %d (the paper's own analysis used nsample = 10000).",
  "This is a deliberate simplification to keep runtime reasonable for a public demo --",
  "posterior estimates here (R-hat, effective sample size) are less reliable than in the",
  "paper. Increase NSAMPLE/BURNIN above for more trustworthy inference.", sep = "\n"
), NSAMPLE))

fit_one_series <- function(df, out_col, family = "negbin") {
  stan_data <- prep_stan_data(df, time_unit = "week", out_col = out_col, family = family,
                               window = WINDOW, min_start = 10)
  stan_lst <- run_stan_EWS(stan_data, nsample = NSAMPLE, burnin = BURNIN)
  organized <- organize_results_data(stan_summary = stan_lst$stan_summary, stan_data = stan_data,
                                      out_col = out_col, lower = "X5.", upper = "X50.")
  outbreaks <- detect_outbreaks(organized$df, organized$eta_est, outbreak_end = "downtrend")
  outbreaks$matched_outbreaks
}

all_outbreaks <- list()

for (fip in FIPS_SUBSET) {
  county_name <- FIPS_COUNTY_NAMES[[fip]]
  cat(sprintf("Fitting Stan model for %s (%s)...\n", fip, county_name))

  df_WW <- ww_cleaned %>% filter(fips == fip) %>% select(date, Y_interp)
  df_cases <- cases_cleaned %>% filter(fips == fip) %>% select(date, new_cases_smoothed)

  ww_outbreaks <- fit_one_series(df_WW, out_col = "Y_interp")
  cases_outbreaks <- fit_one_series(df_cases, out_col = "new_cases_smoothed")

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
