# Main entry point for this repo's demo.
#
# Runs the full pipeline end-to-end for the 10 demo counties (R/00_config.R):
#   1) load + clean wastewater data (CDC NWSS)
#   2) load + clean cases data (JHU)
#   3) peak-detection method, on both series
#   4) Bayesian exponential-growth (Stan) method, on both series
#      -- requires rstan or cmdstanr; see R/04_run_stan_exponential_growth.R
#      -- if neither is installed, this step is skipped with a message,
#         and the rest of the pipeline still runs.
#   5) QC figures for every county
#
# Usage: Rscript R/run_full_pipeline_10_counties.R
# (run R/00_install_packages.R first if you haven't already)

cat("=== [1/5] Loading wastewater data ===\n")
source("R/01_load_ww_data.R")

cat("\n=== [2/5] Loading cases data ===\n")
source("R/02_load_cases_data.R")

cat("\n=== [3/5] Running peak-detection method (WW + cases) ===\n")
source("R/03_run_peak_detection.R")

cat("\n=== [4/5] Running Stan exponential-growth method (WW + cases) ===\n")
if (requireNamespace("rstan", quietly = TRUE) || requireNamespace("cmdstanr", quietly = TRUE)) {
  source("R/04_run_stan_exponential_growth.R")
} else {
  message("Skipping Stan step: neither rstan nor cmdstanr is installed. ",
          "See R/04_run_stan_exponential_growth.R for setup instructions. ",
          "Continuing with the rest of the pipeline.")
}

cat("\n=== [5/5] Making QC figures ===\n")
source("R/05_make_qc_plots.R")

cat("\nDone. See data/processed/, results/, and figures/.\n")
