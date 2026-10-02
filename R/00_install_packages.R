# Installs the CRAN packages used by this pipeline.
# Deliberately avoids packages that require compilation beyond what's in base R
# where a pure-base alternative exists (e.g. stats::approx instead of zoo::na.approx),
# to keep the dependency footprint minimal.
#
# The Stan exponential-growth step (R/stan_exponential_growth.R) additionally
# requires either `rstan` or `cmdstanr`, which are NOT installed here because they
# require a working C++ toolchain and can take a long time to build. See that
# file's header comment for setup instructions.

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "lubridate",
  "ggplot2"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required packages already installed.")
}
