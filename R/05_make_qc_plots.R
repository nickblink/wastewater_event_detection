# Per-county QC/diagnostic figures: raw vs. interpolated wastewater series,
# cases series, and detected peaks overlaid and colored by TP/FP/FN
# classification. If results/stan_results.csv exists (requires having run
# R/04_run_stan_exponential_growth.R, which needs rstan/cmdstanr), outbreak
# windows from the Stan method are shaded on top.
#
# This is a simplified, demo-focused rewrite -- not a port of any single
# function from the original repo's large paper-figure-generation code.

source("R/00_config.R")
library(readr)
library(dplyr)
library(ggplot2)

ww_cleaned <- read_csv(WW_CLEANED_FILE, show_col_types = FALSE)
cases_cleaned <- read_csv(CASES_CLEANED_FILE, show_col_types = FALSE)
classified <- read_csv("results/peak_detection_results.csv", show_col_types = FALSE)

stan_results_file <- "results/stan_results.csv"
stan_results <- if (file.exists(stan_results_file)) {
  read_csv(stan_results_file, show_col_types = FALSE)
} else {
  NULL
}

CLASS_COLORS <- c(TP = "#009E73", FP = "#E69F00", FN = "#0072B2")

for (fip in FIPS_SUBSET) {
  county_name <- FIPS_COUNTY_NAMES[[fip]]

  ww_county <- ww_cleaned %>% filter(fips == fip)
  cases_county <- cases_cleaned %>% filter(fips == fip)
  class_county <- classified %>% filter(fips == fip)

  ww_peaks <- class_county %>%
    filter(class %in% c("TP", "FP")) %>%
    left_join(ww_county %>% select(date, value = Y_interp), by = "date")
  cases_peaks <- class_county %>%
    filter(class == "FN") %>%
    left_join(cases_county %>% select(date, value = new_cases_smoothed), by = "date")

  p_ww <- ggplot(ww_county, aes(x = date)) +
    geom_line(aes(y = Y), color = "gray60", na.rm = TRUE) +
    geom_line(aes(y = Y_interp), color = "steelblue") +
    { if (nrow(ww_peaks) > 0) geom_point(data = ww_peaks, aes(y = value, color = class), size = 2.5) } +
    scale_color_manual(values = CLASS_COLORS, name = "WW peak") +
    labs(title = sprintf("%s -- wastewater concentration", county_name),
         subtitle = "gray = raw, blue = interpolated", x = NULL, y = "copies/L wastewater") +
    theme_minimal()

  p_cases <- ggplot(cases_county, aes(x = date, y = new_cases_smoothed)) +
    geom_line(color = "firebrick") +
    { if (nrow(cases_peaks) > 0) geom_point(data = cases_peaks, aes(y = value, color = class), size = 2.5) } +
    scale_color_manual(values = CLASS_COLORS, name = "Cases peak (unmatched = FN)") +
    labs(title = sprintf("%s -- weekly new cases", county_name), x = NULL, y = "new cases") +
    theme_minimal()

  if (!is.null(stan_results)) {
    stan_county <- stan_results %>% filter(fips == fip)
    ww_windows <- stan_county %>% filter(series == "WW")
    cases_windows <- stan_county %>% filter(series == "cases")
    if (nrow(ww_windows) > 0) {
      p_ww <- p_ww + geom_rect(data = ww_windows, inherit.aes = FALSE,
                                aes(xmin = outbreak_start, xmax = outbreak_end, ymin = -Inf, ymax = Inf),
                                alpha = 0.15, fill = "purple")
    }
    if (nrow(cases_windows) > 0) {
      p_cases <- p_cases + geom_rect(data = cases_windows, inherit.aes = FALSE,
                                      aes(xmin = outbreak_start, xmax = outbreak_end, ymin = -Inf, ymax = Inf),
                                      alpha = 0.15, fill = "purple")
    }
  }

  ggsave(sprintf("figures/%s_wastewater.png", fip), p_ww, width = 9, height = 4, dpi = 150)
  ggsave(sprintf("figures/%s_cases.png", fip), p_cases, width = 9, height = 4, dpi = 150)
}

cat(sprintf("Wrote QC figures for %d counties to figures/\n", length(FIPS_SUBSET)))
