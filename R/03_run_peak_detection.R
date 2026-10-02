# Peak-detection method: for each of the 10 FIPS_SUBSET counties, detects
# peaks independently in the wastewater (WW) and cases series and classifies
# WW peaks as TP/FP against cases peaks (treating cases peaks as the
# reference), using the same window/threshold convention as the original
# wastewater_EWS repo's county-level analysis (w = 3, threshold = 10,
# complete_window_only = TRUE; see code/peak_detection.R in that repo).
#
# Writes per-county classified peak dates and aggregate sensitivity/PPV/F1 to
# results/peak_detection_results.csv and results/peak_detection_summary.csv.

source("R/00_config.R")
source("R/peak_detection_algorithm.R")
source("R/peak_classification_metrics.R")
library(readr)
library(dplyr)

ww_cleaned <- read_csv(WW_CLEANED_FILE, show_col_types = FALSE)
cases_cleaned <- read_csv(CASES_CLEANED_FILE, show_col_types = FALSE)

PEAK_PARAMS <- list(w = 3, threshold = 10, complete_window_only = TRUE)

all_classified <- list()
all_summary <- list()

for (fip in FIPS_SUBSET) {
  county_name <- FIPS_COUNTY_NAMES[[fip]]

  df_WW <- ww_cleaned %>% filter(fips == fip) %>% transmute(date, y = Y_interp)
  df_cases <- cases_cleaned %>% filter(fips == fip) %>% transmute(date, y = new_cases_smoothed)

  # Pre-flight overlap check (get_peaks_WW_cases also checks this internally,
  # but we want to name the offending county rather than fail mid-loop).
  overlap_start <- max(min(df_WW$date), min(df_cases$date))
  overlap_end <- min(max(df_WW$date), max(df_cases$date))
  if (overlap_start > overlap_end) {
    warning(sprintf("Skipping %s (%s): WW and cases date ranges do not overlap.", fip, county_name))
    next
  }

  result <- do.call(classify_peaks_WW_cases, c(
    list(df_WW = df_WW, df_cases = df_cases, left_w = PEAK_PARAMS$w, right_w = PEAK_PARAMS$w),
    PEAK_PARAMS[c("threshold", "complete_window_only")]
  ))

  make_class_df <- function(dates, class_label) {
    if (length(dates) == 0) return(NULL)
    data.frame(fips = fip, county = county_name, class = class_label, date = dates)
  }
  classified <- bind_rows(
    make_class_df(result$TP, "TP"),
    make_class_df(result$FP, "FP"),
    make_class_df(result$FN, "FN")
  )
  all_classified[[fip]] <- classified

  n_tp <- length(result$TP)
  n_fp <- length(result$FP)
  n_fn <- length(result$FN)
  all_summary[[fip]] <- data.frame(
    fips = fip, county = county_name,
    n_WW_peaks = length(result$WW_peaks), n_cases_peaks = length(result$cases_peaks),
    TP = n_tp, FP = n_fp, FN = n_fn,
    sensitivity = n_tp / (n_tp + n_fn),
    PPV = n_tp / (n_tp + n_fp)
  )
}

classified_df <- bind_rows(all_classified)
summary_df <- bind_rows(all_summary)

overall <- data.frame(
  fips = "ALL", county = "All 10 counties combined",
  n_WW_peaks = sum(summary_df$n_WW_peaks), n_cases_peaks = sum(summary_df$n_cases_peaks),
  TP = sum(summary_df$TP), FP = sum(summary_df$FP), FN = sum(summary_df$FN),
  sensitivity = sensitivity(classified_df), PPV = PPV(classified_df)
)
summary_df <- bind_rows(summary_df, overall)

cat("Peak-detection results by county:\n")
print(summary_df)

write_csv(classified_df, "results/peak_detection_results.csv")
write_csv(summary_df, "results/peak_detection_summary.csv")
cat("Wrote results/peak_detection_results.csv and results/peak_detection_summary.csv\n")
