# End-to-end smoke test on a tiny synthetic 2-county fixture: builds
# wastewater and cases series in-memory (bypassing the file-path-based
# R/01.../R/02... loaders, which are wired to the real data/raw files and
# R/00_config.R globals), runs them through the same cleaning + peak-detection
# + classification functions used by the real pipeline, and checks the
# output has the expected shape and sane metric values.

source("R/utils_data_cleaning.R")
source("R/peak_detection_algorithm.R")
source("R/peak_classification_metrics.R")
library(dplyr)

make_fixture <- function(seed, n_weeks = 40) {
  set.seed(seed)
  dates <- seq(as.Date("2021-01-03"), by = "week", length.out = n_weeks)
  base <- 50 + 30 * sin(seq(0, 4 * pi, length.out = n_weeks))
  base[base < 0] <- 1
  ww <- data.frame(date = dates, Y = base + rnorm(n_weeks, sd = 3))
  cases <- data.frame(date = dates, Y = base * 10 + rnorm(n_weeks, sd = 10))
  list(ww = ww, cases = cases)
}

fixtures <- list(county_A = make_fixture(1), county_B = make_fixture(2))

all_classified <- list()
for (county in names(fixtures)) {
  ww_df <- interp_wrapper(fixtures[[county]]$ww, interp_threshold = 3, outcome = "Y")
  cases_df <- interp_wrapper(fixtures[[county]]$cases, interp_threshold = 3, outcome = "Y")

  df_WW <- ww_df %>% transmute(date, y = Y_interp)
  df_cases <- cases_df %>% transmute(date, y = Y_interp)

  result <- classify_peaks_WW_cases(df_WW, df_cases, left_w = 3, right_w = 3,
                                     threshold = 10, complete_window_only = TRUE)

  make_class_df <- function(dates, class_label) {
    if (length(dates) == 0) return(NULL)
    data.frame(county = county, class = class_label, date = dates)
  }
  all_classified[[county]] <- bind_rows(
    make_class_df(result$TP, "TP"), make_class_df(result$FP, "FP"), make_class_df(result$FN, "FN")
  )
}

classified_df <- bind_rows(all_classified)

stopifnot(nrow(classified_df) > 0)
stopifnot(all(c("county", "class", "date") %in% names(classified_df)))
stopifnot(all(classified_df$class %in% c("TP", "FP", "FN")))

sens <- sensitivity(classified_df)
ppv <- PPV(classified_df)
stopifnot(sens >= 0, sens <= 1)
stopifnot(ppv >= 0, ppv <= 1)

cat(sprintf("test_pipeline_smoke.R: ran on %d synthetic counties, sensitivity=%.2f, PPV=%.2f\n",
            length(fixtures), sens, ppv))
cat("test_pipeline_smoke.R: all tests passed\n")
