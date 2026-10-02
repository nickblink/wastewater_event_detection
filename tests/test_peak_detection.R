# Unit tests for R/peak_detection_algorithm.R.

source("R/peak_detection_algorithm.R")

## A single clean spike produces exactly one peak, at the right index.
y <- c(1, 2, 3, 20, 3, 2, 1, 1, 2, 3, 1, 1, 1)
peaks <- find_simple_peaks(y, w = 3, threshold = 10, complete_window_only = FALSE)
stopifnot(identical(as.integer(peaks), 4L))

## A flat series has no peaks.
y_flat <- rep(5, 20)
peaks_flat <- find_simple_peaks(y_flat, w = 3, threshold = 1, complete_window_only = FALSE)
stopifnot(length(peaks_flat) == 0)

## Values below threshold are not flagged, even if locally maximal.
y_below <- c(1, 2, 3, 5, 3, 2, 1, 1, 2, 3, 1, 1, 1)
peaks_below <- find_simple_peaks(y_below, w = 3, threshold = 10, complete_window_only = FALSE)
stopifnot(length(peaks_below) == 0)

## A window containing a zero is never flagged as a peak (no_zero rule).
y_zero <- c(1, 2, 0, 20, 3, 2, 1, 1, 2, 3, 1, 1, 1)
peaks_zero <- find_simple_peaks(y_zero, w = 3, threshold = 10, complete_window_only = FALSE)
stopifnot(length(peaks_zero) == 0)

## w_error_check: must specify exactly one of w or (left_w & right_w).
err1 <- tryCatch({ w_error_check(w = NULL, left_w = NULL, right_w = NULL); FALSE }, error = function(e) TRUE)
err2 <- tryCatch({ w_error_check(w = 3, left_w = 2, right_w = 2); FALSE }, error = function(e) TRUE)
stopifnot(err1, err2)
resolved <- w_error_check(w = 3, left_w = NULL, right_w = NULL)
stopifnot(resolved$left_w == 3, resolved$right_w == 3)

cat("test_peak_detection.R: all tests passed\n")
