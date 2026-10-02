# Unit tests for R/utils_data_cleaning.R.
# Uses base R stopifnot() rather than testthat (kept out of this repo's
# dependencies -- see R/00_install_packages.R).

source("R/utils_data_cleaning.R")

## detect_and_average_zeros_balanced: preserves the total sum.
x <- c(1, 2, 20, 0, 0, 0, 30, 1, 1)
adjusted <- detect_and_average_zeros_balanced(x, cutoff = 5)
stopifnot(sum(adjusted) == sum(x))
stopifnot(length(adjusted) == length(x))

## ...and redistributes a post-zero-run spike across the zero run.
x2 <- c(0, 0, 0, 10, 0, 0, 0, 40)
adjusted2 <- detect_and_average_zeros_balanced(x2, cutoff = 5)
stopifnot(sum(adjusted2) == sum(x2))
stopifnot(all(adjusted2 >= 0))

## No zero runs after the cutoff is reached -> unchanged.
x3 <- c(1, 2, 10, 11, 12, 13)
stopifnot(identical(detect_and_average_zeros_balanced(x3, cutoff = 5), x3))

## adjust_negative_cases_backwards: no negative values in the output, sum preserved.
y <- c(10, 20, -5, 15, 30, -40, 50)
adjusted_y <- adjust_negative_cases_backwards(y)
stopifnot(all(adjusted_y >= 0))
stopifnot(sum(adjusted_y) == sum(y))

## A series with no negatives is returned unchanged.
y2 <- c(1, 2, 3, 4, 5)
stopifnot(identical(adjust_negative_cases_backwards(y2), y2))

## interpolate_if_few_missing: interpolates short gaps, leaves long gaps as NA.
df <- data.frame(
  date = seq(as.Date("2022-01-01"), by = "week", length.out = 10),
  Y = c(10, NA, 30, 40, NA, NA, NA, NA, 90, 100)
)
result <- interpolate_if_few_missing(df, interp_threshold = 1)
stopifnot(!is.na(result$Y_interp[2]))              # 1-week gap: interpolated
stopifnot(result$Y_interp[2] == 20)                 # linear interpolation between 10 and 30
stopifnot(all(is.na(result$Y_interp[5:8])))         # 4-week gap, threshold=1: left as NA

cat("test_cleaning_utils.R: all tests passed\n")
