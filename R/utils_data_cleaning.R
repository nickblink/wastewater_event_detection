# Source-agnostic time-series cleaning helpers, ported from the original
# wastewater_EWS repo's WW_EWS_functions.R and data_preprocessing.R. These
# operate on generic date/value columns or bare numeric vectors and have no
# dependency on any particular data source.

library(dplyr)

### Linearly interpolate across short gaps of consecutive NAs in a date-indexed
### series, leaving longer gaps untouched.
### df: data frame with `date` and `Y` columns.
### interp_threshold: maximum run length of consecutive NAs to interpolate across.
interpolate_if_few_missing <- function(df, interp_threshold) {
  df <- df %>% arrange(date)
  na_runs <- rle(is.na(df$Y))
  na_replacement <- rep(na_runs$lengths, na_runs$lengths) <= interp_threshold
  interp_vals <- approx(x = df$date[!is.na(df$Y)], y = df$Y[!is.na(df$Y)],
                        xout = df$date, rule = 1)$y
  df$Y_interp <- ifelse(is.na(df$Y) & na_replacement, interp_vals, df$Y)
  return(df)
}

### Fill in a complete weekly date grid for one location's series, interpolate
### short gaps, and drop rows that are still missing afterward.
### df: single-location data frame with a `date` column and an `outcome` column.
### cols_to_keep: metadata columns that are constant within this df and should
###   be carried through onto the filled-in dates.
interp_wrapper <- function(df, interp_threshold, outcome, cols_to_keep = character(0)) {
  df$Y <- df[[outcome]]
  full_dates <- data.frame(date = seq(min(df$date), max(df$date), by = "week"))
  df <- merge(full_dates, df, by = "date", all.x = TRUE)
  for (col in cols_to_keep) {
    col_val <- unique(df[[col]]) %>% na.omit()
    if (length(col_val) == 1) {
      df[[col]] <- rep(col_val, nrow(df))
    } else {
      stop(sprintf("Column '%s' has multiple values within one location's data.", col))
    }
  }
  if (interp_threshold > 0) {
    df <- interpolate_if_few_missing(df, interp_threshold)
  } else {
    df$Y_interp <- df$Y
  }
  df <- df[!is.na(df$Y_interp), ]
  return(df)
}

### Redistribute a post-zero-run reporting spike evenly backward across the
### zero run, preserving the exact total sum. Handles the common artifact where
### a health department batches several days/weeks of counts into one report
### after a run of zero-reporting days.
### outcome_values: numeric vector of non-negative counts, in chronological order.
detect_and_average_zeros_balanced <- function(outcome_values, cutoff = 5) {
  adjusted_values <- outcome_values
  n <- length(outcome_values)
  cutoff_reached <- FALSE
  original_sum <- sum(outcome_values)
  if (n >= 3) {
    i <- 1
    while (i <= n) {
      if (!cutoff_reached && outcome_values[i] >= cutoff) {
        cutoff_reached <- TRUE
      }
      if (cutoff_reached) {
        while (i <= n && outcome_values[i] != 0) { i <- i + 1 }
        j <- i
        while (j <= n && outcome_values[j] == 0) { j <- j + 1 }
        if (j <= n && outcome_values[j] > 0) {
          avg_value <- outcome_values[j] / (j - i + 1)
          rounded_values <- rep(floor(avg_value), j - i + 1)
          rounded_values[length(rounded_values)] <- outcome_values[j] - sum(rounded_values[-length(rounded_values)])
          adjusted_values[i:(j - 1)] <- rounded_values[-length(rounded_values)]
          adjusted_values[j] <- rounded_values[length(rounded_values)]
        }
        i <- j
      } else {
        i <- i + 1
      }
    }
  }
  adjusted_sum <- sum(adjusted_values)
  difference <- original_sum - adjusted_sum
  if (difference != 0) {
    max_index <- which.max(abs(adjusted_values - outcome_values))
    adjusted_values[max_index] <- adjusted_values[max_index] + difference
  }
  return(adjusted_values)
}

### Walk backward through a daily count series and absorb negative values
### (data corrections) into preceding day(s) instead of reporting a negative
### count, flooring at zero and carrying further back if needed.
### new_cases: numeric vector of daily (possibly negative) count diffs, in
###   chronological order.
adjust_negative_cases_backwards <- function(new_cases) {
  adjusted <- new_cases
  carry_over <- 0
  n <- length(new_cases)
  for (i in seq_along(new_cases)) {
    index <- n - i + 1
    adjusted[index] <- new_cases[index] + carry_over
    if (adjusted[index] < 0) {
      carry_over <- adjusted[index]
      adjusted[index] <- 0
    } else {
      carry_over <- 0
    }
  }
  return(adjusted)
}
