# Peak-classification and evaluation metrics, ported from the original
# wastewater_EWS repo's WW_EWS_functions.R. These already operate generically
# on date/y data frames and have no data-source coupling.
library(dplyr)

### Function to get peaks for WW and cases
get_peaks_WW_cases <- function(df_WW, df_cases, null_cases_score_thresh = F, null_cases_quantile_threshold = F, align_NAs = F, cases_peak_dates = NULL, ...) {
  # Error checking: ensure data frames have required columns
  if (!"date" %in% names(df_WW) || !"y" %in% names(df_WW)) {
    stop("df_WW must have 'date' and 'y' columns")
  }
  if (!"date" %in% names(df_cases) || !"y" %in% names(df_cases)) {
    stop("df_cases must have 'date' and 'y' columns")
  }
  
  # Ungroup if grouped (to avoid warnings when selecting columns)
  df_WW <- df_WW %>% ungroup()
  df_cases <- df_cases %>% ungroup()
  
  # Filter to just date and y columns
  df_WW <- df_WW %>% select(date, y)
  df_cases <- df_cases %>% select(date, y)
  
  # Ensure dates are Date objects
  if (!inherits(df_WW$date, "Date")) df_WW$date <- as.Date(df_WW$date)
  if (!inherits(df_cases$date, "Date")) df_cases$date <- as.Date(df_cases$date)
  
  # Error checking: ensure time frames have overlap and filter to overlapping range
  min_date_WW <- min(df_WW$date, na.rm = TRUE)
  max_date_WW <- max(df_WW$date, na.rm = TRUE)
  min_date_cases <- min(df_cases$date, na.rm = TRUE)
  max_date_cases <- max(df_cases$date, na.rm = TRUE)
  
  # Find overlapping date range
  min_date <- max(min_date_WW, min_date_cases)
  max_date <- min(max_date_WW, max_date_cases)
  
  if (min_date > max_date) {
    stop(paste0("Time frames do not overlap. WW: [", min_date_WW, ", ", max_date_WW, 
                "], Cases: [", min_date_cases, ", ", max_date_cases, "]"))
  }
  
  # Filter to overlapping date range
  df_WW <- df_WW %>% filter(date >= min_date, date <= max_date)
  df_cases <- df_cases %>% filter(date >= min_date, date <= max_date)
  
  # Check if dates align after filtering
  dates_WW <- unique(df_WW$date)
  dates_cases <- unique(df_cases$date)
  
  if (length(dates_WW) != length(dates_cases) || !all(dates_WW == dates_cases)) {
    warning(paste0("After filtering to overlapping date range, data frames have different dates. WW: ",
                   length(dates_WW), " unique dates, Cases: ", length(dates_cases), " unique dates. Aligning dates by filling missing values with NA."))
    
    # Use full join to align dates - this will create NA values where dates are missing
    df_merged <- merge(df_WW, df_cases, by = "date", all = TRUE, suffixes = c("_WW", "_cases"))
    
    # if aligning the NAs between datasets.
    if(align_NAs){
      # Identify dates where either dataset has NA
      missing_in_either <- is.na(df_merged$y_WW) | is.na(df_merged$y_cases)
      
      # Set both to NA where either was missing
      df_merged$y_WW[missing_in_either] <- NA_real_
      df_merged$y_cases[missing_in_either] <- NA_real_
    }
    
    # Rename columns back to y for each dataset
    df_WW <- df_merged %>%
      select(date, y = y_WW) %>%
      arrange(date)
    df_cases <- df_merged %>%
      select(date, y = y_cases) %>%
      arrange(date)
  } else {
    # Dates already align, just ensure they're in the same order
    df_WW <- df_WW %>% arrange(date)
    df_cases <- df_cases %>% arrange(date)
  }
  # df_WW <- df_WW %>% arrange(date)
  # df_cases <- df_cases %>% arrange(date)
  
  # Extract y values and dates
  y_WW <- as.numeric(df_WW$y)
  y_cases <- as.numeric(df_cases$y)
  dates <- df_WW$date
  
  # Call find_simple_peaks for WW (all parameters passed through ...)
  peaks_WW_idx <- find_simple_peaks(y = y_WW, x = dates, ...)
  
  # Convert WW peak indices to dates
  peaks_WW_dates <- if (length(peaks_WW_idx) > 0) dates[peaks_WW_idx] else as.Date(character(0))
  
  # Cases peaks can be supplied externally (curated), otherwise compute from df_cases
  if (!is.null(cases_peak_dates)) {
    if (!inherits(cases_peak_dates, "Date")) {
      cases_peak_dates <- as.Date(cases_peak_dates)
    }
    peaks_cases_dates <- sort(unique(cases_peak_dates[!is.na(cases_peak_dates)]))
    peaks_cases_dates <- peaks_cases_dates[peaks_cases_dates >= min_date & peaks_cases_dates <= max_date]
  } else {
    # Call find_simple_peaks for cases (all parameters passed through ...)
    if(null_cases_score_thresh || null_cases_quantile_threshold){
      # Capture ... into a list, modify it, then use do.call
      dots_list <- list(...)
      if(null_cases_score_thresh){
        dots_list$score_thresh <- NULL
      }
      if(null_cases_quantile_threshold){
        dots_list$quantile_threshold <- NULL
        dots_list$quantile_lookback <- NULL
      }
      peaks_cases_idx <- do.call(find_simple_peaks, c(list(y = y_cases, x = dates), dots_list))
    } else {
      peaks_cases_idx <- find_simple_peaks(y = y_cases, x = dates, ...)
    }
    
    peaks_cases_dates <- if (length(peaks_cases_idx) > 0) dates[peaks_cases_idx] else as.Date(character(0))
  }
  
  # Return list with peaks as dates, plus date range
  return(list(
    WW = peaks_WW_dates,
    cases = peaks_cases_dates,
    min_date = min_date,
    max_date = max_date
  ))
}

### Function to classify peaks for WW and cases
classify_peaks_WW_cases <- function(df_WW, df_cases, left_w, right_w, left_TP_w = 3, right_TP_w = 3, 
                                    quantile_lookback = NULL, censor_left_w = NULL, censor_right_w = NULL, cases_peak_dates = NULL, ...) {
  
  # Determine censoring windows
  # If censor_left_w/censor_right_w are provided, use them for censoring
  # Otherwise, use the peak detection windows (left_w/right_w)
  censor_left <- if (!is.null(censor_left_w)) censor_left_w else left_w
  censor_right <- if (!is.null(censor_right_w)) censor_right_w else right_w
  
  # Adjust censor_left for quantile_lookback if provided
  # quantile_lookback requires at least that many weeks of data before a peak
  effective_censor_left <- if (!is.null(quantile_lookback)) max(quantile_lookback, censor_left) else censor_left
  
  # Convert weeks to days for date calculations
  left_TP_w_days <- left_TP_w * 7
  right_TP_w_days <- right_TP_w * 7
  left_w_days <- effective_censor_left * 7  # Use effective_censor_left for censoring
  right_w_days <- censor_right * 7  # Use censor_right for censoring
  
  # Get peaks using get_peaks_WW_cases (pass original left_w and quantile_lookback through ...)
  peaks_result <- get_peaks_WW_cases(
    df_WW, df_cases,
    left_w = left_w,
    right_w = right_w,
    quantile_lookback = quantile_lookback,
    cases_peak_dates = cases_peak_dates,
    ...
  )
  
  WW_peaks <- peaks_result$WW
  cases_peaks <- peaks_result$cases
  min_date <- peaks_result$min_date
  max_date <- peaks_result$max_date
  
  # Calculate censor values (using effective_left_w)
  LCV_WW <- min_date + left_w_days + right_TP_w_days   # Left censor value for WW
  RCV_WW <- max_date - right_w_days - left_TP_w_days   # Right censor value for WW
  
  LCV_cases <- min_date + left_w_days + left_TP_w_days # Left censor value for cases (unchanged)
  RCV_cases <- max_date - right_w_days - right_TP_w_days # Right censor value for cases (unchanged)
  
  # Initialize classification vectors
  TP <- as.Date(character(0))  # True Positives: WW peaks that match cases peaks
  FP <- as.Date(character(0))  # False Positives: WW peaks that don't match any cases peak
  FN <- as.Date(character(0))  # False Negatives: Cases peaks without matching WW peaks
  FP_excluded <- as.Date(character(0))  # False Positives excluded due to censoring
  FN_excluded <- as.Date(character(0))  # False Negatives excluded due to censoring
  TP_excluded <- as.Date(character(0))  # True Positives excluded due to censoring
  TP_excluded2 <- as.Date(character(0))  # True Positives excluded due to censoring
  
  # Track which cases peaks have been matched and which WW peak matched them
  cases_matched <- logical(length(cases_peaks))
  ww_to_cases_match <- list()  # Map WW peak to matched cases peak
  
  # For each WW peak, check if it's within the window of any cases peak
  for (i in seq_along(WW_peaks)) {
    ww_peak <- WW_peaks[i]
    matched <- FALSE
    
    # Check against each cases peak
    for (j in seq_along(cases_peaks)) {
      cases_peak <- cases_peaks[j]
      
      # Check if WW peak is within [cases_peak - left_TP_w_days, cases_peak + right_TP_w_days]
      if (ww_peak >= (cases_peak - left_TP_w_days) && ww_peak <= (cases_peak + right_TP_w_days)) {
        TP <- c(TP, ww_peak)
        cases_matched[j] <- TRUE
        ww_to_cases_match[[format(as.Date(ww_peak, origin = "1970-01-01"))]] <- cases_peak
        matched <- TRUE
        break  # Each WW peak can only match one cases peak
      }
    }
    
    # If no match found, it's a false positive
    if (!matched) {
      FP <- c(FP, ww_peak)
    }
  }
  # Cases peaks without matching WW peaks are false negatives
  # But we need to check: a cases peak that wasn't matched might still fall within
  # a WW peak's window that was already matched to another cases peak.
  # Only cases peaks that don't fall within ANY WW peak's window are true FNs.
  unmatched_cases_peaks <- cases_peaks[!cases_matched]
  
  for (i in seq_along(unmatched_cases_peaks)) {
    cases_peak <- unmatched_cases_peaks[i]
    falls_within_any_ww_window <- FALSE
    
    # Check if this cases peak falls within the window of any WW peak
    for (j in seq_along(WW_peaks)) {
      ww_peak <- WW_peaks[j]
      
      # A cases peak falls within a WW peak's window if:
      # cases_peak is within [ww_peak - right_TP_w_days, ww_peak + left_TP_w_days]
      # (This is the inverse of: ww_peak is within [cases_peak - left_TP_w_days, cases_peak + right_TP_w_days])
      if (cases_peak >= (ww_peak - right_TP_w_days) && cases_peak <= (ww_peak + left_TP_w_days)) {
        falls_within_any_ww_window <- TRUE
        break  # Found a WW peak window, no need to check others
      }
    }
    
    # Only add to FN if it doesn't fall within any WW peak's window
    if (!falls_within_any_ww_window) {
      FN <- c(FN, cases_peak)
    }else{
      TP_excluded2 <- c(TP_excluded2, cases_peak)
    }
  }
  
  # Apply censoring rules for left side
  # WW peaks before LCV
  WW_before_LCV <- WW_peaks[WW_peaks < LCV_WW]
  for (i in seq_along(WW_before_LCV)) {
    ww_peak <- WW_before_LCV[i]
    # Check if this WW peak is a TP for a cases peak on or after LCV
    # Use the same date-key format as when ww_to_cases_match is populated.
    # (ww_peak may not always be a Date; coercing makes the lookup robust.)
    matched_cases_peak <- ww_to_cases_match[[format(as.Date(ww_peak, origin = "1970-01-01"))]]
    
    if (!is.null(matched_cases_peak) && matched_cases_peak >= LCV_cases) {
      # Keep as TP (matched to cases peak on or after LCV)
      # Do nothing
    } else {
      # Either matched to cases peak before LCV, or no match (FP)
      # Remove from TP if it was there (matched to cases peak before LCV)
      if (ww_peak %in% TP) {
        TP <- TP[TP != ww_peak]
        TP_excluded <- c(TP_excluded, ww_peak)
      }
      # Remove from FP if it was there, add to FP_excluded
      if (ww_peak %in% FP) {
        FP <- FP[FP != ww_peak]
        FP_excluded <- c(FP_excluded, ww_peak)
      }
    }
  }
  
  # Cases peaks before LCV without matching WW peaks
  cases_before_LCV <- cases_peaks[cases_peaks < LCV_cases]
  for (i in seq_along(cases_before_LCV)) {
    cases_peak <- cases_before_LCV[i]
    # If this cases peak is in FN (no matching WW peak), exclude it
    if (cases_peak %in% FN) {
      FN <- FN[FN != cases_peak]
      FN_excluded <- c(FN_excluded, cases_peak)
    }
  }
  
  # Apply censoring rules for right side
  # WW peaks after RCV
  WW_after_RCV <- WW_peaks[WW_peaks > RCV_WW]
  for (i in seq_along(WW_after_RCV)) {
    ww_peak <- WW_after_RCV[i]
    # Check if this WW peak is a TP for a cases peak on or before RCV
    # Use the same date-key format as when ww_to_cases_match is populated.
    matched_cases_peak <- ww_to_cases_match[[format(as.Date(ww_peak, origin = "1970-01-01"))]]
    
    if (!is.null(matched_cases_peak) && matched_cases_peak <= RCV_cases) {
      # Keep as TP (matched to cases peak on or before RCV)
      # Do nothing
    } else {
      # Either matched to cases peak after RCV, or no match (FP)
      # Remove from TP if it was there (matched to cases peak after RCV)
      if (ww_peak %in% TP) {
        TP <- TP[TP != ww_peak]
        TP_excluded <- c(TP_excluded, ww_peak)
      }
      # Remove from FP if it was there, add to FP_excluded
      if (ww_peak %in% FP) {
        FP <- FP[FP != ww_peak]
        FP_excluded <- c(FP_excluded, ww_peak)
      }
    }
  }
  
  # Cases peaks after RCV without matching WW peaks
  cases_after_RCV <- cases_peaks[cases_peaks > RCV_cases]
  for (i in seq_along(cases_after_RCV)) {
    cases_peak <- cases_after_RCV[i]
    # If this cases peak is in FN (no matching WW peak), exclude it
    if (cases_peak %in% FN) {
      FN <- FN[FN != cases_peak]
      FN_excluded <- c(FN_excluded, cases_peak)
    }
  }
  
  # Return classification results
  return(list(
    TP = TP,
    FP = FP,
    FN = FN,
    FP_excluded = FP_excluded,
    FN_excluded = FN_excluded,
    TP_excluded = TP_excluded,
    TP_excluded2 = TP_excluded2,
    WW_peaks = WW_peaks,
    cases_peaks = cases_peaks,
    min_date = min_date,
    max_date = max_date,
    LCV_WW = LCV_WW,
    RCV_WW = RCV_WW,
    LCV_cases = LCV_cases,
    RCV_cases = RCV_cases,
    ww_to_cases_match = ww_to_cases_match
  ))
}

### Function to quantify classification components 
quantify_classification_components <- function(classification_result, df_WW, df_cases, left_w, right_w, left_TP_w = 3, right_TP_w = 3) {
  # Error checking
  if (is.null(classification_result$TP) || is.null(classification_result$FP) || is.null(classification_result$FN)) {
    stop("classification_result must contain TP, FP, and FN vectors")
  }
  if (!"date" %in% names(df_WW) || !"y" %in% names(df_WW)) {
    stop("df_WW must have 'date' and 'y' columns")
  }
  if (!"date" %in% names(df_cases) || !"y" %in% names(df_cases)) {
    stop("df_cases must have 'date' and 'y' columns")
  }
  
  # Convert weeks to days
  left_TP_w_days <- left_TP_w * 7
  right_TP_w_days <- right_TP_w * 7
  left_w_days <- left_w * 7
  right_w_days <- right_w * 7
  
  # Get date range from classification result
  min_date <- classification_result$min_date
  max_date <- classification_result$max_date
  
  # Prepare WW and cases data (ensure they're ordered by date and filtered to relevant range)
  df_WW_ordered <- df_WW %>% 
    filter(date >= min_date, date <= max_date) %>%
    arrange(date)
  
  df_cases_ordered <- df_cases %>% 
    filter(date >= min_date, date <= max_date) %>%
    arrange(date)
  
  # Get all values for quantile calculations (non-NA values)
  all_cases_values <- df_cases_ordered$y[!is.na(df_cases_ordered$y)]
  all_WW_values <- df_WW_ordered$y[!is.na(df_WW_ordered$y)]
  
  # Helper function to calculate proportion of missing WW data in a window
  calc_observed_prop <- function(start_date, end_date, df_WW) {
    # Get all dates in the window
    window_dates <- seq(start_date, end_date, by = "week")
    # Filter to dates that exist in df_WW
    window_data <- df_WW %>% filter(date %in% window_dates)
    
    # Calculate proportion missing
    if (nrow(window_data) == 0) {
      return(NA_real_)
    }
    
    n_observed <- nrow(window_data)
    n_total <- length(window_dates)
    
    if(!(n_total %in% c(7, 13))){stop('ERROR on date length')}
    return(n_observed / n_total)
  }
  
  # Helper function to get value at a specific date
  get_value_at_date <- function(df, target_date) {
    value <- df$y[df$date == target_date]
    if (length(value) > 0) {
      return(value[1])  # Return first match if multiple
    } else {
      return(NA_real_)
    }
  }
  
  # Quantify TP: proportion of observed WW data in WW window (no quantile)
  TP <- classification_result$TP
  TP_observed_prop <- numeric(length(TP))
  
  for (i in seq_along(TP)) {
    ww_peak_date <- TP[i]
    window_start <- ww_peak_date - left_w_days
    window_end <- ww_peak_date + right_w_days
    TP_observed_prop[i] <- calc_observed_prop(window_start, window_end, df_WW_ordered)
  }
  
  # Quantify FP: proportion of observed WW data in WW window + quantiles of cases data
  FP <- classification_result$FP
  FP_observed_prop <- numeric(length(FP))
  FP_quantile_same <- numeric(length(FP))  # Quantile of cases data within detection window around WW FP
  FP_quantile_next <- numeric(length(FP))   # Quantile of cases one week later
  
  for (i in seq_along(FP)) {
    ww_peak_date <- FP[i]
    window_start <- ww_peak_date - left_w_days
    window_end <- ww_peak_date + right_w_days
    FP_observed_prop[i] <- calc_observed_prop(window_start, window_end, df_WW_ordered)
    
    # Quantile of cases data at WW FP peak date relative to cases values in detection window
    # Get all cases values in the detection window
    cases_window_data <- df_cases_ordered %>% 
      filter(date >= window_start, date <= window_end)
    cases_window_values <- cases_window_data$y[!is.na(cases_window_data$y)]
    
    # Get cases value at the WW FP peak date
    cases_value_at_peak <- get_value_at_date(df_cases_ordered, ww_peak_date)
    
    if (!is.na(cases_value_at_peak) && length(cases_window_values) > 0) {
      # Calculate quantile: what proportion of values in the window are <= the value at peak date
      FP_quantile_same[i] <- mean(cases_window_values <= cases_value_at_peak, na.rm = TRUE)
    } else {
      FP_quantile_same[i] <- NA_real_
    }
    
    # Quantile of cases data one week later (relative to cases values in shifted window)
    next_week_date <- ww_peak_date + 7
    # Window shifted one week forward: [WW FP - left_W + 1 week, WW FP + right_W + 1 week]
    window_start_shifted <- ww_peak_date - left_w_days + 7
    window_end_shifted <- ww_peak_date + right_w_days + 7
    cases_window_data_shifted <- df_cases_ordered %>% 
      filter(date >= window_start_shifted, date <= window_end_shifted)
    cases_window_values_shifted <- cases_window_data_shifted$y[!is.na(cases_window_data_shifted$y)]
    
    cases_value_next <- get_value_at_date(df_cases_ordered, next_week_date)
    if (!is.na(cases_value_next) && length(cases_window_values_shifted) > 0) {
      # Calculate quantile: what proportion of values in the shifted window are <= the value one week later
      FP_quantile_next[i] <- mean(cases_window_values_shifted <= cases_value_next, na.rm = TRUE)
    } else {
      FP_quantile_next[i] <- NA_real_
    }
  }
  
  # Quantify FN: proportion of observed WW data in extended window + quantiles of WW data
  FN <- classification_result$FN
  FN_observed_prop <- numeric(length(FN))
  FN_quantile_same <- numeric(length(FN))  # Quantile of WW data within detection window around cases FN
  FN_quantile_prior <- numeric(length(FN))   # Quantile of WW one week prior
  
  for (i in seq_along(FN)) {
    cases_peak_date <- FN[i]
    # Extended window for observed data calculation
    window_start_extended <- cases_peak_date - left_w_days - left_TP_w_days
    window_end_extended <- cases_peak_date + right_w_days + right_TP_w_days
    FN_observed_prop[i] <- calc_observed_prop(window_start_extended, window_end_extended, df_WW_ordered)
    
    # Detection window for quantile calculation: cases_peak_date - left_w to cases_peak_date + right_w
    window_start_detection <- cases_peak_date - left_w_days
    window_end_detection <- cases_peak_date + right_w_days
    
    # Get all WW values in the detection window
    ww_window_data <- df_WW_ordered %>% 
      filter(date >= window_start_detection, date <= window_end_detection)
    ww_window_values <- ww_window_data$y[!is.na(ww_window_data$y)]
    
    # Get WW value at the cases FN peak date
    ww_value_at_peak <- get_value_at_date(df_WW_ordered, cases_peak_date)
    
    if (!is.na(ww_value_at_peak) && length(ww_window_values) > 0) {
      # Calculate quantile: what proportion of values in the window are <= the value at peak date
      FN_quantile_same[i] <- mean(ww_window_values <= ww_value_at_peak, na.rm = TRUE)
    } else {
      FN_quantile_same[i] <- NA_real_
    }
    
    # Quantile of WW data one week prior (relative to WW values in shifted window)
    prior_week_date <- cases_peak_date - 7
    # Window shifted one week backward: [cases FN - left_W - 1 week, cases FN + right_W - 1 week]
    window_start_shifted <- cases_peak_date - left_w_days - 7
    window_end_shifted <- cases_peak_date + right_w_days - 7
    ww_window_data_shifted <- df_WW_ordered %>% 
      filter(date >= window_start_shifted, date <= window_end_shifted)
    ww_window_values_shifted <- ww_window_data_shifted$y[!is.na(ww_window_data_shifted$y)]
    
    ww_value_prior <- get_value_at_date(df_WW_ordered, prior_week_date)
    if (!is.na(ww_value_prior) && length(ww_window_values_shifted) > 0) {
      # Calculate quantile: what proportion of values in the shifted window are <= the value one week prior
      FN_quantile_prior[i] <- mean(ww_window_values_shifted <= ww_value_prior, na.rm = TRUE)
    } else {
      FN_quantile_prior[i] <- NA_real_
    }
  }
  
  # Return quantification results
  return(list(
    TP_observed_prop = TP_observed_prop,
    FP_observed_prop = FP_observed_prop,
    FP_quantile_same = FP_quantile_same,
    FP_quantile_next = FP_quantile_next,
    FN_observed_prop = FN_observed_prop,
    FN_quantile_same = FN_quantile_same,
    FN_quantile_prior = FN_quantile_prior
  ))
}

### Compute sensitivity of a df with a "class" column (TP/FP/FN).
sensitivity <- function(df){
  return(sum(df$class == 'TP')/(sum(df$class == 'TP') + sum(df$class == 'FN')))
}

### Compute positive predictive value of a df with a "class" column.
PPV <- function(df){
  return(sum(df$class == 'TP')/(sum(df$class == 'TP') + sum(df$class == 'FP')))
}

### Compute the F1 score of a df with a "class" column.
F1 <- function(df){
  return(2*sum(df$class == 'TP')/(2*sum(df$class == 'TP') + sum(df$class == 'FP') + sum(df$class == 'FN')))
}
