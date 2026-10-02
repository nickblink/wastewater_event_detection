# Peak-detection algorithm, ported verbatim from the original wastewater_EWS
# repo's WW_EWS_functions.R (find_simple_peaks + its w_error_check helper).
# Pure numeric-vector algorithm with no data-source coupling.
#
# rolling_score() from the original file is intentionally not ported: it is a
# separate exported utility never called internally by find_simple_peaks, and
# this repo's demo only exercises the w/threshold/complete_window_only path.

#### Peak detection functions ####

### Helper function for window parameter error checking 
w_error_check <- function(w, left_w, right_w) {
  # Error checking: must specify either w OR (left_w AND right_w), but not both
  has_w <- !is.null(w)
  has_left_right <- !is.null(left_w) && !is.null(right_w)
  
  if (has_w && has_left_right) {
    stop("Cannot specify both 'w' and 'left_w'/'right_w'. Use either 'w' (for symmetric windows) or 'left_w' and 'right_w' (for asymmetric windows).")
  }
  if (!has_w && !has_left_right) {
    stop("Must specify either 'w' (for symmetric windows) or both 'left_w' and 'right_w' (for asymmetric windows).")
  }
  if (has_w) {
    left_w <- w
    right_w <- w
  }
  
  return(list(left_w = left_w, right_w = right_w))
}

### Peak finder 
find_simple_peaks <- function(y, x = NULL, left_w = NULL, right_w = NULL, w = NULL, threshold = 10, complete_window_only = FALSE, centered_sd_diff_threshold = NULL, score_thresh = NULL, lambda = 1.0, weight_mean = FALSE, include_current = TRUE, quantile_threshold = NULL, quantile_lookback = NULL, piecewise_slope_test = NULL, piecewise_slope_alpha = NULL, ... ) {
  # Error checking and window parameter resolution
  window_params <- w_error_check(w, left_w, right_w)
  left_w <- window_params$left_w
  right_w <- window_params$right_w
  
  stopifnot(is.numeric(y), left_w >= 1, right_w >= 1, length(y) > left_w + right_w)
  
  # Piecewise slope test parameter validation
  if (!is.null(piecewise_slope_test)) {
    valid_options <- c('linear', 'log', 'both')
    if (!(piecewise_slope_test %in% valid_options)) {
      stop("piecewise_slope_test must be one of: NULL, 'linear', 'log', 'both'")
    }
  }
  
  # Check that both parameters are either NULL or not NULL together
  if (is.null(piecewise_slope_test) != is.null(piecewise_slope_alpha)) {
    stop("piecewise_slope_test and piecewise_slope_alpha must both be NULL or both be specified")
  }
  
  if (!is.null(piecewise_slope_alpha) && !is.numeric(piecewise_slope_alpha)) {
    stop("piecewise_slope_alpha must be numeric")
  }
  if (!is.null(x)) {
    stopifnot(length(x) == length(y))
    # Convert x to Date if not already
    if (!inherits(x, "Date")) x <- as.Date(x)
  }
  n <- length(y)
  peaks <- logical(n)
  
  for (i in seq_len(n)) {
    # Only consider points where full window is available (at least left_w points before and right_w points after)
    # For left_w=8, right_w=3, this means indices from 9 to (n-3), ensuring full window on both sides
    full_window_available <- (i > left_w) && (i <= n - right_w)
    if (!full_window_available) {
      peaks[i] <- FALSE
      next
    }
    
    lo <- max(1, i - left_w)
    hi <- min(n, i + right_w)
    window <- y[lo:hi]
    
    left  <- if (i > lo) y[lo:(i - 1)] else numeric(0)
    right <- if (i < hi) y[(i + 1):hi] else numeric(0)
    
    left_strict  <- length(left)  == 0 || (all(y[i] >  left, na.rm = T))
    right_nonstr <- length(right) == 0 || (all(y[i] >= right, na.rm = T))
    
    no_zero      <- all(window != 0, na.rm = T)
    above_thresh <- !is.na(y[i]) && y[i] > threshold
    
    # Check for complete window if requested
    complete_window <- TRUE
    if (complete_window_only && !is.null(x)) {
      window_x <- x[lo:hi]
      if (length(window_x) > 1) {
        # Check that all consecutive points are no more than 1 week (7 days) apart
        diffs <- diff(window_x)
        complete_window <- all(!is.na(diffs)) && all(diffs <= 7) && all(diffs > 0)
      }
    }
    
    # Check centered SD difference threshold if requested
    above_sd_thresh <- TRUE
    if (!is.null(centered_sd_diff_threshold)) {
      # Calculate window excluding current point (include_current = FALSE)
      window_excl_current <- window
      current_idx_in_window <- i - lo + 1
      if (current_idx_in_window >= 1 && current_idx_in_window <= length(window_excl_current)) {
        window_excl_current <- window_excl_current[-current_idx_in_window]
      }
      # Remove NA values
      window_excl_current <- window_excl_current[!is.na(window_excl_current)]
      
      if (length(window_excl_current) > 1 && !is.na(y[i])) {
        sd_window <- sd(window_excl_current)
        if (!is.na(sd_window) && sd_window > 0) {
          mean_window <- mean(window_excl_current)
          centered_sd_diff <- (y[i] - mean_window) / sd_window
          above_sd_thresh <- !is.na(centered_sd_diff) && centered_sd_diff >= centered_sd_diff_threshold
          
        } else {
          above_sd_thresh <- FALSE
        }
      } else {
        above_sd_thresh <- FALSE
      }
    }
    
    # Check score threshold if requested
    above_score_thresh <- TRUE
    if (!is.null(score_thresh)) {
      # Compute score for current point (inline calculation)
      score_window <- window
      score_window_indices <- lo:hi
      
      # Optionally exclude the current point from the window
      if (!include_current) {
        current_idx_in_window <- i - lo + 1
        if (current_idx_in_window >= 1 && current_idx_in_window <= length(score_window)) {
          score_window <- score_window[-current_idx_in_window]
          score_window_indices <- score_window_indices[-current_idx_in_window]
        }
      }
      
      # Remove NA values
      na_mask <- !is.na(score_window)
      score_window <- score_window[na_mask]
      score_window_indices <- score_window_indices[na_mask]
      
      if (length(score_window) > 1 && !is.na(y[i])) {
        sd_score_window <- sd(score_window)
        if (!is.na(sd_score_window) && sd_score_window > 0) {
          # Calculate mean(window) - either regular or weighted (triangle)
          if (weight_mean) {
            distances <- abs(score_window_indices - i)
            weights <- distances
            if (sum(weights) == 0) {
              mean_score_window <- mean(score_window)
            } else {
              mean_score_window <- weighted.mean(score_window, weights)
            }
          } else {
            mean_score_window <- mean(score_window)
          }
          
          mean_diff_div_sd <- (y[i] - mean_score_window) / sd_score_window
          
          # Calculate noise_left: SD of first derivatives / range for left window [i-left_w, i] (includes current point)
          lo_left <- max(1, i - left_w)
          hi_left <- i
          left_window_score <- y[lo_left:hi_left]
          left_window_score <- left_window_score[!is.na(left_window_score)]
          
          noise_left <- 0
          if (length(left_window_score) >= 2) {
            left_first_derivs <- diff(left_window_score)
            if (length(left_first_derivs) >= 2) {
              left_range <- max(left_window_score) - min(left_window_score)
              if (!is.na(left_range) && left_range > 0) {
                left_sd_derivs <- sd(left_first_derivs)
                if (!is.na(left_sd_derivs) && left_sd_derivs > 0) {
                  noise_left <- left_sd_derivs / left_range
                }
              }
            }
          }
          
          # Calculate noise_right: SD of first derivatives / range for right window [i, i+right_w] (includes current point)
          lo_right <- i
          hi_right <- min(n, i + right_w)
          right_window_score <- y[lo_right:hi_right]
          right_window_score <- right_window_score[!is.na(right_window_score)]
          
          noise_right <- 0
          if (length(right_window_score) >= 2) {
            right_first_derivs <- diff(right_window_score)
            if (length(right_first_derivs) >= 2) {
              right_range <- max(right_window_score) - min(right_window_score)
              if (!is.na(right_range) && right_range > 0) {
                right_sd_derivs <- sd(right_first_derivs)
                if (!is.na(right_sd_derivs) && right_sd_derivs > 0) {
                  noise_right <- right_sd_derivs / right_range
                }
              }
            }
          }
          
          left_weight <- left_w / (left_w + right_w)
          right_weight <- right_w / (left_w + right_w)
          score_val <- lambda * mean_diff_div_sd - (left_weight * noise_left + right_weight * noise_right)
          above_score_thresh <- !is.na(score_val) && score_val >= score_thresh
        } else {
          above_score_thresh <- FALSE
        }
      } else {
        above_score_thresh <- FALSE
      }
    }
    
    # Check quantile threshold if requested
    above_quantile_thresh <- TRUE
    if (!is.null(quantile_threshold) && !is.null(quantile_lookback)) {
      if (is.null(x)) {
        stop("quantile_lookback requires date vector 'x' to be provided for week-based lookback")
      }
      
      # Calculate cutoff date: quantile_lookback weeks before current date
      quantile_lookback_days <- quantile_lookback * 7
      current_date <- x[i]
      cutoff_date <- current_date - quantile_lookback_days
      
      # Find all previous points within the lookback window (excluding current point)
      # Points where date >= cutoff_date and date < current_date
      previous_indices <- which(x < current_date & x >= cutoff_date)
      
      if (length(previous_indices) > 0) {
        # Get values for those previous points
        previous_values <- y[previous_indices]
        
        # Remove NA values
        previous_values <- previous_values[!is.na(previous_values)]
        
        # Check if we have at least some values and current value is not NA
        if (length(previous_values) > 0 && !is.na(y[i])) {
          quantile_val <- quantile(previous_values, probs = quantile_threshold, na.rm = TRUE)
          above_quantile_thresh <- !is.na(quantile_val) && y[i] >= quantile_val
        } else {
          # Not enough non-NA values in the lookback window
          above_quantile_thresh <- FALSE
        }
      } else {
        # No previous points within the lookback window
        above_quantile_thresh <- FALSE
      }
    }
    
    # Check piecewise slope test if requested
    passes_slope_test <- TRUE
    if (!is.null(piecewise_slope_test) && !is.null(piecewise_slope_alpha)) {
      # Build data frame for the window around the peak
      slope_lo <- max(1, i - left_w)
      slope_hi <- min(n, i + right_w)
      slope_indices <- slope_lo:slope_hi
      slope_y <- y[slope_indices]
      
      # Check we have enough non-NA data
      if (sum(!is.na(slope_y)) >= 4) {
        slope_df <- data.frame(
          y_val = slope_y,
          idx_centered = slope_indices - i,
          right = as.numeric(slope_indices > i)
        )
        slope_df <- slope_df[!is.na(slope_df$y_val), ]
        
        p_values <- c()
        
        # Linear model
        if (piecewise_slope_test %in% c('linear', 'both')) {
          tryCatch({
            fit_linear <- lm(y_val ~ idx_centered + idx_centered:right, data = slope_df)
            coef_summary <- summary(fit_linear)$coefficients
            # The interaction term (slope change) is the last row
            if (nrow(coef_summary) >= 3) {
              p_val_linear <- coef_summary[3, 4]
              p_values <- c(p_values, p_val_linear)
            }
          }, error = function(e) NULL)
        }
        
        # Log model
        if (piecewise_slope_test %in% c('log', 'both')) {
          # Only fit log model if all y values are positive
          if (all(slope_df$y_val > 0)) {
            tryCatch({
              slope_df$log_y <- log(slope_df$y_val)
              fit_log <- lm(log_y ~ idx_centered + idx_centered:right, data = slope_df)
              coef_summary <- summary(fit_log)$coefficients
              if (nrow(coef_summary) >= 3) {
                p_val_log <- coef_summary[3, 4]
                p_values <- c(p_values, p_val_log)
              }
            }, error = function(e) NULL)
          }
        }
        
        # Check if min p-value is below alpha
        if (length(p_values) > 0) {
          passes_slope_test <- min(p_values) < piecewise_slope_alpha
        } else {
          passes_slope_test <- FALSE
        }
      } else {
        passes_slope_test <- FALSE
      }
    }
    
    peaks[i] <- left_strict && right_nonstr && no_zero && above_thresh && complete_window && above_sd_thresh && above_score_thresh && above_quantile_thresh && passes_slope_test
  }
  
  which(peaks)
}
