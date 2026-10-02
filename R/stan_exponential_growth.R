# Bayesian exponential-growth outbreak detection, ported from the original
# wastewater_EWS repo's WW_EWS_functions.R. Requires rstan or cmdstanr (NOT
# installed by R/00_install_packages.R -- see note below) because fitting
# modular_outbreak_detection.stan requires a working C++ toolchain.
#
# To use this file:
#   install.packages('rstan')   # or: install cmdstanr per https://mc-stan.org/cmdstanr/
#
# Huisman/Rt-comparator and diagnostics-plotting code from the original file
# are intentionally not ported here (out of scope for this minimal demo).
library(dplyr)

prep_stan_data <- function(df, time_unit, out_col, family, remove_NAs = T, min_start = 0, rescale_data = T, constant_family_param = 1, window = 5, theta_prior_shape = .001, theta_prior_rate = .001, sigma_prior_sd = 100, beta_prior_shape = .001, beta_prior_rate = .001, sigma_eta = 1, sigma_beta = 1, eta_RW = 0, beta_RW = 0, debug = 0, transformation = NULL, ...){
  
  # remove NAs
  if(remove_NAs){
    print('removing NAs in df')
    df <- df[!is.na(df[,out_col]),]
  }
  
  if(min_start > 0){
    ind_start <- min(which(df[, out_col] > min_start))
    print(sprintf('filtering data to start at index %s (%s); value = %s', ind_start, df$date[ind_start], df[ind_start,out_col]))
    df <- df[ind_start:nrow(df),]
  }
  
  if(!is.null(transformation)){
    outcome <- transform_col(df[,out_col,drop=T], transformation)
  }else{
    outcome <- df[,out_col,drop = T]
  }
  
  # rescaling data to either have a minimum value of 1 or maximum value of 1000.
  if(rescale_data){
    # scale_factor <- min(min(outcome[outcome > 0]),
    #                     max(outcome)/1000)
    scale_factor <- max(outcome, na.rm = T)/1000
    print(sprintf('scaling data down by %s to set max value to 1000', round(scale_factor, 2)))
    Y <- outcome/scale_factor
  }else{
    scale_factor = NA
    Y <- outcome
  }
  
  # rounding outcomes to fit a discrete model.
  Y_int <- round(Y)
  warning('rounding Y for integer outcomes')
  
  # getting the time points as integers.
  T <- difftime(df$date, min(df$date), units = time_unit) %>%
    as.numeric()
  T <- T + 1
  df$T <- T
  
  # getting the family index
  # if(family == 'normal0'){
  #   family_index = 0
  if(family == 'normal'){
    family_index = 1
  }else if(family == 'poisson'){
    family_index = 2
    if(constant_family_param == 0){
      stop('Constant family param should be 1 for Poisson.')
    }
  }else if(family == 'negbin'){
    family_index = 3
  }else{
    stop('Please input a proper family')
  }
  
  # making stan data object.
  stan_data <- list(N = nrow(df),
                    N_estimated = sum(T > window),
                    window = window,
                    T = T,
                    Y = Y,
                    Y_int = Y_int,
                    family = family_index,
                    constant_family_param = constant_family_param,
                    theta_prior_shape = theta_prior_shape,
                    theta_prior_rate = theta_prior_rate,
                    #sigma2_prior_shape = sigma2_prior_shape,
                    #sigma2_prior_rate = sigma2_prior_rate,
                    sigma_prior_sd = sigma_prior_sd,
                    beta_prior_shape = beta_prior_shape,
                    beta_prior_rate = beta_prior_rate,
                    sigma_eta = sigma_eta,
                    sigma_beta = sigma_beta,
                    eta_RW = eta_RW,
                    beta_RW = beta_RW,
                    debug = debug)
  
  # returning stan data object and original data frame.
  return(list(input = stan_data, original = df, scale_factor = scale_factor, out_col = out_col))
}

run_stan_EWS <- function(stan_data, nsample = 1000, burnin = 500, use_cmdstanr = F, ...){
  print('in run_stan_EWS')
  stan_path <- 'stan/modular_outbreak_detection.stan'
  
  # set up initialization.
  init_list <- function(){list(eta = rep(0, stan_data$input$N_estimated),
                               beta = rep(mean(stan_data$input$Y), stan_data$input$N_estimated))}

  
  # run stan model!
  if(use_cmdstanr){
    library(cmdstanr)
    stan_m <- cmdstan_model(stan_path)
    print('model prepped')
    # Run sampling
    stan_fit_initial <- stan_m$sample(
      data = stan_data$input,
      init = init_list,
      iter_warmup = burnin,
      iter_sampling = nsample - burnin,  # total iterations = warmup + sampling
      chains = 1,
      seed = 1,
      thin = 10,
      save_warmup = TRUE,
      show_messages = TRUE,
      refresh = 100  # adjusts how often progress is printed; similar to verbosity
    )
  }else{
    stan_m <- rstan::stan_model(stan_path)
    print('model prepped')
    stan_fit_initial <- rstan::sampling(object = stan_m,
                                        data = stan_data$input, 
                                        init = init_list,
                                        iter = nsample, 
                                        warmup = burnin,
                                        chains = 1,
                                        seed = 1,
                                        show_messages = TRUE,
                                        verbose = FALSE, 
                                        thin = 10,
                                        save_warmup = TRUE)
  }
  
  #browser()
  
  # Extract divergence info
  sampler_params_initial <- rstan::get_sampler_params(stan_fit_initial, inc_warmup = FALSE)[[1]]
  n_iter <- nrow(sampler_params_initial)
  n_divergent_initial <- sum(sampler_params_initial[, "divergent__"])
  divergence_rate_initial <- n_divergent_initial / n_iter
  
  cat(sprintf("Initial divergence rate: %.1f%% (%d of %d iterations)\n", 
              100 * divergence_rate_initial, n_divergent_initial, n_iter))
  
  # Check if we need to rerun
  if (divergence_rate_initial > 0.01) {
    stan_model_rerun <- T
    cat("⚠️  Divergence rate above threshold (1%). Re-running with adapt_delta = 0.99 and max_treedepth = 12...\n")
    
    stan_fit_rerun <- rstan::sampling(object = stan_m,
                                      data = stan_data$input, 
                                      init = init_list,
                                      iter = nsample, 
                                      warmup = burnin,
                                      chains = 1,
                                      seed = 1,
                                      show_messages = TRUE,
                                      verbose = FALSE, 
                                      thin = 10,
                                      save_warmup = TRUE,
                                      control = list(adapt_delta = 0.99, max_treedepth = 12))
    
    sampler_params_rerun <- rstan::get_sampler_params(stan_fit_rerun, inc_warmup = FALSE)[[1]]
    n_divergent_rerun <- sum(sampler_params_rerun[, "divergent__"])
    divergence_rate_rerun <- n_divergent_rerun / n_iter
    
    cat(sprintf("✅  Rerun divergence rate: %.1f%% (%d of %d iterations)\n",
                100 * divergence_rate_rerun, n_divergent_rerun, n_iter))

    # Replace fit object with rerun version
    stan_fit <- stan_fit_rerun
  } else {
    stan_model_rerun <- F
    cat("✅  Divergence rate acceptable. No rerun needed.\n")
    
    # If not rerun, use the initial fit
    stan_fit <- stan_fit_initial
    divergence_rate_rerun <- NA
  }
  
  # get summary statistics and save.
  stan_summary <- summary(stan_fit, probs = c(0.025,0.05,0.5,0.95,0.975))

  # Compile the results - this crashes Rstudio on my desktop.
  res_lst <- list(stan_fit = stan_fit, 
                  stan_summary = stan_summary, 
                  stan_extract = rstan::extract(stan_fit),
                  stan_model_rerun = stan_model_rerun,
                  divergence_rate_initial = divergence_rate_initial,
                  divergence_rate_rerun = divergence_rate_rerun)
  
  return(res_lst)
}

get_time_interval <- function(dates, verbose = T){
  # get the intervals between successive dates.
  intervals <- dates %>% 
    sort() %>%
    diff()
  
  # get unique interval lengths.
  uni_int <- unique(intervals)
  
  if(length(uni_int) > 1){
    if(verbose){
      warning(sprintf('multiple interval lengths found. Taking the minimum (%s)', min(uni_int)))  
    }
    return(min(uni_int))
  }else{
    if(verbose){
      print(sprintf('Using time interval of length %s', uni_int))
    }
    return(uni_int)
  }
}

transform_predictions <- function(stan_summary, stan_data, out_col, time_interval = NULL, scale_down = F, verbose = T){
  # Pull out the stan data and list of predictions.
  stan_res <- stan_summary$summary %>% data.frame()
  stan_pred <- stan_res[grep('predictions', rownames(stan_res)),] 
  index_list <- lapply(rownames(stan_pred), function(text){
    numbers <- stringr::str_extract_all(text, "\\d+")
    num_list <- as.numeric(numbers[[1]])
    num_list
  })
  
  # get the time interval
  if(is.null(time_interval)){
    time_interval <- get_time_interval(stan_data$original$date)
  }
  
  # pull out the fitting date index and the adjust by the starting fit time point.
  fit_index <- sapply(index_list, function(x) x[1]) + stan_data$input$N - stan_data$input$N_estimated
  fit_id <- match(stan_data$input$T[fit_index], stan_data$original$T)
  stan_pred$fit_date <- stan_data$original$date[fit_id]
  
  # get the predict date index.
  fit_diff <- sapply(index_list, function(x) x[2])
  fit_diff <- fit_diff - max(fit_diff)
  stan_pred$pred_date <- stan_pred$fit_date + lubridate::days(fit_diff)*time_interval
  
  pred_id <- match(stan_pred$pred_date, stan_data$original$date)
  
  if(scale_down){
    stan_pred$Y <- stan_data$original[pred_id, out_col]/stan_data$scale_factor
    if(verbose){
      print(sprintf('scaling down observed Y values in model fitting by a factor of %s', stan_data$scale_factor))
    }
  }else{
    stan_pred$Y <- stan_data$original[pred_id, out_col, drop = T]
  }
  
  return(stan_pred)
}

transform_original_data <- function(df, scale_factor, out_col, time_interval = NULL, verbose = T, ...){

  # get the time interval
  if(is.null(time_interval)){
    time_interval <- get_time_interval(df$date, verbose)
  }

  # make data frame with all dates.
  all_dates <- data.frame(date = seq(min(df$date), max(df$date), by = time_interval))
  
  # merge in, left joining to add NAs.
  df_extended <- all_dates %>%
    left_join(df, by = "date")

  # scale the Y value by the scaling factor.
  df_extended$Y <- df_extended[,out_col]/scale_factor
  
  if(verbose){
    print(sprintf('scaling down observed Y values in model fitting by a factor of %s', scale_factor))
  }
  
  return(df_extended %>% select(any_of(c('date','Y','data_pts'))))
}

filter_outbreaks_and_downtrends <- function(outbreak_dates, downtrend_dates) {
  
  if(length(outbreak_dates) == 0){
    # Return filtered outbreaks and downtrend dates as a list
    return(list(filtered_outbreaks = c(), filtered_downtrends = c()))
  }
  
  # Sort outbreaks and downtrend_dates
  outbreak_dates <- sort(outbreak_dates)
  downtrend_dates <- sort(downtrend_dates)
  
  # Initialize filtered outbreaks and downtrend dates
  current_outbreak <- outbreak_dates[1]
  filtered_outbreaks <- current_outbreak # Always keep the first outbreak
  filtered_downtrends <- as.Date(character())
  
  if(length(downtrend_dates) > 0){
    while(T){
      # find the downtrend of the outbreak.
      if(any(downtrend_dates > current_outbreak)){
        current_downtrend <- min(downtrend_dates[downtrend_dates > current_outbreak])
        filtered_downtrends <- c(filtered_downtrends, current_downtrend)
      }else{
        break
      }
      
      # find the next outbreak.
      if(any(outbreak_dates > current_downtrend)){
        current_outbreak <- min(outbreak_dates[outbreak_dates > current_downtrend])
        filtered_outbreaks <- c(filtered_outbreaks, current_outbreak)
      }else{
        break
      }
    }
  }
  
  # error checking.
  if(length(filtered_downtrends) > length(filtered_outbreaks)){browser()}
  if(any(filtered_downtrends < filtered_outbreaks[1:length(filtered_downtrends)])){browser()}
  
  # Return filtered outbreaks and downtrend dates as a list
  return(list(filtered_outbreaks = filtered_outbreaks, filtered_downtrends = filtered_downtrends))
}

posterior_downtrends <- function(eta_est, upper = 'upper', ...){
  ends <- eta_est$date[eta_est[,upper] < 0]
  return(ends)
}

organize_results_data <- function(stan_summary, stan_data, out_col, date_min = NULL, date_max = NULL, lower = 'X5.', upper = 'X95.', ...){
  
  # get raw data
  df <- stan_data$original 
  
  # filter by date
  if(!is.null(date_min)){
    df <- df %>%
      filter(date >= as.Date(date_min))
  }
  if(!is.null(date_max)){
    df <- df %>%
      filter(date <= as.Date(date_max))
  }
  
  # make the extended df dataset.
  df <- transform_original_data(df, scale_factor = 1, out_col = stan_data$out_col, ...) # set scale factor = 1 to keep original scale.

  # get the eta estimates by fit date
  stan_res <- stan_summary$summary %>% data.frame()
  eta_est <- stan_res[grep('^eta', rownames(stan_res)),] 
  
  if(!(lower %in% colnames(stan_res))){
    stop('lower limit not found in summary data')
  }else if(!(upper %in% colnames(stan_res))){
    stop('upper limit not found in summary data')
  }
  
  # make upper and lower columns.
  eta_est$lower = eta_est[,lower]
  eta_est$upper = eta_est[,upper]
  
  # get the fit dates
  eta_est <- add_eta_date(stan_data, eta_est)
  # index_list <- lapply(rownames(eta_est), function(text){
  #   numbers <- stringr::str_extract_all(text, "\\d+")
  #   num_list <- as.numeric(numbers[[1]])
  #   num_list
  # })
  # fit_index <- sapply(index_list, function(x) x[1]) + stan_data$input$N - stan_data$input$N_estimated
  # fit_id <- match(stan_data$input$T[fit_index], stan_data$original$T)
  # eta_est$date <- stan_data$original$date[fit_id]
  
  return(list(df = df, 
              eta_est = eta_est))
}

add_eta_date <- function(stan_data, eta_est){
  # get the fit dates
  index_list <- lapply(rownames(eta_est), function(text){
    numbers <- stringr::str_extract_all(text, "\\d+")
    num_list <- as.numeric(numbers[[1]])
    num_list
  })
  fit_index <- sapply(index_list, function(x) x[1]) + stan_data$input$N - stan_data$input$N_estimated
  fit_id <- match(stan_data$input$T[fit_index], stan_data$original$T)
  eta_est$date <- stan_data$original$date[fit_id]
  return(eta_est)
}

match_outbreaks <- function(starts, ends, max_date = NULL){
  max_length <- max(length(starts$date), length(ends$date))
  
  df <- data.frame(
    outbreak_start = c(starts$date, rep(NA, max_length - length(starts$date))),
    outbreak_end = c(ends$date, rep(NA, max_length - length(ends$date)))
  )
  
  if(nrow(df) == 0){return(df)}
  
  if(nrow(df) > 1){
    if(any(df$outbreak_end[1:(nrow(df) - 1)] > df$outbreak_start[2:nrow(df)], na.rm = T)){
      browser()
      stop('date mismatch')
    }
    if(any(df$outbreak_end[1:nrow(df) - 1] < df$outbreak_start[1:nrow(df) - 1], na.rm = T)){
      print('shifting back df outbreak ends.')
      # try shifting df end back one.
      df$outbreak_end <- c(df$outbreak_end[2:nrow(df)], NA)
      # check if that works.
      if(any(df$outbreak_end[1:nrow(df) - 1] < df$outbreak_start[1:nrow(df) - 1], na.rm = T) | any(df$outbreak_end[1:(nrow(df) - 1)] > df$outbreak_start[2:nrow(df)], na.rm = T)){
        stop('shifting back didnt work. Outbreak mismatch.')
      }
      # remove rows with all NAs.
      df <- df[!(rowSums(!is.na(df)) == 0),]
    }
  }
  if(!is.null(max_date)){
    df$outbreak_end[is.na(df$outbreak_end)] <- max_date
  }

  return(df)
}

detect_outbreaks <- function(df, eta_est, outbreak_end = 'downtrend', start_date = NULL, end_date = NULL, outbreak_update = NULL, ...){
  # convert date formats.
  eta_est$date <- as.Date(eta_est$date)
  df$date <- as.Date(df$date)
  date_range <- range(c(eta_est$date, df$date))
  
  if(!is.null(start_date)){
    eta_est <- eta_est %>% 
      filter(date >= start_date)
  }
  
  if(!is.null(end_date)){
    eta_est <- eta_est %>% 
      filter(date <= end_date)
  }

  # make outbreaks data set.
  outbreaks <- eta_est$date[eta_est$lower > 0]
  outbreaks_df <- df %>% 
    filter(date %in% outbreaks) %>%
    arrange(date)
  
  if(outbreak_end == 'downtrend'){
    # get the downtrends.
    ends <- posterior_downtrends(eta_est, ...)
    
    # get the outbreaks and downtrends so there is only one in a row.
    tmp <- filter_outbreaks_and_downtrends(outbreaks, ends)
    
    # create the data for plotting.
    outbreaks_df <- df %>% 
      filter(date %in% tmp$filtered_outbreaks) %>%
      arrange(date)
    ends_df <- df %>% 
      filter(date %in% tmp$filtered_downtrends) %>%
      arrange(date)
  }else{
    stop('havent implemented other outbreak end methods.')
  }
  
  # match the outbreak dates.
  matched_outbreaks <- match_outbreaks(outbreaks_df, ends_df, max_date = max(df$date))
  
  if(!is.null(outbreak_update)){
    matched_outbreaks <- update_cases_outbreaks(matched_outbreaks, outbreak_update)
    outbreaks_df <- df %>% 
      filter(date %in% matched_outbreaks$outbreak_start) %>%
      arrange(date)
    ends_df <- df %>% 
      filter(date %in% matched_outbreaks$outbreak_end) %>%
      arrange(date)
  }
  
  return(list(outbreaks_df = outbreaks_df,
              ends_df = ends_df,
              matched_outbreaks = matched_outbreaks))
  
}

