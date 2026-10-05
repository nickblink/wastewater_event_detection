data {
  int<lower=1> N; // Number of data points.
  int<lower=1> N_estimated; // Number of data points the model is fit on.
  int<lower=1> window; // Time window in days for fitting exponential growth.
  vector[N] T; // Time points in days, starting at 1.
  array[N] real Y; // Outcomes at each time point.
  array[N] int<lower=0> Y_int; // Integer version of Y
  int<lower=0, upper=3> family; // 1 = Normal (constant variance), 2 = Poisson, 3 = Negative Binomial.
  int<lower=0,upper=1> constant_family_param; // 1 = keep family parameter constant. 0 = independently estimate.
  real theta_prior_shape; // Prior value for theta shape.
  real theta_prior_rate; // Prior value for theta rate.
  //real sigma2_prior_shape; // Prior value for sigma2 shape.
  //real sigma2_prior_rate; // Prior value for sigma2 rate.
  real sigma_prior_sd; // Prior standard deviation on the sigma.
  real beta_prior_shape; // Prior value for beta shape.
  real beta_prior_rate; // Prior value for beta rate.
  real sigma_eta; // Standard deviation on the eta prior.
  int<lower=0,upper=1> debug; // Whether to print out values for the sake of debugging.
}
transformed data {
  array[N] int start_idx; // Precompute the start index for each window
  int not_est; // Number of not estimated points.
  int integer_outcome; // 0 = real, continuous outcome. 1 = integer outcome.
  row_vector[window + 1] predict_window; // window to estimate on.
  for (n in 1:N) {
    int idx = 1;
    while (T[n] > window && T[idx] < T[n] - window) {
      idx += 1;
    }
    start_idx[n] = idx;
  }
  not_est = N - N_estimated;
  for(i in 0:window){
    predict_window[i+1] = i;
  }
  if(family == 1){
    integer_outcome = 0;
  }else{
    integer_outcome = 1;
  }
}
parameters{
  array[constant_family_param ? 1 : N_estimated] real<lower=0> theta; // Negative binomial dispersion.
  array[constant_family_param ? 1 : N_estimated] real<lower=0> sigma; // Standard deviation of normal.
  array[N_estimated] real<lower=0> beta; // The intercept at each time window.
  vector[N_estimated] eta; // The growth rate at each time window.
}
transformed parameters{
  array[N_estimated] real<lower=0> theta_vec; // theta_vec used.
  array[N_estimated] real<lower=0> sigma_vec; // sigma_vec used.
  
  // make theta vec
  if(constant_family_param){
    for(t in 1:N_estimated){
	  theta_vec[t] = theta[1];
	}
  }else{
    theta_vec = theta;
  }
  
  // make sigma vec
  if(constant_family_param){
    for(t in 1:N_estimated){
	  sigma_vec[t] = sigma[1];
	}
  }else{
    sigma_vec = sigma;
  }
}
model{
  // priors 
  theta ~ gamma(theta_prior_shape, theta_prior_rate);
  //sigma2 ~ gamma(sigma2_prior_shape, sigma2_prior_rate);
  sigma ~ normal(0, sigma_prior_sd); 
  //beta ~ normal(Y[not_est + 1:N], sigma2_OLD); // centering each intercept at the initial value of that time series.
  beta ~ gamma(beta_prior_shape, beta_prior_rate);
  eta ~ normal(0, sigma_eta);

  // likelihood
  for (n in 1:N) {
    int start = start_idx[n];
    int count = n - start + 1; // Number of data points in the window
	int t = n - not_est; // 
    
	// only fit model if past initial window and at least one datapoint in window.
    if (T[n] > window && count > 1) {
	  // debugging
	  if(t < 0){
	    print("n = ", n, ", t = ", t, ", not_est = ", not_est, ", count = ", count);  // Debug output
	  }
	  if(beta[t] < 0){
	    print("beta[t] = ", beta[t]);
	  }
	  
      vector[count] T_window =  T[start:n] - T[n] + window; // Counting indices from the start of the window.
      //real Y_window[count] = Y[start:n];
	  
      // Run the model for the current window of data (vectorized)
	  vector[count] mu_t = beta[t]*exp(eta[t]*T_window) + 1; // Get the mean
	  
	  if(family == 1){ // normal
	    array[count] real Y_window = Y[start:n];
	    Y_window ~ normal(mu_t, sigma_vec[t]);
	  }else if(family == 2){ // poisson
	    array[count] int Y_window = Y_int[start:n];
	    Y_window ~ poisson(mu_t); // the intended model.
	  }else if(family == 3){ // negative binomial
	    array[count] int Y_window = Y_int[start:n];
	    Y_window ~ neg_binomial_2(mu_t, theta_vec[t]); // the intended model.
		if(debug == 1){
		  print("theta = ", theta); 
		}
	  }
    }
  } 
}
generated quantities{
  matrix[N_estimated, window + 1] estimates; // matrix of estimates. Rows are the starting points for estimation. Columns are the predictions at each time point for that fit.
  matrix[N_estimated*(1-integer_outcome), (window + 1)*(1-integer_outcome)] predictions; // matrix of predictions for real, continous outcomes. Dimensions are [0,0] if using integer outcomes.
  array[N_estimated*integer_outcome, (window + 1)*integer_outcome] int predictions_int; // matrix of predictions for integer outcomes.
  for (t in 1:N_estimated){
    estimates[t] = beta[t]*exp(eta[t]*predict_window);
	if(family == 1){
	  predictions[t] = to_row_vector(normal_rng(estimates[t], sigma_vec[t]));
	}else if(family == 2){
	  predictions_int[t] = poisson_rng(estimates[t]);
	}else if(family == 3){
	  predictions_int[t] = neg_binomial_2_rng(estimates[t], theta_vec[t]);
	}
  }
}
