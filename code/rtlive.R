rtlive_pymc3 <- function(
  python_path,
  onset, # Time series of onset cases
  cumulative_p_delay, # Adjustment for right-censoring
  step_size_prior_sd, # Prior for theta random walk, original uses 0.03
  theta_init_mean, theta_init_sd, # Prior for initial theta
  si_shape, si_rate,
  chains = 1, tune = 2000, draws = 2000
) {
  library(jsonlite)
  
  input <- toJSON(list(
    model_data = list(
      onset = onset,
      cumulative_p_delay = cumulative_p_delay,
      step_size_prior_sd = unbox(step_size_prior_sd),
      theta_init_mean = unbox(theta_init_mean),
      theta_init_sd = unbox(theta_init_sd),
      si_shape = unbox(si_shape),
      si_rate = unbox(si_rate)
    ),
    mcmc_params = list(
      chains = unbox(chains),
      tune = unbox(tune),
      draws = unbox(draws)
    )
  ))
  
  json_output <- system(
    sprintf('%s rtlive.py', python_path),
    input = input,
    intern = TRUE
  )
  
  fromJSON(json_output)
}

rtlive_stan <- function(
  mod = 'rtlive.stan',
  onset, # Time series of onset cases
  cumulative_p_delay, # Adjustment for right-censoring
  step_size_prior_sd, # Prior for theta random walk, original uses 0.03
  theta_init_mean, theta_init_sd, # Prior for initial theta
  si_shape, si_rate,
  chains = 1, # Number of MCMC chains to run
  cores = parallel::detectCores(logical = FALSE),
  ...
) {
  library(rstan)
  
  model <- stan_model(mod)
  
  T <- length(onset)
  model_input_data <- list(
    T = T,
    
    onset = onset,
    cumulative_p_delay = cumulative_p_delay,
    
    step_size_prior_sd = step_size_prior_sd,
    theta_init_mean = theta_init_mean,
    theta_init_sd = theta_init_sd,
    si_shape = si_shape,
    si_rate = si_rate
  )
  
  fit <- sampling(
    model, model_input_data, chains = chains, cores = cores,
    ...
  )
}
