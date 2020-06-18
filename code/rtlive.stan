data {
  int<lower=1> T;
  
  int onset[T];
  vector[T] cumulative_p_delay;
  
  real<lower=0> step_size_prior_sd;
  real<lower=0> theta_init_mean;
  real<lower=0> theta_init_sd;
  
  real<lower=0> si_shape;
  real<lower=0> si_rate;
}

transformed data {
  vector[T-1] inferred_yesterday;
  
  for(t in 1:(T-1)) {
    inferred_yesterday[t] = onset[t] / cumulative_p_delay[t];
  }
}

parameters {
  real<lower=0> step_size;
  real theta_init;
  vector[T-2] theta_steps_unscaled;
  
  real<lower=0> serial_interval;
}

transformed parameters {
  vector[T-1] theta = cumulative_sum(
    append_row(
      rep_vector(theta_init, 1),
      step_size * theta_steps_unscaled
    )
  );
}

model {
  step_size ~ normal(0, step_size_prior_sd);
  theta_init ~ normal(theta_init_mean, theta_init_sd);
  theta_steps_unscaled ~ normal(0, 1);
  
  serial_interval ~ gamma(
    si_shape,
    si_rate
  );
  
  {
    vector[T-1] expected_today = inferred_yesterday
      .* cumulative_p_delay[2:T]
      .* exp(theta);
    
    // Continuous approximation to Poisson as per Li/Dushoff/Bolker
    // The rt.live hack to make mean >= 0.1 MAKES A HUGE DIFFERENCE.
    for(t in 2:T) {
      onset[t] ~ poisson(fmax(expected_today[t-1], 0.1));
    }
  }
}

generated quantities {
  vector[T-1] Rt = serial_interval * theta + 1.0;
}
