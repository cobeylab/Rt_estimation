#' is.wholenumber
#'
#' For parameter checking
is.wholenumber <- function(x) {
  x == round(x)
}

check_args <- function(
  arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t, n_steps_per_t
) {
  # Check all input parameters
  stopifnot(is.wholenumber(N) && length(N) == 1 && N >= 1)
  stopifnot(is.wholenumber(n_t) && length(n_t) == 1 && n_t >= 1)
  stopifnot(
    is.wholenumber(n_steps_per_t) && length(n_steps_per_t) == 1 &&
      n_steps_per_t >= 1
  )
  stopifnot(
    is.numeric(arnaught) && arnaught > 0 &&
      (length(arnaught) == 1 || length(arnaught) == n_t + 1)
  )
  stopifnot(
    is.numeric(t_E) && length(t_E) == 1 && t_E >= 0
  )
  stopifnot(
    is.numeric(t_I) && length(t_I) == 1 && t_I > 0
  )
}

construct_beta <- function(arnaught, t_I, n_t) {
  beta_t_all <- arnaught / t_I
  if(length(arnaught) == 1) {
    function(t) beta_t_all
  } else {
    approxfun(0:n_t, beta_t_all)
  }
}

simulate_seir <- function(
  arnaught, t_E, t_I,
  N, S_init, E_init, I_init,
  n_t, n_steps_per_t = 1,
  method = 'stochastic'
) {
  func <- if(method == 'stochastic') simulate_seir_stochastic else simulate_seir_ode
  func(
    arnaught, t_E, t_I,
    N, S_init, E_init, I_init,
    n_t, n_steps_per_t
  )
}

#' Simulate a discrete-time approximation of a continuous-time, discrete-state stochastic S(E)IR model.
#' 
#' Ed Baskerville
#' 15 April 2020
#' 
#' No age structure.
#' 
#' @param n_t Number of units of time to simulate.
#' @param n_steps_per_t Number of discrete simulation steps to take per
#' unit of simulation time.
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0. Average number of new infections produced by an infection in a susceptible population. A scalar or a vector of length `n_t + 1`, which specifies R0 at the start of each timestep. R0 is linearly interpolated between timesteps.
#' @param t_E Mean latent period. If set to 0, the model reduces to an SIR.
#' @param t_I Mean duration of infectiousness.
# # for debugging
# arnaught
# N = N
# E_init = 0          # Initial in E
# I_init = 1e-5 * N     # Initial in I
# t_E = 3               # time from exposed to infected
# t_I = 5               # time from infected to recovery
# n_t = 1000            # total simulation time
# max_R0 = 2.0          # Initial, max R0
# min_R0 = 0.9         # Minimum with intervention
# end_max_time = 50     # time of intervention
# start_min_time = 70  #time when R0 hits lowest value
# n_steps_per_t = 10
# CONTINUOUS = TRUE       # If TRUE, R decreases according to arctan after intervention. If FALSE, R


simulate_seir_stochastic <- function(
  arnaught, 
  t_E, 
  t_I,
  N, 
  S_init, 
  E_init, 
  I_init,
  n_t, 
  n_steps_per_t
) {
  check_args(
    arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t, n_steps_per_t
  ) 
  
  # Precompute a few things
  delta_t <- 1 / n_steps_per_t
  
  # Draws a binomial based on a rate
  draw <- function(n, rate) {
    p <- 1 - exp(-rate * delta_t)
    rbinom(1, n, p)
  }
  
  # Function to compute beta at a particular time
  beta <- construct_beta(arnaught, t_I, n_t)
  
  # Step forward from t to t + delta_t
  step <- function(t, S_prev, E_prev, I_prev) {
    dS <- draw(S_prev, beta(t) * I_prev / N)
    dIR <- draw(I_prev, 1 / t_I)
    
    if(t_E > 0) {
      # SEIR model
      dEI <- draw(E_prev, 1 / t_E)
      list(
        S = S_prev - dS,
        E = E_prev + dS - dEI,
        I = I_prev + dEI - dIR,
        dS = dS,
        dEI = dEI,
        dIR = dIR
      )
    }
    else {
      # SIR model
      list(
        S = S_prev - dS,
        E = 0,
        I = I_prev + dS - dIR,
        dS = dS,
        dEI = dS,
        dIR = dIR
      )
    }
  }
  
  # Set up state vectors over time
  S <- numeric(n_t + 1)
  S[1] <- S_init
  
  E <- numeric(n_t + 1)
  I <- numeric(n_t + 1)
  if(t_E > 0) {
    # SEIR model
    E[1] <- E_init
    I[1] <- I_init
  }else {
    # SIR model: all initial E get dumped in I
    E[1] <- 0
    I[1] <- E_init + I_init
  }
  
  # Track transitions over time
  dS <- rep(NA, n_t + 1)
  dEI <- rep(NA, n_t + 1)
  dIR <- rep(NA, n_t + 1)
  
  # Simulate
  for(tt in 1:(n_t-1)) {
    S_prev <- S[tt]
    E_prev <- E[tt]
    I_prev <- I[tt]
    
    dS[tt+1] <- 0
    dEI[tt+1] <- 0
    dIR[tt+1] <- 0
    for(i in 1:n_steps_per_t) {
      state_next <- step(tt + delta_t * (i - 1), S_prev, E_prev, I_prev)
      S_prev <- state_next$S
      E_prev <- state_next$E
      I_prev <- state_next$I
      dS[tt+1] <- dS[tt+1] + state_next$dS
      dEI[tt+1] <- dEI[tt+1] + state_next$dEI
      dIR[tt+1] <- dIR[tt+1] + state_next$dIR
    }
    
    S[tt+1] <- S_prev
    E[tt+1] <- E_prev
    I[tt+1] <- I_prev
  }
  
  # Return each compartment over time
  data.frame(
    time = 0:n_t,
    S = S,
    E = E,
    I = I,
    R = N - S - E - I,
    dS = dS,
    dEI = dEI,
    dIR = dIR
  )
}

#' Simulate a deterministic ODE approximation of a continuous-time, discrete-state stochastic S(E)IR model.
#' 
#' Ed Baskerville
#' 15 April 2020
#' 
#' No age structure.
#' 
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0. Average number of new infections produced by an infection in a susceptible population. A scalar or a vector of length `n_t + 1`, which specifies R0 at the start of each timestep. R0 is linearly interpolated between timesteps.
#' @param t_E Mean latent period. If set to 0, the model reduces to an SIR.
#' @param t_I Mean duration of infectiousness.
#' @param n_t Number of units of time to simulate.
simulate_seir_ode <- function(
  arnaught, t_E, t_I,
  N, S_init, E_init, I_init,
  n_t,
  n_steps_per_t = 1 # Ignored; included so the function signature matches stochastic version
) {
  library(deSolve)
  
  check_args(
    arnaught, t_E, t_I, N, S_init, E_init, I_init, n_t, n_steps_per_t
  )
  
  beta <- construct_beta(arnaught, t_I, n_t)
  d_dt <- function(t, y, params) {
    dS <- y['S'] * beta(t) * y['I'] / N
    dIR <- y['I'] / t_I
    
    if(t_E > 0) {
      # SEIR model
      dEI <- y['E'] / t_E
      list(c(
        S = -dS,
        E = dS - dEI,
        I = dEI - dIR,
        R = dIR,
        cum_dS = dS,
        cum_dEI = dEI
      ), NULL)
    }
    else {
      # SIR model
      list(c(
        S = -dS,
        E = 0,
        I = dS - dIR,
        R = dIR,
        cum_dS = dS,
        cum_dEI = dS
      ), NULL)
    }
  }
  
  y_init <- c(
    S = S_init,
    E = if(t_E > 0) E_init else 0,
    I = if(t_E > 0) I_init else E_init + I_init,
    R = 0,
    cum_dS = 0,
    cum_dEI = 0
  )
  #automatic ode solver is lsoda, an "automatic stiff/non-stiff solver"
  as.data.frame(ode(y_init, 0:n_t, d_dt, NULL)) %>%
    mutate(dS = cum_dS - lag(cum_dS, 1)) %>%
    mutate(dEI = cum_dEI - lag(cum_dEI, 1)) %>%
    mutate(dIR = R - lag(R, 1))
}




# ## Test
# simulate_seir_stochastic(arnaught = 2.0, t_E = 4, t_I = 4, N = 1e6, S_init = N-10, E_init = 0, I_init = 10, n_t = 1000, n_steps_per_t = 10) -> stochastic
# simulate_seir_ode(arnaught = 2.0, t_E = 4, t_I = 4, N = 1e6, S_init = N-10, E_init = 0, I_init = 10, n_t = 1000, n_steps_per_t = 10) -> ode
# 
# bind_rows(stochastic, ode, .id = c('method')) %>%
#   mutate(method = factor(method, labels = c('stochastic', 'ode'))) %>%
#   pivot_longer(S:dIR) %>%
#   ggplot()+
#   geom_line(aes(x = time, y = value, color = method)) +
#   facet_wrap(.~name, scales = 'free_y') +
# #   xlim(c(0, 200))
# 
# 
# lapply(1:1000, function(xx) simulate_seir_example(arnaught = 2.0, t_E = 4, t_I = 4, N = 1e6, E_init = 0, I_init = 10, n_t = 200, n_steps_per_t = 10, method = 'stochastic')) %>%
#   bind_rows() %>%
#   group_by(time) %>%
#   summarise_all(mean) -> stochastic
#   
# simulate_seir_example(arnaught = 2.0, t_E = 4, t_I = 4, N = 1e6, E_init = 0, I_init = 10, n_t = 200, n_steps_per_t = 10, method = 'ode') -> ode
# bind_rows(stochastic, ode, .id = c('method')) %>%
#   mutate(method = factor(method, labels = c('stochastic', 'ode'))) %>%
#   pivot_longer(S:dIR) %>%
#   ggplot()+
#   geom_line(aes(x = time, y = value, color = method)) +
#   facet_wrap(.~name, scales = 'free_y') +
#   xlim(c(0, 200))
