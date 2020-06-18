## Functions to calculate exact Rt
## Code by Ed Baskerville
## June 13, 2020

## For an explanation of the calculations shown here, see Rc_math.Rmd
library(deSolve)



## Function to compute p(u) the probability that an individual is infectious u days after the moment of infection.
make_p_infectious <- function(sigma, ## rate of leaving compartment E
                              gamma  ## rate of leaving compartment I
                              ) {
  if(sigma == gamma) {
    function(s) {
      gamma * s * exp(-gamma * s)
    }
  } else {
    function(s) {
      sigma * (exp(-sigma * s) - exp(-gamma * s)) / (gamma - sigma)
    }
  }
}


## Function to integrate SEIR
integrate_seir <- function(beta_t, ## Function that returns beta (t) for any time t
                           sigma,  ## Rate of leaving latent class 
                           gamma,  ## Rate of leaving infected class
                           N,      ## Total population size 
                           T,      ## Final timestep
                           E0,     ## Initial exposed
                           I0,     ## Initial infeced
                           dt = 0.01 ## Size of timestep over which to integrate SEIR
                           ) {
  S0 <- N - E0 - I0
  ddt_ode <- function(t, y, params) {
    S <- y[1]
    E <- y[2]
    I <- y[3]
    
    dSE <- beta_t(t) * S * I / N
    dEI <- sigma * E
    dIR <- gamma * I
    
    list(c(S = -dSE, E = dSE - dEI, I = dEI - dIR))
  }
  
  as.data.frame(ode(
    c(
      S = S0, E = E0, I = I0
    ), seq(0, T, dt), ddt_ode, list(), method = 'rk4'
  ))
}


# # Example
# arnaught <- 2
# t_E <- 2
# t_I <- 4
# beta_t <- function(t) { arnaught / t_I } # Can replace this with a beta(t) with changes
# sigma <- 1 / t_E
# gamma <- 1 / t_I
# N <- 1
# T <- 80
# E0 <- 0.001
# I0 <- 0
# ode_output <- as.data.frame(integrate_seir(beta_t, sigma, gamma, N, T, E0, 0))
# ggplot(ode_output, aes(x = time, y = S / N)) + geom_line()



## Calculate the instaneous case reproductive number using the deterministic SEIR output
integrate_Rt <- function(beta_t,  ## Function that returns beta (t) for any time t
                         sigma,   ## Rate of leaving latent class 
                         gamma,   ## Rate of leaving infected class
                         N,       ## Total population size
                         T,       ## Final time
                         E0,      ## Initial latent
                         I0,      ## Initial infected
                         dt_int = 0.1,  ## Integration timesetp
                         dt_Rt = 1,     ## Rc output timestep
                         T_inf = 1000   ## Final time in SEIR integration. Should be at least one generation interval higher than T
                         ) {
  S0 <- (N - E0 - I0)
  T_inf = T+100
  if(any(is.na(beta_t(0:T_inf)))){stop('beta_t must return values from 0:(T+100). Either decrease T or increase the range of beta_t.')}
  
  ode_output <- as.data.frame(
    integrate_seir(beta_t, sigma, gamma, N, T_inf, E0, I0)
  )
  SoverN <- approxfun(ode_output$time, ode_output$S / N)
  
  p_infectious <- make_p_infectious(sigma, gamma)
  
  times <- seq(0, T, dt_Rt)
  Rt <- sapply(times, function(t) {
    integrate(function(u) {
      beta_t(u) * p_infectious(u - t) * SoverN(u)
    }, t, T_inf)$value
  })
  
  data.frame(
    time = times,
    R_case = Rt
  )
}


# ## Test - See how it looks (dotted line is the instantaneous reproductive number):
# Rt_output <- integrate_Rt(beta_t, sigma, gamma, N, 300, E0, I0)
# ggplot() +
#   geom_line(data = Rt_output, aes(x = time, y = Rt)) +
#   geom_line(
#     data = ode_output %>%
#       mutate(Rt = beta_t(time) / gamma * S / N),
#     aes(x = time, y = Rt),
#     linetype = 2
#   )


