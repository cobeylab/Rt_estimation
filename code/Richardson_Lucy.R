## This function inputs observed incidence (numbers of cases, deaths, hospitalizations, etc. per day), and uses an adaptation of  Richardson-Lucy deconvolution to recover the underlying incidence curve. 

## ------- Inputs
##  This function follows methods in Goldstein et al. (https://www.pnas.org/content/pnas/106/51/21825.full.pdf)
#' @param observed - a vector of the number of observed cases, deaths, hospitalizations, etc. per day
#' @param times - a numeric vector of times corresponding to the entries in observed. Must be the same length as observed.
#' @param p_delay - a numeric vector whose entries give the probability that the delay from infection to observation is exactly 0, 1, ... N days. p_delay must sum to 1 and to to facilitate vectorization should be the same length as observed, and times (the last several entries will most likely be 0s).
#' @param max_iter - maximum number of times to iterate. Following Goldstien et al., the algorithm continues to run until the normalized chi2 statistic comparing the observed and expected number of deaths per day falls below 1, or until the maximum number of iterations is reached.
#' @param out_col_name - string giving the name of the column in which to output the imputed times of infection
## --------- Outputs
#' @value - the function returns a data frame with columns time and out_col_name (which gives the imputed number of infections per timestep)



get_RL <- function(observed, ## Vector of observed cases. Let L = length(observed)
                   times,    ## Vector of times of observations. Must be length L.
                   p_delay,  ## Vector of probabilities that the delay from infection to observation is 0, 1, ... days. Can be any length <= L.
                   max_iter = 50,
                   out_col_name = 'RL_result',
                   right_censor = TRUE, # If TRUE, upscale the observations at the end of the time series based on the cumulative probability that an infection occurring on that date would have been observed.
                   verbose = TRUE){
  
  ## Check inputs
  stopifnot(is.vector(observed) & is.vector(times) & is.vector(p_delay))
  stopifnot(length(p_delay) <= length(observed))
  if(length(p_delay)<length(observed)){
    p_delay <- c(p_delay, rep(0, length(observed)-length(p_delay))) ## Pad the end of p_delay with 0s so it's the same length as the observations
  }
  stopifnot(length(observed) == length(times) & length(observed) == length(p_delay))
  if(abs(sum(p_delay)-1) > 1e-6 & verbose){warning('p delay vector does not sum to 1. Extened max delay?')}
  observed <- na_to_0(observed)
  
  
  ## From the delay multinomial, extract the mean and max delay from infection to observation
  mean_delay <- sum(p_delay * 0:(length(p_delay)-1))
  max_delay <- tail((0:length(p_delay))[p_delay>0], 1)
  orig_min_t <- min(times)
  orig_ceiling <- observed[1]*10
  
  ## Get the cumulative probability of a delay <= x days, which we'll use later to adjust for right censoring
  if(right_censor){
    cum_p_delay <- cumsum(p_delay)
  }else{
    cum_p_delay <- p_delay*0+1
  }
  
  ## Infections could have occured up to max_delay days prior to the start of the observed time series. Pad the beginning of the observed time   series and the times vector with zeros as appropriate.
  observed <- c(
    rep(0, max_delay), ## Pad with 0s @ beginning
    observed
  )
  times <- c(
    min(times)-(max_delay:1), ## Insert times prior to first observed
    times
  )
  stopifnot(length(observed) == length(times))
  
  
  ## The initial guess for lambda is the observed time series, shifted back by the mean dealy to observation
  lambda <- c(
    observed[ceiling(mean_delay):length(observed)],  ## Observed vector, with first mean_delay days removed
    rep(tail(observed, 1), floor(mean_delay)) ## Replace discarded values with 1s at the end of the vector
  )
  
  # Write a function to get the expected deaths per day, as a function of lambda:
  #  Note that by convention, p_delay[1] corresponds to a 0-day delay, p_delay[2] corresponds to a 1-day delay, etc.
  get_expected_D_i <- function(ii, lambda, p_delay){
    sum(p_delay[1:ii]*lambda[ii:1], na.rm = TRUE)
  }
  
  # Write a function to get the updated value of lambda_j, given the previous value
  get_lambda_j_update <- function(jj, lambda, p_delay, obs, expected){
    N = length(obs)
    lambda[jj]/cum_p_delay[(N-jj+1)] * sum(
      obs[jj:N]/expected[jj:N] * p_delay[1:(N-jj+1)],
      na.rm = TRUE
    )
  }
  
  ## Stopping condition: Write a funtion to get the chi2 statistic comparing observed and expected deaths per day
  get_chisq <- function(obs, exp){
    ## If right censoring, drop dates on which the probability of observation is <50% from the stopping condition
    drop_me <- sum(cum_p_delay < 0.5)
    nn = length(obs)-drop_me
    obs = obs[1:nn]
    exp = exp[1:nn]
    ## Get chi2 statistic using reliable values
    1/length(obs) * sum(ifelse(exp == 0, 1, (obs-exp)^2/(exp)))
  }
  
  ## Set an intial value for the expected vector
  expected_D <- sapply(1:(length(lambda)), get_expected_D_i, lambda = lambda, p_delay = p_delay)
  ## Iterate to solve for lambda (the inferred time series of infections)
  iter = 1
  while(get_chisq(observed, expected_D) > 1 & iter < max_iter){
    expected_D <- sapply(1:(length(lambda)), get_expected_D_i, lambda = lambda, p_delay = p_delay)
    lambda <- sapply(1:(length(lambda)), get_lambda_j_update,
                     lambda = lambda, 
                     p_delay = p_delay, 
                     obs = observed, 
                     expected = expected_D)
    iter = iter+1
  }
  ## Clean
  clean_lambda <- function(xx){ifelse(xx>orig_ceiling, NA, xx)}
  lambda[times<=orig_min_t] = clean_lambda(lambda[times<=orig_min_t])
  ## Return
  if(verbose){print(sprintf('Returning RL after %.0f iterations', iter))}
  data.frame(time = times, imputed = lambda) %>% setNames(c('time', out_col_name))
}





## From supplement of Cori et al.
get_discrete_lognormal <- function(k, mu, sigma){
  ff <- function(xx, mu, sigma){xx*dlnorm(xx, mu, sigma)}
  sapply(X = k, FUN = function(k, mu, sigma){
    (1+k)*plnorm(k+1, mu, sigma) -
      2*k*plnorm(k, mu, sigma) +
      (k-1)*plnorm(k-1, mu, sigma)+
      integrate(f = ff, lower = k-1, upper = k, mu = mu, sigma = sigma)$value-
      integrate(f = ff, lower = k, upper = k+1, mu = mu, sigma = sigma)$value
  }, mu = mu, sigma = sigma)
}