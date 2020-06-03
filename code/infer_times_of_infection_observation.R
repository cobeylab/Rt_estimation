## Impute times of observagtion from times of infection


#' This function draws a random sample of size n from the delay distribution (the distribution of times from infection to case observation)
#' Specify the incubation period distribution as gamma(shape = 2, scale = 3) (mean of 6d)
#' Specify the time from symtpom onset to case observation as an exponential distribution with mean 3d
#' These choices don't necessarily reflect COVID biology, but are convenient for testing
#' @param nn is an integer specifying the number of observations to draw
arbitrary_delay_dist <- function(nn){
  r_inc_dist <- function(n){rgamma(n, shape = 2, scale = 3)} # Incubation period (infection -> symptoms)
  r_sym_to_obs_dist <- function(n){rexp(n, 1/3)} # Additional delay from symptoms -> observation
  r_inc_dist(nn) + r_sym_to_obs_dist(nn)
}

SEIR_delay_dist <- function(nn){
  r_inc_dist <- function(n){rgamma(n, shape = 1, scale = 4)} # Incubation period (infection -> symptoms)
  r_sym_to_obs_dist <- function(n){0} # Additional delay from symptoms -> observation
  r_inc_dist(nn) + r_sym_to_obs_dist(nn)
}


#' This function imputes times of case observation from true, underlying times of infection output by the S->E transition in the SEIR model.
#' @param n_dS is a vector of n infections incident on day t
#' @param times is a vector or times, corresponding to the entries of n_dS
#' @param r_delay_dist is a function that draws random samples from the distrubition of time delays form infection to observation
#' #' @param  return_times if true, return a data frame with time. If false, return only a vector of the number observed at each time.
get_tObs_from_tInf <- function(n_dS, 
                               times, 
                               r_delay_dist,
                               return_times = FALSE){
  stopifnot(length(n_dS) == length(times))
  n_dS <- ifelse(is.na(n_dS), 0, n_dS)
  stopifnot(n_dS == round(n_dS))

  
  # This function draws times of observation for each of the ndS cases incident at a given time point
  get_obs_times_for_one_timestep <- function(ndS, tt){
    (tt - 
       runif(ndS) +  ## Subtract a uniform between 0 and 1 to get exact time of day of onset on the previous day
       r_delay_dist(ndS))  %>% 
      ceiling()  # Round up: infections are recorded at the end of the day in progress
  } 
  
  ## Get a vector of imputed observation times for each infection
  obs_time_vec <- mapply(FUN = get_obs_times_for_one_timestep, ndS = n_dS, tt = times) %>%
    unlist() 
  stopifnot(length(obs_time_vec) == sum(n_dS))
  
  ## Reformat: count the number of observed infections at each time, and output
  data.frame(time = obs_time_vec) %>%
    group_by(time) %>%
    summarise(n = n()) %>%
    # Pad with 0s at times with no observed cases
    complete(time = times, fill = list(n = 0)) -> out

  if(return_times) out else out$n
}


#' This function imputes times of infection from times of observation, given an known delay distribution
#' @param n_obs is a vector of n infections newly observed on day t
#' @param times is a vector or times, corresponding to the entries of n_obs
#' @param r_delay_dist is a function that draws random samples from the distrubition of time delays form infection to observation
#' 
# n_obs <- obs_df$imputed_outpatient[1:mt]
# times = obs_df$time[1:mt]
# r_delay_dist = outpatient_delay_dist

get_tInf_from_tObs <- function(n_obs,
                               times,
                               r_delay_dist,
                               adjust_right_censor = T){
  stopifnot(length(n_obs) == length(times))
  n_obs <- ifelse(is.na(n_obs), 0, n_obs)
  stopifnot(n_obs == round(n_obs))

  set.seed(33)
  
  ## Draw times of infection for each of the n_obs cases observed at time tt
  get_inftimes_for_one_timestep <- function(nn, tt){
    (tt  - 
       runif(nn) -
       r_delay_dist(nn))  %>%
      ceiling()
  } 
  
  ## Get a vector of observation times for each infection
  inf_time_vec <- mapply(FUN = get_inftimes_for_one_timestep, nn = n_obs, tt = times) %>%
    unlist() 
  
  ## Reformat: count the number of observations at each timestep
  data.frame(inf_time = inf_time_vec) %>%
    group_by(inf_time) %>%
    summarise(n = n()) %>%
    ## Drop imputed infection times <0
    filter(inf_time >= min(times)) %>%
    # Pad with 0s at times with no observed cases
    complete(inf_time = times, fill = list(n = 0)) -> raw_df
  
  if(adjust_right_censor){
  ## Adjust for right censoring
  #   xx = 10000
  # delay_sample = r_delay_dist(xx) ## Get the cumulative probability that an infection occurring n days ago is included in the data, given the maximum possible date of observation
  # cum_p_delay <- sapply((rev(times)-min(times)+1), function(tt){sum(delay_sample< tt)/xx})
  cum_p_delay <- get_cum_p_delay(rdd = r_delay_dist, mt = max(times))
  raw_df$n <- ifelse(cum_p_delay == 0, 0, raw_df$n/cum_p_delay)
  #raw_df$n <- ifelse(is.na(raw_df$n), 0, raw_df$n)
  }
  raw_df %>% pull(n)
}


move_back_in_time <- function(n_obs,
                              mean_delay){
  c(
    n_obs[ceiling(mean_delay):length(n_obs)],
    rep(NA, floor(mean_delay))
  ) 
}

get_cum_p_delay <- function(rdd, mt){
  get_tObs_from_tInf(n_dS = 10000, times = 0, r_delay_dist = rdd, return_times = T) %>% 
    mutate(cump = cumsum(n)/sum(n)) %>% 
    complete(time = 0:mt, fill = list(cump = 1)) %>%
    pull(cump) %>% 
    rev()
  
}
