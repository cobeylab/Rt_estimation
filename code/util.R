## Utility functions for simulating data and estimating R0


## For generating synthetic data ----------------------------------------

## Output underlying R0 values for input into the GIR or SEIR model across once change (increase or decrease)
get_arnaught_step <- function(start_R0, 
                                  final_R0, 
                                  start_change,  ## Time at which R0 first departs from its initial value, start_R0
                                  end_change,   ## Time at which R0 first reaches its final value, final_R0
                                  n_t){  ## total timesteps
  arnaught <- if(start_R0==final_R0){
    #constant arnaught
    start_R0
  } else if(start_change==end_change){
    #step function arnaught
    swap_time <- start_change + 1
    c(rep(start_R0, swap_time), rep(final_R0, n_t-swap_time))

  } else {
    #linearly decreasing intervention arnaught
    c(rep(start_R0, start_change), seq(start_R0, final_R0, length.out=(end_change-start_change+1)), rep(final_R0, n_t - end_change))
  } 
}

## Repeatedly call get_arnaught step to specify R0 across an arbitrary number of changes
specify_arnaught <- function(R0_vec, ## Vector of equlibrium R0 values. Should be 1 greater in length than the desired number of changes.
                             change_start_vec, ## Vector of timepoints at which R0 starts to change. Should be length(R0_vec) - 1
                             change_end_vec, ## Vector of timepoints at which R0 first reaches its new value. Should be length(R0_vec) - 1
                             NT){ ## Scalar - total number of timepoints.
  stopifnot(length(change_end_vec) == length(change_start_vec))
  stopifnot(length(change_start_vec) == length(R0_vec)-1)
  stopifnot(all(diff(change_start_vec)>0))
  stopifnot(all(diff(change_end_vec)>0))
  arnaught <- NULL
  n.changes <- length(R0_vec)-1
  breakpoints <- c(0, change_start_vec[-1]-1, NT+1)
  stopifnot(all(breakpoints[-1]>change_end_vec))
  for(ii in 1:n.changes){
    arnaught <- c(arnaught,
                  get_arnaught_step(start_R0 = R0_vec[ii], 
                                  final_R0 = R0_vec[ii+1], 
                                  start_change = change_start_vec[ii]-breakpoints[ii], 
                                  end_change = change_end_vec[ii]-breakpoints[ii], 
                                  n_t = diff(breakpoints)[ii]-1)
    )
  }
  #cbind(arnaught, 0:NT)
  arnaught
}
# # Test 1
# specify_arnaught(R0_vec = c(2.5, .7, 1.1), change_start_vec = 45+c(0, 45), change_end_vec = 45+c(0, 45)+7, NT = 150)
# 
# ## Test 2
# specify_arnaught(R0_vec = c(2.5, .7, 1.1, 2.0), change_start_vec = 45+c(0, 45, 45+30), change_end_vec = 45+c(0, 45, 45+30)+7, NT = 150) -> test
# plot(0:150, test)
# abline(v = c(45+c(0, 45, 45+30)), lty = 2, col = 'red')
# abline(v = c(45+7+c(0, 45, 45+30)), lty = 2, col = 'red')
# 
# # Test 3 - decreasing change points, should throw error
# specify_arnaught(R0_vec = c(2.5, .7, 1.1, 2.0), change_start_vec = 45+c(0, 45, 30), change_end_vec = 45+c(0, 45, 45+30)+7, NT = 150)
# 
# #Test 4 - endpoint is after breakpoint - should throw error
# specify_arnaught(R0_vec = c(2.5, .7, 1.1, 2.0), change_start_vec = 45+c(0, 45, 90), change_end_vec = 45+c(0, 45, 45+30)+50, NT = 150)



## Wrappers that call simulate_sir defined in simulation.R
simulate_sir_example <- function(
  arnaught, t_I, N, I_init, n_t, n_steps_per_t = 10,
  method = 'stochastic'
) {
  simulate_seir(
    arnaught = arnaught,
    t_E = 0,
    t_I = t_I,
    N = N,
    S_init = N - I_init,
    E_init = 0,
    I_init = I_init,
    n_t = n_t, n_steps_per_t = n_steps_per_t,
    method = method
  )
}


## Wrappers that call simulate_seir defined in simulation.R
simulate_seir_example <- function(
  arnaught, t_E, t_I, N, E_init, I_init, n_t, n_steps_per_t = 10,
  method = 'stochastic'
) {
  simulate_seir(
    arnaught = arnaught,
    t_E = t_E,
    t_I = t_I,
    N = N,
    S_init = N - E_init - I_init,
    E_init = E_init,
    I_init = I_init,
    n_t = n_t, n_steps_per_t = n_steps_per_t,
    method = method
  )
}


## Function to load saved outputs form a simulation run
load_sims_for_one_R0 <-  function(arnaught, model_type = 'seir', method = 'stochastic'){
  data.frame(fns = list.files(path = sprintf('R0-%.1f', arnaught))) %>%
    filter(grepl(fns, pattern = method)) %>%
    pull(fns) %>%
    lapply(FUN = function(xx){
      readRDS(paste0(sprintf('R0-%.1f/', arnaught), xx)) -> tmp
      tmp$sim_df %>% 
        mutate(int_time = as.character(tmp$intervention_time), 
               dec_dur = as.character(tmp$decrease_duration))
    }) %>%
    bind_rows()
}



## Write a function to extract the simulation results as a data frame
get_sim_df <- function(method = 'ode'){  ## can also be 'stochastic'
  readRDS(sprintf('R0-%.1f/seir_%s_dec%.0f-%.0f_sim.rds', 
                  parlist$pre_intervention_R0, 
                  method,
                  parlist$intervention_time_1, 
                  parlist$days_intervention_to_min))$sim_df %>%
    mutate_all(.funs = function(xx){ifelse(is.na(xx), 0, xx)}) %>%
    mutate(incidence = round(dS))
}


## Function to replace NAs with 0s in simulation output
na_to_0 <- function(vec){
  if(any(is.na(vec))){
    warning(sprintf('WARNING: Replacing NAs in %s with 0s\n', deparse(substitute(vec))))
    vec[is.na(vec)] = 0
  }
  vec
}






## For estimating Rt ----------------------------------------

## Output cori estimate with mean, CI and times given an input df, and the name of the incidence column
get_WT <- function(df.in, 
                   icol_name, 
                   outcol_name = 'WT',
                   window = 1, 
                   GI_mean=parlist$true_mean_GI, 
                   GI_var=2*(parlist$true_mean_GI/2)^2,
                   wend = FALSE){
  idat <- df.in %>%
    filter(get(icol_name) > 0 & !is.na(get(icol_name))) %>%
    complete(time = 2:max(time))%>%
    mutate_all(.funs = function(xx){ifelse(is.na(xx), 0, xx)}) %>%
    arrange(time)
  
  ts <- idat$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te<- ts+(window-1)
  
  wallinga_teunis(
    incid = pull(idat, eval(icol_name)),
    method = "parametric_si",
    config = list(
      mean_si = GI_mean,
      std_si = sqrt(GI_var),
      t_start=ts,
      t_end=te,
      n_sim=100
    )
  ) -> outs
  
  outs$R %>%
    mutate(time = if(wend == TRUE) t_end else ceiling((t_end+t_start)/2) ) %>%
    select(time, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
    setNames(c('time', paste0(outcol_name, '.mean'), paste0(outcol_name, '.025'), paste0(outcol_name, '.975')))
}



# df.in = rtdf
# icol_name = 'incidence'
# out_name = 'Cori'
# window = 1
# GI_mean = 8
# GI_var = 2*(parlist$true_mean_GI/2)^2
## Output cori estimate with mean, CI and times given an input df, and the name of the incidence column
get_cori <- function(df.in, 
                     icol_name, 
                     out_name = 'Cori',
                     window = 1, 
                     GI_mean=parlist$true_mean_GI, 
                     GI_var=2*(parlist$true_mean_GI/2)^2,
                     wend = TRUE){
  
  max.obs.time <- df.in %>% filter(!is.na(!!sym(icol_name))) %>% pull(time) %>% tail(1)
  
  
  idat <- df.in %>%
    #filter(get(icol_name) > 0 & !is.na(get(icol_name))) %>%
    complete(time = 2:max.obs.time) %>%
    arrange(time) %>%
    filter(time <= max.obs.time)
  idat[icol_name] <- na_to_0(idat[icol_name])
    #mutate(cleaned = ifelse(is.na(!!sym(icol_name)) & time <= max.obs.time, 0, !!sym(icol_name)))

  
  ts <- idat$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te <- ts+(window-1)
  
  estimate_R(
    incid = pull(idat, !!icol_name),
    method = "uncertain_si",
    config = make_config(
      list(
        mean_si = GI_mean,
        min_mean_si = GI_mean -1,
        max_mean_si = GI_mean + 1,
        std_mean_si = 1.5,
        std_std_si = 1.5,
        std_si = sqrt(GI_var),
        min_std_si = sqrt(GI_var)*.8,
        max_std_si = sqrt(GI_var)*1.2,
        n1 = 50,
        n2 = 100, 
        t_start=ts,
        t_end=te
      )
    )
  ) -> outs
  
  outs$R %>%
    mutate(time = if(wend == TRUE) t_end else ceiling((t_end+t_start)/2) ) %>%
    select(time, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
    setNames(c('time', paste0(out_name, '.mean'), paste0(out_name, '.025'), paste0(out_name, '.975')))
}


## Write a function to parse the BR outputs for plotting
parse_fits <- function(fits, max.time = parlist$n_t){
  
  if(typeof(fits) == 'list'){
    params <- fits
  }else{
    params <- rstan::extract(fits)
  }
  
  Rt <- params$Rt
  theta <- params$theta
  
  get_quantile <- function(mat, qq){
    apply(mat, 2, FUN = function(mm) quantile(mm, qq))
  }
  
  
  tibble(
    time = 1:max.time,
    BR.025 = get_quantile(Rt, 0.025),
    BR.mean = colMeans(Rt),
    BR.975 = get_quantile(Rt, 0.975)
  )
}


## Function to estimate BR 
get_BR <- function(df.in, filename, parlist, reset){
  ## Set up for stan fits
  source('../code/rtlive.R')
  onset_frac <- round(df.in$incidence %>% ifelse(is.na(.), 0, .) * 0.5)
  onset_draw <- df.in$incidence %>% ifelse(is.na(.), 0, .) %>% sapply(FUN = function(xx){rbinom(n = 1, size = xx, prob = .5)})
  onset <- df.in$incidence %>% ifelse(is.na(.), 0, .)
  cumulative_p_delay <- rep(1, length(onset_frac)) ## All SEIR infections are observed at the moment they occur, with no delay
  theta_init_mean <- 0.1
  theta_init_sd <- 0.1
  step_size_prior_sd <- 0.03
  
  ## These are the paramters for GI ~ gamma(2, 1/4)
  get_shape_rate <- function(mean, vv){
    rate = mean/vv
    shape = mean*rate
    stopifnot(shape/rate == mean)
    stopifnot(shape/rate^2 == vv)
    c(shape, rate)
  }
  ### Here, the si variance is the variance OF THE MEAN, so assume it's low.
  si_shape <- get_shape_rate(parlist$true_mean_GI, .5)[1]
  si_rate <- get_shape_rate(parlist$true_mean_GI, .5)[2]
  
  
  ## Estimate Rt using Rstan
  ## Compare against a version of the script that uses a smaller fudge factor (1e-6 instead of 0.1)
  if(!file.exists(filename)|reset){
    fit_best <- rtlive_stan(
      mod = '../code/rtlive.stan',
      onset,
      cumulative_p_delay,
      step_size_prior_sd,
      theta_init_mean, theta_init_sd,
      si_shape, si_rate,
      chains = 1
    )
    saveRDS(fit_best, filename)
  }else{
    fit_best <- readRDS(filename)
  }
  fit_best
}


## Wrapper to save a png using ggsave, without having to specify units and dpi
gg_png <- function(ww, ## width (in)
                   hh, ## height (in)
                   fn, ## filename
                   pp = last_plot()){ ## name of plot to save. default is the last plot in the working device
  ggsave(filename = fn, width =  ww, height = hh, plot = pp, units = 'in', dpi = 300, device = png())
}





# ## Get the first time at which cumulative case count is >12, or start at time 2 if cumulative incidence is immediately >12
# st.time <- max(
#   3,
#   df.in %>% 
#     mutate(cs = cumsum(get(icol_name))) %>%
#     filter(cs > 12) %>%
#     slice(1) %>%
#     pull(time) 
#   )
# 
# ## Get the last time at which incidence is >0 
# ed.time <- min(
#   max(df.in$time),
#   df.in %>%
#     filter(get(icol_name)==0 & time >= st.time) %>%
#     slice(1) %>%
#     pull(time)
#   
# )
# 
# idat <- df.in %>%
#   filter(time >= st.time & time <= ed.time)