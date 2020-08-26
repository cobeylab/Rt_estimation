#!/usr/bin/env Rscript
source('../code/simulation.R')
## This script defines wrapper functions and utility functions to generate, load and parse synthetic data.
## The underlying simulation functions are defined in ../../simulation.R


## MAIN FUNCTION - 
## Simulate data for various intervention times, speed of change in Rt and initial R0 values
## Save outputs to a new directory
## Input PARAMS is a list with elements:

sim_sweep <- function(PARAMS, path_out= NULL){
  
  dirname <- sprintf("R0-%.1f", PARAMS$pre_intervention_R0)
  dirname <- paste0(path_out, dirname)
  if(!dir.exists(dirname)){
    sprintf('creating new output directory, %s', dirname)
    dir.create(dirname)
  }
  setwd(dirname)
  # saveRDS(PARAMS, paste0('PARAMS_df.Rds'))
  for(intervention_time_1 in PARAMS$intervention_time_1){
    for(decrease_duration in PARAMS$days_intervention_to_min){
      for(intervention_time_2 in PARAMS$intervention_time_2){
        for(increase_duration in PARAMS$days_to_Rt_rise){
          ## Specify time-varying R0
          if(PARAMS$pre_intervention_R0 == PARAMS$intervention_R0 & PARAMS$intervention_R0 == PARAMS$partially_lifeted_R0){
            arnaught <- rep(PARAMS$pre_intervention_R0, PARAMS$n_t+1)
          }else{
            arnaught <- specify_arnaught(R0_vec = c(PARAMS$pre_intervention_R0, PARAMS$intervention_R0, PARAMS$partially_lifeted_R0), 
                                         change_start_vec = c(intervention_time_1, intervention_time_2), 
                                         change_end_vec =  c(intervention_time_1+decrease_duration,intervention_time_2+increase_duration), 
                                         NT = PARAMS$n_t)
          }
          do_one_arnaught(arnaught, intervention_time_1, decrease_duration, PARAMS)
        }
      }
    }
  }
  setwd('..')
}

do_one_arnaught <- function(arnaught, intervention_time, decrease_duration, PARAMS) {
  for(method in PARAMS$methods) {
    for(model_type in PARAMS$model_types) {
      sim_df <- if(model_type == 'sir') {
        simulate_sir_example(
          arnaught = arnaught,
          t_I = PARAMS$t_I,
          N = PARAMS$N, I_init = PARAMS$I_init,
          n_t = PARAMS$n_t,
          method = method
        ) %>%
          mutate(
            incidence = round(dS),
            obs_cases = NA
          )
      } else {
        simulate_seir_example(
          arnaught = arnaught,
          t_E = PARAMS$t_E, t_I = PARAMS$t_I,
          N = PARAMS$N, E_init = PARAMS$E_init, I_init = PARAMS$I_init,
          n_t = PARAMS$n_t,
          method = method
        ) %>%
          mutate(
            incidence = round(dS),
            obs_cases = round(dEI)
          )
      }
      
      saveRDS(
        list(
          sim_df = sim_df %>% 
            mutate(true_r0 = arnaught,
                   true_rt = arnaught*S/(S+E+I+R)),
          arnaught = arnaught,
          intervention_time = intervention_time,
          decrease_duration = decrease_duration,
          
          params = PARAMS,
          method = method,
          model_type = model_type
        ),
        sprintf('%s_%s_dec%.0f-%.0f_sim.rds', model_type, method, intervention_time, decrease_duration)
      )
    }
  }
}
  
 
#sim_sweep(PARAMS)
#END UP WITH one folder for each R0, which contains: SEIR stoc {decreases}; SEIR ode {decreases}, parameters



## Check outputs
testplots <- function(PARAMS) {
  for(arnaught in PARAMS$pre_intervention_R0) {
    for(model_type in c('seir')) {
      for(method in PARAMS$methods) {
        # Load results from all interventions applied to a given R0, and bind into a single data frame
       sim_results <- data.frame(fns = list.files(path = sprintf('R0-%.1f', arnaught)), stringsAsFactors = FALSE) %>%
         filter(grepl(fns, pattern = method)) %>%
         filter(grepl(fns, pattern = '_sim')) %>%
         pull(fns) %>%
         lapply(FUN = function(xx){
           readRDS(paste0(sprintf('R0-%.1f/', arnaught), xx)) -> tmp
           tmp$sim_df %>% 
             mutate(int_time = as.character(tmp$intervention_time), 
                    dec_dur = as.character(tmp$decrease_duration))
         }) %>%
         bind_rows()
        
            cat(sprintf('Plotting %s, %s', model_type, method))
            # sim_result <- readRDS(sprintf('R0-%.1f/%s_%s_dec%.0f-%.0f_sim.rds', 
            #                               arnaught, model_type, method, intervention_time, decrease_duration))
            
            ## Plot prevalence in each compartment
            sim_results %>%
              select(time:R, int_time, dec_dur) %>%
              pivot_longer(S:R, names_to = 'Compartment', values_to = 'Prevalence') %>%
              ggplot() +
              geom_line(aes(x = time, y = Prevalence, color = int_time, lty = dec_dur)) +
              facet_wrap(.~Compartment, scales = 'free_y') +
              ggtitle(sprintf('%s - %s - R0=%.1f', model_type, method, arnaught))
            if(!dir.exists('simplots/')){ dir.create('simplots/')}
            ggsave(filename = sprintf('simplots/prevalence_R0-%.1f_%s_%s.png', arnaught, model_type, method), width = 11, height = 8.5)
            
            ## Plot incidence in each compartment
            sim_results %>%
              select(time, int_time, dec_dur, dEI, dIR, incidence) %>%
              pivot_longer(dEI:incidence, names_to = 'transition', values_to = 'incidence') %>%
              ggplot() +
              geom_line(aes(x = time, y = incidence, color = int_time, linetype = dec_dur), size = 1, alpha = .5) +
              theme_bw()+
              facet_wrap(.~transition) +
              ggtitle(sprintf('%s - %s - R0=%.1f', model_type, method, arnaught))
            ggsave(filename = sprintf('simplots/incidence_R0-%.1f_%s_%s.png', arnaught, model_type, method), width = 11, height = 8.5)
          }
        }
      }
    }

#testplots(PARAMS)




## ADDITIONAL UTILITY FUNCTIONS
## Output underlying R0 values for input into the SIR or SEIR model across once change (increase or decrease)
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
    n_t = n_t, n_steps_per_t = 1,
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
get_sim_df <- function(method = 'ode'){ # can also input 'sotchastic'
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

