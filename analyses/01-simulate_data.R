## Simulate synthetic data
rm(list = ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


## Set parameters to simulate SEIR epidemic
parlist <- {
  list(
    N = 2e6, #total population size
    E_init = 0,
    I_init = 60,
    t_E = 4, # mean time in E (latent period)
    t_I = 4, # mean time in I (duration of infectiousness)
    n_t = 300, # total timesteps
    pre_intervention_R0 = 2.0, # Initial R0 before interventions
    intervention_R0 = 0.8, # Final R0 after interventions
    partially_lifeted_R0 = 1.15,
    intervention_time_1 = 60, # Timepoint at which intervention starts (at which underlying transmission rate begins to fall)
    intervention_time_2 = 60+30,
    days_intervention_to_min = c(7), # Days from intervention start until transmission rate hits min_R0
    days_to_Rt_rise = 7,
    model_types = c('seir'), # Can also choose sir
    methods = c('ode', 'stochastic') # could also choose ode
  )
}

## Derive the mean and variance of the serial interval from the input parameters
parlist$true_mean_SI = (parlist$t_E+parlist$t_I)
parlist$true_var_SI = 2*parlist$true_mean_SI^2



## Simulate SEIR data using a stochastic (ode) model
## Results are saved to a subdirectory called 'R0-xx.xx/'
source('../code/simulation-sweep.R')
sim_sweep(parlist)
testplots(parlist)


write_rds(parlist, path = 'true_pars.rds')
