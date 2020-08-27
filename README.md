# Code for analyses and figures shown in:

Practical considerations for detecting changes in the effective reproductive number, Rt. 

Last updated 27-Aug-2020. 



## [Code](https://github.com/cobeylab/Rt_estimation/tree/master/code)
This directory contains functions and wrappers used to perform analyses and generate figures:

* `simulation.R` - code to generate synthetic data using an SIR or SEIR-type model, deterministic or stochastic.
* `funs_simulation-sweep.R` - wrapper functions for epidemic simulation.
* `infer_times_of_infection_observation.R` - functions to infer times of observation from SEIR times of infection, and to infer times of infection from times of observation by (1) drawing samples from a known delay distribution, or (2) shifting back in time by the mean delay to observation.
* `rtlive.R` and `rtlive.stan` together provide code to reproduce an [adaptation](https://github.com/k-sys/covid-19/blob/master/Realtime%20Rt%20mcmc.ipynb) of the [Bettencourt & Ribeiro](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0002185) method for Rt estimation popularized by [rt.live](https://rt.live).
* `util.R` - various utility functions, including wrappers to estimate Rt using the methods of [Cori et al.](https://academic.oup.com/aje/article/178/9/1505/89262), [Wallinga & Teunis](https://academic.oup.com/aje/article/160/6/509/79472), and using methods adapted from [Bettencourt & Ribeiro](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0002185) by [rt.live](https://rt.live). The first two methods are implemented in the package [EpiEstim](https://CRAN.R-project.org/package=EpiEstim). The final method uses the rstan implementation above.
* `caseR.R` - Functions to calculate the exact case reproductive number within the synthetic data (dashed black lines shown in Fig. 2 and Fig. B.2).
* `Rc_math.Rmd` - Notes on the math used to calculate the case reproductive number exactly.

## [Analyses](https://github.com/cobeylab/Rt_estimation/tree/master/analyses)
This directory contains scripts and notebooks used to run analyses and generate figures:

### Workflow:

* `01-simulate_data.R` - Specify inputs, generate synthetic data and save to a directory called `R0-xx/`.
* `02-...` - Various notebooks estimate Rt from synthetic data and generate plots.
* `Run_all_scripts.R` - Runs the entire workflow. Comments within indicate which notebooks generate which figures.

