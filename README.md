# Code for analyses and figures shown in:

Practical considerations for detectingchanges in the effective reproductive number, Rt. 

Code by Katie Gostic, Ed Baskerville and Lauren McGough. 

Last updated 18-June-2020. 

## [Code](https://github.com/cobeylab/Rt_estimation/tree/master/code)
This directory contains functions and wrappers used to perform analyses:

* `simulation.R` - code to generate synthetic data using an SIR or SEIR-type model, deterministic or stochastic.
* `funs_simulation-sweep.R` - wrapper functions for epidemic simulation.
* `infer_times_of_infection_observation.R` - functions to infer times of observation from SEIR times of infection, and to infer times of infection from times of observation by (1) drawing samples from a known delay distribution, or (2) shifting back in time by the mean delay to observation.
* `rtlive.R` and `rtlive.stan` together provide code to reporduce an [adaptation](https://github.com/k-sys/covid-19/blob/master/Realtime%20Rt%20mcmc.ipynb) of the [Bettencourt & Ribeiro](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0002185) method for Rt estimation popularized by [rt.live](https://rt.live).
* `util.R` - various utility functions, including wrappers to estimate Rt using the methods of [Cori et al.](https://academic.oup.com/aje/article/178/9/1505/89262), [Wallinga & Teunis](https://academic.oup.com/aje/article/160/6/509/79472), and using methods adapted from [Bettencourt & Ribeiro](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0002185) by [rt.live](https://rt.live). The first two methods are implemented in the package [EpiEstim](https://CRAN.R-project.org/package=EpiEstim). The final method uses the rstan implementation above.
* `caseR.R` - Functions to calculate the exact case reproductive number within the synthetic data (dahsed black lines shown in Fig. 2 and Fig. B.2).
* `Rc_math.Rmd` - Notes on the math used to calculate the case reproductive number exactly.

## [Analyses](https://github.com/cobeylab/Rt_estimation/tree/master/analyses)
This directory contains scripts and notebooks used to run analyses and generate figures:

### Workflow:

* `01-simulate_data.R` - Specify inputs, generate synthetic data and save to a directory called `R0-xx/`.
* `02-Analyze_synthetic_data.Rmd` - Annotated code performs the analyses and makes all plots shown in Fig. 2, 5 and 6, and appendix Fig. B.3. 
* `02-make_Fig_2_appendix.Rmd` - Code performs analyses and generates Fig. B2, an alternate version of Fig. 2 in the main text, which compares the performance of Cori, Wallinga and Teunis, and Bettencourt and Ribeiro when the time series is truncated as Rt rises or falls.
* `03-Make_Fig_4.Rmd` - Annotated code performs the analyses and makes all plots shown in Fig. 4.
* `04-Make_Fig_deconvolution.Rmd` - Annotated code performs the analyses and makes the plot shown in the appendix, a tutorial on why deconvolution methods are needed.

