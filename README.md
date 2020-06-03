# Code for analyses and figures shown in:

Practical considerations for detectingchanges in the effective reproductive number, $R_t$
Code by Katie Gostic, Ed Baskerville and Lauren McGough
\\
Last updated 3-June-2020

## [Code]()
This directory contains functions and wrappers used to perform analyses:

* `simulation.R` - code to generate synthetic data using an SIR or SEIR-type model, deterministic or stochastic.
* `funs_simulation-sweep.R` - wrapper functions for epidemic simulation.
* `infer_times_of_infection_observation.R` - functions to infer times of observation from SEIR times of infection, and to infer times of infection from times of observation by (1) drawing samples from a known delay distribution, or (2) shifting back in time by the mean delay to observation.
* `rtlive.R` and `rtlive.stan` together provide code to reporduce an [adaptation]() of the [Bettencourt & Ribeiro]() method for $R_t$ estimation popularized by [rt.live](https://rt.live).
* `util.R` - various utility functions, including wrappers to estimate $R_t$ using the methods of [Cori et al.](), [Wallinga & Teunis](), and using methods adapted from [Bettencourt & Ribeiro]() by [rt.live](https://rt.live). The first two methods are implemented in the package [EpiEstim](https://CRAN.R-project.org/package=EpiEstim). The final method uses the rstan implementation above.

## [Analyses]()
This directory contains scripts and notebooks used to run analyses and generate figures:

### Workflow:

* `01-simulate_data.R` - Specify inputs, generate synthetic data and save to a directory called `R0-xx/`.
* `02-Analyze_synthetic_data.Rmd` - Annotated code performs the analyses and makes all plots shown in Fig. 1, 4 and 5. 
* `03-Make_Fig_3.Rmd` - Annotated code performs the analyses and makes all plots shown in Fig. 3.
* `04-make_FigS1_deconvolution.Rmd` - Annotated code performs the analyses and makes the plots shown in Fig. S1, a tutorial on why deconvolution methods are needed.

