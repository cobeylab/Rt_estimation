## Run all analyses

## Generate synthetic data
set.seed(32)
source('01-simulate_data.R')

## Fig 2
knitr::knit('02-make_fig_compare_estimators.Rmd')

## Fig 3 
knitr::knit('02-make_fig_GI.Rmd')

## Fig 4 is a diagram drawn in .ppt

## Fig 5
knitr::knit('02-Make_Fig_backcalc.Rmd')

## Fig 6
knitr::knit('02-make_fig_smoothing_window.Rmd')

## Fig S1 (deconvolution)
knitr::knit('02-Make_Fig_deconvolution.Rmd')

## Fig S2 (compare estimators appendix)
knitr::knit('02-make_compare_appendix.Rmd')

## Fig S3 (smoothed Cori vs WT)
knitr::knit('02-make_fig_Sup_smooth_WT_Cori.Rmd')

## Save session info
sessionInfo() %>%
  write_rds(path = 'Session_Info.rds')
