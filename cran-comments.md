# Cran comments

* Updated title to title case
  
* The goal of this package is to reduce the compilation time, of
  'rxode2' as requested by Prof Brian Ripley. The first attempt only
  split of the likelihoods ('rxode2ll').  This did not quite make the
  required install/check time of 10 minutes so the additional pieces
  were split off in this package.
  
* The first part that was split off in this package is the parsing
  (from dparser) and large linear compartment compiles (with
  derivatives) from 'stan' to allow nonlinear mixed effects models use
  solved linear compartmental models in 'nlmixr2'

* These stan functions will also remove the interactions between Eigen and 
  Armadillo that 'rxode2' was working around by putting all the Eigen/stan
  pieces outside of the 'rxode2' core

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
