# flocker: flexible ocʞupancy estimation in R
<!-- badges: start -->
[![R-CMD-check](https://github.com/jsocolar/flocker/workflows/R-CMD-check/badge.svg)](https://github.com/jsocolar/flocker/actions?workflow=R-CMD-check)
[![Coverage
Status](https://codecov.io/gh/jsocolar/flocker/branch/main/graph/badge.svg)](https://codecov.io/gh/jsocolar/flocker)
<!-- badges: end -->


The `flocker` R package enables users to fit flexible occupancy models using 
the same extended `lme4` formula syntax employed by R package `brms`. Also 
included are functions to format data for occupancy modeling and to extract 
estimates and predictions from fitted occupancy models.

Compared to existing R packages for fitting occupancy models, `flocker` is 
substantially more flexible in the variety of single-season models that can 
be fit, including 
* random effects, optionally with sophisticated covariance structures spanning 
both the occupancy and detection sub-models
* phylogenetic models
* generalized additive models
* spatially autoregressive models
* measurement error in the predictors
* and more!

### Installation
Install `flocker` with 
```
# install.packages("remotes")
remotes::install_github("jsocolar/flocker")
```
`flocker` requires a working installation of cmdstan, which is an interface to
the probabilistic programming language Stan. We recommend installing cmdstan 
from R with
```
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```
Users occasionally report difficulties installing cmdstan. To 
troubleshoot, consult https://mc-stan.org/cmdstanr/articles/cmdstanr.html. 
For more about Stan and cmdstan, see https://mc-stan.org/. If you really 
can't figure out installation, feel free to ask for help at 
https://discourse.mc-stan.org/.

Finally, `flocker` currently requires the development version of `brms` (this
is a short-term situation until `brms` is updated on CRAN). Install via
```
remotes::install_github("paul-buerkner/brms")
```

### Getting started
Soon we'll have a vignette here and some links to `brms` vignettes.  In the mean
time, try running
```
library(flocker)
example_data <- example_flocker_data()
fd <- make_flocker_data(example_data$obs, example_data$site_covs, example_data$visit_covs)
ff <- flock(f_occ = ~ sc1 + sc2 + (1|grp),
              f_det = ~ sc1 + vc1 + vc2 + (1|grp),
              flocker_data = fd,
              refresh = 50, chains = 1, iter_warmup = 5, iter_sampling = 200,
              adapt_engaged = F, step_size = .05, max_treedepth = 5, seed = 123)
summary(ff)
```

### Citing flocker
Please cite `flocker` as:
* Socolar, J.B. & Mills, S.C. (2021). "flocker: flexible occupancy estimation in 
R." R package version XXX, <URL: https://github.com/jsocolar/flocker/>.

Additionally, please consider citing [Stan](https://mc-stan.org/users/citations/)
and R package [brms](https://mc-stan.org/users/interfaces/brms).
