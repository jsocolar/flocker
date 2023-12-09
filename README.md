<img src="man/figures/flocker_sticker.png" width = 200 alt="flocker logo" align = "right">

# flocker: flexible occupancy estimation in R
<!-- badges: start -->
[![R-CMD-check](https://github.com/jsocolar/flocker/workflows/R-CMD-check/badge.svg)](https://github.com/jsocolar/flocker/actions?workflow=R-CMD-check)
[![Coverage
Status](https://codecov.io/gh/jsocolar/flocker/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jsocolar/flocker)
<!-- badges: end -->


The `flocker` R package enables users to fit flexible occupancy models using 
the extended `lme4` formula syntax provided by R package `brms`. Also 
included are functions to format data for occupancy modeling and functions 
for model post-processing, including posterior predictions, the posterior
distribution of the latent true occupancy state, and model comparison. 
`flocker` is under active development: development priorities include better 
integrated functionality for posterior predictive checking and full
cross-validation.

Compared to existing R packages for fitting occupancy models, `flocker` is 
substantially more flexible in the variety of models that can 
be fit, and contains advanced functionality for model comparison and 
checking.

### Getting started
To get started, check out our 
[tutorial vignette, available here](https://jsocolar.github.io/flocker/articles/flocker_tutorial.html) and our
[introductory paper, available here](https://www.biorxiv.org/content/10.1101/2023.10.26.564080v1). 
For installation instructions, see below.

### Installation
Install the latest CRAN release of `flocker` with
```
install.packages("flocker")
```

Install the current development version of `flocker` with
```
# install.packages("remotes")
remotes::install_github("jsocolar/flocker")
```
`flocker` requires a working version of either `rstan` or `cmdstan`, which are 
interfaces to [Stan](https://mc-stan.org/), a state-of-the-art the 
probabilistic programming language. We recommend using `cmdstan` if possible.
To do so, first install R package `cmdstanr` with
```
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```
You must additionally install `cmdstan` itself. We strongly recommend using
`cmdstanr` to manage the cmdstan installation:
```
cmdstanr::install_cmdstan()
```
Both `rstan` and `cmdstan` require a working C++ toolchain, which has posed occasional 
complications for Stan users. 
[See here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for 
`rstan` focused advice on installing the toolchain, and 
[see here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) for `cmdstan` focused 
advice. If you encounter toolchain issues that you are unable to troubleshoot, 
ask for help at https://discourse.mc-stan.org/.

### Bugs, features, contributions
To request a feature or report a bug (much appreciated!) please 
[open an issue at the GitHub repository](https://github.com/jsocolar/flocker/issues).
To contribute to `flocker` (very much appreciated!) have a look at existing open 
issues, or open a new issue to discuss.

### Citing flocker
When using `flocker`, please cite the [companion manuscript](https://www.biorxiv.org/content/10.1101/2023.10.26.564080v1):
* Socolar, J.B. & Mills, S.C. (2023). "Introducing flocker: an R package for flexible occupancy modeling via brms and Stan." https://doi.org/10.1101/2023.10.26.564080

Additionally, please consider citing [Stan](https://mc-stan.org/users/citations/)
and R package [brms](https://mc-stan.org/users/interfaces/brms).

<img src="man/figures/logo2.png" width = 200 alt="AI generated" align = "">