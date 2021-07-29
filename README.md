# flocker: flexible ocʞupancy estimation in R

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
Install `flocker` directly from GitHub with 
```
# install.packages("remotes")
remotes::install_github("jsocolar/flocker")
```
`flocker` requires a working installation of cmdstan, an interface to the 
probabilistic programming language Stan. We recommend installing cmdstan from 
R with
```
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```
Occasionally users encounter difficulties in installing cmdstan. To 
troubleshoot, consult https://mc-stan.org/cmdstanr/articles/cmdstanr.html. 
For more about Stan and cmdstan, see https://mc-stan.org/. If you really 
can't figure out installation, feel free to ask for help at 
https://discourse.mc-stan.org/.

### Getting started
This vignette shows how to fit several classes of model in flocker.  For 
additional classes of model, consult the relevant `brms` vignettes.

**add vignette links to a flocker vignette and to relevant brms vignettes**

### Citing flocker
Include citation details and suggestion to also cite `brms` and Stan.
