# flocker: FLexible OCupancy Estimation in R

This package uses `brms` to generate Stan code for the linear predictors in
an single-season occupancy models (including multi-species models), then fits 
the model using an appropriate likelihood term. Also included are functions to 
format data for occupancy modeling.

Compared to existing packages for fitting occupancy models, `flocker`
distinguishes itself by virtue of its flexibility. Via `brms`, `flocker` 
supports all manner of random terms, including 
* random effects with arbitrary covariance structures spanning both the 
occupancy and detection sub-models
* phylogenetic models
* generalized additive models
* autoregressive models (ARMA)
* spatial autoregressive models (CAR)
* measurement error in the predictors
* and more!

When the model does not involve visit-specific detection covariates, it is fit 
natively in `brms` as a zero-inflated binomial regression. A `brmsfit` object is
returned that is compatible with additional packages in the Stan ecosystem just 
like any other `brmsfit`. When the model does involve visit-specific covariates, 
`brms` is used only to generate Stan code for the linear predictor, which is then
combined with the appropriate likelihood term to yield the desired model. The 
model is fit directly via `cmdstanr` or `rstan`.
