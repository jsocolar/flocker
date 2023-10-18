## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----library, echo = FALSE----------------------------------------------------
library(flocker)
library(brms)

## ----simulate-----------------------------------------------------------------
d <- simulate_flocker_data()

## ----format-------------------------------------------------------------------
fd_rep_varying <- make_flocker_data(
  obs = d$obs, 
  unit_covs = d$unit_covs, 
  event_covs = d$event_covs
  )

## ----rep-varying--------------------------------------------------------------
rep_varying <- flock(
  f_occ = ~ uc1 + (1 + uc1 | species),
  f_det = ~ uc1 + ec1 + (1 + uc1 + ec1 | species),
  flocker_data = fd_rep_varying,
  cores = 4,
  silent = 2,
  refresh = 0
  )

## ----rep-constant-------------------------------------------------------------
fd_rep_constant <- make_flocker_data(
  obs = d$obs, 
  unit_covs = d$unit_covs
  )
rep_constant <- flock(
  f_occ = ~ uc1 + (1 + uc1 | species),
  f_det = ~ uc1 + (1 + uc1 | species),
  flocker_data = fd_rep_constant,
  save_pars = save_pars(all = TRUE), # for loo with moment matching
  silent = 2,
  refresh = 0,
  cores = 4
  )

## ----rep-constant-parallel, eval = FALSE--------------------------------------
#  rep_constant <- flock(
#    f_occ = ~ uc1 + (1 + uc1 | species),
#    f_det = ~ uc1 + (1 + uc1 | species),
#    flocker_data = fd_rep_constant,
#    silent = 2,
#    refresh = 0,
#    chains = 2,
#    cores = 2,
#    threads = 2
#    )

## ----multiseason-data---------------------------------------------------------
multi_data <- simulate_flocker_data(
  n_season = 3,
  n_pt = 300,
  n_sp = 1,
  multiseason = "colex", 
  multi_init = "explicit",
  seed = 1
  )

fd_multi <- make_flocker_data(
  multi_data$obs, 
  multi_data$unit_covs, 
  multi_data$event_covs, 
  type = "multi",
  quiet = TRUE
  )

## ----multi-colex--------------------------------------------------------------
multi_colex <- flock(
  f_occ = ~ uc1,
  f_det = ~ uc1 + ec1,
  f_col = ~ uc1,
  f_ex = ~ uc1,
  flocker_data = fd_multi,
  multiseason = "colex",
  multi_init = "explicit",
  cores = 4,
  silent = 2,
  refresh = 0
)

## ----multi-colex-eq-----------------------------------------------------------
multi_colex_eq <- flock(
  f_det = ~ uc1 + ec1,
  f_col = ~ uc1,
  f_ex = ~ uc1,
  flocker_data = fd_multi,
  multiseason = "colex",
  multi_init = "equilibrium",
  cores = 4,
  silent = 2,
  refresh = 0
  )

## ----multi-auto---------------------------------------------------------------
multi_auto <- flock(
  f_occ = ~ uc1,
  f_det = ~ uc1 + ec1,
  f_col = ~ uc1,
  f_auto = ~ 1,
  flocker_data = fd_multi,
  multiseason = "autologistic",
  multi_init = "explicit",
  cores = 4,
  silent = 2,
  refresh = 0
  )

## ----multi-auto-eq------------------------------------------------------------
multi_auto_eq <- flock(
  f_det = ~ uc1 + ec1,
  f_col = ~ uc1,
  f_auto = ~ 1,
  flocker_data = fd_multi,
  multiseason = "autologistic",
  multi_init = "equilibrium",
  cores = 4,
  silent = 2,
  refresh = 0
)

## ----data-augmented-----------------------------------------------------------
augmented_data <- simulate_flocker_data(
    augmented = TRUE
    )
fd_augmented <- make_flocker_data(
  augmented_data$obs, augmented_data$unit_covs, augmented_data$event_covs,
  type = "augmented", n_aug = 100,
  quiet = TRUE
  )
augmented <- flock(
  f_occ = ~ (1 | ff_species),
  f_det = ~ uc1 + ec1 + (1 + uc1 + ec1 | ff_species),
  augmented = TRUE,
  flocker_data = fd_augmented,
  cores = 4,
  silent = 2,
  refresh = 0
  )

## ----rep-constant-processing, eval = FALSE------------------------------------
#  predictions_rep_constant <- brms::posterior_predict(rep_constant)
#  loo_rep_constant <- brms::loo(rep_constant, moment_match = TRUE)
#  brms::conditional_effects(rep_constant)

## ----fitted-flocker, eval = FALSE---------------------------------------------
#  fitted_flocker(rep_constant)
#  fitted_flocker(rep_varying)
#  fitted_flocker(multi_colex)
#  fitted_flocker(augmented)

## ----get_Z, eval=FALSE--------------------------------------------------------
#  get_Z(rep_varying)

## ----predict_flocker, eval=FALSE----------------------------------------------
#  predict_flocker(rep_varying)

## -----------------------------------------------------------------------------
loo_compare_flocker(
  list(rep_constant, rep_varying)
)

## ----multi-season-loo---------------------------------------------------------
loo_compare_flocker(
  list(multi_colex, multi_colex_eq, multi_auto, multi_auto_eq)
)

## ----prior-summary------------------------------------------------------------
get_flocker_prior(
  f_occ = ~ uc1 + (1 + uc1 | species),
  f_det = ~ uc1 + ec1 + (1 + uc1 + ec1 | species),
  flocker_data = fd_rep_varying
  )
brms::prior_summary(rep_varying)

## ----priors_2, eval=FALSE-----------------------------------------------------
#  get_flocker_prior(
#    f_occ = ~ uc1 + (1|species),
#    f_det = ~ uc1 + ec2 + (1|species),
#    flocker_data = fd_rep_varying
#    )

## ----priors_3, eval=FALSE-----------------------------------------------------
#  user_prior <- c(brms::set_prior("normal(0, 1)", coef = "uc1"),
#                  brms::set_prior("normal(0, 3)", coef = "uc1", dpar = "occ"))

## ----priors1------------------------------------------------------------------
rep_varying_prior1 <- flock(
  f_occ = ~ uc1,
  f_det = ~ ec1,
  flocker_data = fd_rep_varying,
  prior = 
    brms::set_prior("logistic(0,1)", class = "Intercept") +
    brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "occ") +
    brms::set_prior("normal(0,2)", class = "b") +
    brms::set_prior("normal(0,2)", class = "b", dpar = "occ"),
  cores = 4,
  silent = 2,
  refresh = 0
  )
brms::prior_summary(rep_varying_prior1)

## ----priors2------------------------------------------------------------------
rep_varying_prior2 <- flock(
  f_occ = ~ 0 + Intercept + uc1,
  f_det = ~ 0 + Intercept + ec1,
  flocker_data = fd_rep_varying,
  prior = 
    brms::set_prior("normal(0,2)", class = "b") +
    brms::set_prior("normal(0,2)", class = "b", dpar = "occ") +
    brms::set_prior("normal(1, 1)", class = "b", coef = "Intercept") +
    brms::set_prior("normal(-1, 1)", class = "b", coef = "Intercept", dpar = "occ"),
  cores = 4,
  silent = 2,
  refresh = 0
  )
brms::prior_summary(rep_varying_prior2)

## ----lm, eval=FALSE-----------------------------------------------------------
#  mod1 <- flock(
#    f_occ = ~ uc1 + (1|species),
#    f_det = ~ 1,
#    flocker_data = fd_rep_constant
#    )

## ----lme4, eval=FALSE---------------------------------------------------------
#  mod2 <- flock(
#    f_occ = ~ uc1 + (1|species),
#    f_det = ~ 1,
#    flocker_data = fd_rep_constant
#    )

## ----lme4_2, eval=FALSE-------------------------------------------------------
#  mod3 <- flock(
#    f_occ = ~ uc1 + (1 + uc1 | species),
#    f_det = ~ ec1 + (1 | species),
#    flocker_data = fd_rep_varying
#    )

## ----lme4_3, eval=FALSE-------------------------------------------------------
#  mod4 <- flock(
#    f_occ = ~ uc1 + (1 |g1| species) + (0 + uc1 |g2| species),
#    f_det = ~ ec1 + (1 |g1| species) + (0 + ec1 |g2| species),
#    flocker_data = fd_rep_varying
#    )

## ----spline-gp, eval=FALSE----------------------------------------------------
#  mod5 <- flock(
#    f_occ = ~ s(uc1),
#    f_det = ~ t2(uc1, ec1),
#    flocker_data = fd_rep_varying
#    )
#  
#  mod6 <- flock(
#    f_occ = ~ 1,
#    f_det = ~ gp(uc1, ec1),
#    flocker_data = fd_rep_varying
#    )

## ----phylo, eval=FALSE--------------------------------------------------------
#  # simulate an example phylogeny
#  phylogeny <- ape::rtree(30, tip.label = paste0("sp_", 1:30))
#  
#  # calculate covariance matrix
#  A <- ape::vcv.phylo(phylogeny)
#  
#  mod8 <- flock(
#    f_occ = ~ 1 + (1|gr(species, cov = A)),
#    f_det = ~  1 + ec1 + (1|species),
#    flocker_data = fd_rep_varying,
#    data2 = list(A = A)
#    )
#  
#  mod9 <- flock(
#    f_occ = ~ 1 + (1|gr(species, cov = A)),
#    f_det = ~  1 + ec1 + (1|gr(species, cov = A)),
#    flocker_data = fd_rep_varying,
#    data2 = list(A = A)
#    )

