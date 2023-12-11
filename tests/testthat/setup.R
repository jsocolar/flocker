## code to prepare `example_flocker_model_xx` datasets

set.seed(1)

#### single-season rep-constant model ####
suppressWarnings({
  fd <- simulate_flocker_data(n_pt = 20, n_sp = 5)
  mfd <- make_flocker_data(
    fd$obs, 
    fd$unit_covs
  )
  example_flocker_model_single_C <- flock(
    f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
    f_det = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
    flocker_data = mfd,
    cores = 1,
    iter = 8,
    save_warmup = FALSE
  )

#### data augmented model ####
  fd <- simulate_flocker_data(augmented = TRUE,
                              n_sp = 10, n_pt = 20)
  mfd <- make_flocker_data(
    fd$obs, 
    fd$unit_covs, 
    fd$event_covs,
    type = "augmented",
    n_aug = 10
  )
  example_flocker_model_aug <- flock(
    f_occ = ~ 0 + Intercept + (1 | ff_species),
    f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | ff_species),
    flocker_data = mfd,
    augmented = TRUE,
    cores = 1,
    iter = 8,
    save_warmup = FALSE
  )

#### multiseason colex with explicit inits ####
  fd <- simulate_flocker_data(
    n_pt = 20, n_sp = 1, n_season = 6, 
    multiseason = "colex", multi_init = "explicit",
    ragged_rep = TRUE, missing_seasons = TRUE
  )
  mfd <- make_flocker_data(
    fd$obs, fd$unit_covs, fd$event_covs, type = "multi"
  )
  example_flocker_model_multi_colex_ex <- flock(
    f_occ = ~ 0 + Intercept + uc1,
    f_det = ~ 0 + Intercept + uc1 + ec1,
    f_col = ~ 0 + Intercept + uc1,
    f_ex = ~ 0 + Intercept + uc1,
    flocker_data = mfd,
    multiseason = "colex",
    multi_init = "explicit",
    cores = 1,
    iter = 8,
    save_warmup = FALSE
  )

#### multiseason colex with equilibrium inits ####
  fd <- simulate_flocker_data(
    n_pt = 20, n_sp = 1, n_season = 6, 
    multiseason = "colex", multi_init = "equilibrium",
    ragged_rep = TRUE, missing_seasons = TRUE
  )
  mfd <- make_flocker_data(
    fd$obs, fd$unit_covs, fd$event_covs, type = "multi"
  )
  example_flocker_model_multi_colex_eq <- flock(
    f_det = ~ 0 + Intercept + uc1 + ec1,
    f_col = ~ 0 + Intercept + uc1,
    f_ex = ~ 0 + Intercept + uc1,
    flocker_data = mfd,
    multiseason = "colex",
    multi_init = "equilibrium",
    cores = 1,
    iter = 8,
    save_warmup = FALSE
  )

#### multiseason autologistic with explicit inits ####
  fd <- simulate_flocker_data(
    n_pt = 20, n_sp = 1, n_season = 6, 
    multiseason = "autologistic", multi_init = "explicit",
    ragged_rep = TRUE, missing_seasons = TRUE
  )
  mfd <- make_flocker_data(
    fd$obs, fd$unit_covs, fd$event_covs, type = "multi"
  )
  example_flocker_model_multi_auto_ex <- flock(
    f_occ = ~ 0 + Intercept + uc1,
    f_det = ~ 0 + Intercept + uc1 + ec1,
    f_col = ~ 0 + Intercept + uc1,
    f_auto = ~ 0 + Intercept + uc1,
    flocker_data = mfd,
    multiseason = "autologistic",
    multi_init = "explicit",
    cores = 1,
    iter = 8,
    save_warmup = FALSE
  )

#### multiseason autologistic with equilibrium inits ####
  fd <- simulate_flocker_data(
    n_pt = 20, n_sp = 1, n_season = 6, 
    multiseason = "autologistic", multi_init = "equilibrium",
    ragged_rep = TRUE, missing_seasons = TRUE
  )
  mfd <- make_flocker_data(
    fd$obs, fd$unit_covs, fd$event_covs, type = "multi"
  )
  example_flocker_model_multi_auto_eq <- flock(
    f_det = ~ 0 + Intercept + uc1 + ec1,
    f_col = ~ 0 + Intercept + uc1,
    f_auto = ~ 0 + Intercept + uc1,
    flocker_data = mfd,
    multiseason = "autologistic",
    multi_init = "equilibrium",
    cores = 1,
    iter = 8,
    save_warmup = FALSE
  )
})