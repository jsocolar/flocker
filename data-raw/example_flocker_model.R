## code to prepare `example_flocker_model_xx` datasets

set.seed(1)

#### single-season model ####
mfd <- make_flocker_data(
  example_flocker_data$obs, 
  example_flocker_data$unit_covs, 
  example_flocker_data$event_covs
)
example_flocker_model_single <- flock(
  f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
  f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | species),
  flocker_data = mfd,
  backend = "cmdstanr",
  cores = 4,
  iter = 1000,
  warmup = 990
)

#### single-season rep-constant model ####
mfd <- make_flocker_data(
  example_flocker_data$obs, 
  example_flocker_data$unit_covs
)
example_flocker_model_single_C <- flock(
  f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
  f_det = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
  flocker_data = mfd,
  backend = "cmdstanr",
  cores = 4,
  iter = 1000,
  warmup = 990
)

#### data augmented model ####
fd <- simulate_flocker_data(augmented = TRUE,
                            n_sp = 50)
mfd <- make_flocker_data(
  fd$obs, 
  fd$unit_covs, 
  fd$event_covs,
  type = "augmented",
  n_aug = 100
)
example_flocker_model_aug <- flock(
  f_occ = ~ 0 + Intercept + (1 | ff_species),
  f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | ff_species),
  flocker_data = mfd,
  augmented = TRUE,
  backend = "cmdstanr",
  cores = 4,
  iter = 1000,
  warmup = 990
  )

#### multiseason colex with explicit inits ####
fd <- simulate_flocker_data(
  n_pt = 200, n_sp = 1, n_season = 6, 
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
  cores = 4,
  backend = "cmdstanr",
  iter = 1000,
  warmup = 990
)

#### multiseason colex with equilibrium inits ####
fd <- simulate_flocker_data(
  n_pt = 200, n_sp = 1, n_season = 6, 
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
  cores = 4,
  backend = "cmdstanr",
  iter = 1000,
  warmup = 990
)

#### multiseason autologistic with explicit inits ####
fd <- simulate_flocker_data(
  n_pt = 200, n_sp = 1, n_season = 6, 
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
  cores = 4,
  backend = "cmdstanr",
  iter = 1000,
  warmup = 990
)

#### multiseason autologistic with equilibrium inits ####
fd <- simulate_flocker_data(
  n_pt = 200, n_sp = 1, n_season = 6, 
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
  cores = 4,
  backend = "cmdstanr",
  iter = 1000,
  warmup = 990
)

#### single-season fp model ####
fd <- simulate_flocker_data(n_sp = 1, n_pt = 200, fp = .9)
mfd <- make_flocker_data(
  fd$obs, 
  fd$unit_covs,
  fd$event_covs,
  fp = TRUE
)
example_flocker_model_single_fp <- flock(
  f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
  f_det = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
  flocker_data = mfd,
  fp = TRUE,
  backend = "cmdstanr",
  cores = 4,
  iter = 1000,
  warmup = 990
)

#### multiseason colex with explicit inits and fp ####
fd <- simulate_flocker_data(
  n_pt = 200, n_sp = 1, n_season = 6, 
  multiseason = "colex", multi_init = "explicit",
  ragged_rep = TRUE, missing_seasons = TRUE,
  fp = .9
)
mfd <- make_flocker_data(
  fd$obs, fd$unit_covs, fd$event_covs, type = "multi", fp = TRUE
)
example_flocker_model_multi_colex_ex_fp <- flock(
  f_occ = ~ 0 + Intercept + uc1,
  f_det = ~ 0 + Intercept + uc1 + ec1,
  f_col = ~ 0 + Intercept + uc1,
  f_ex = ~ 0 + Intercept + uc1,
  flocker_data = mfd,
  multiseason = "colex",
  multi_init = "explicit",
  fp = TRUE,
  cores = 4,
  backend = "cmdstanr",
  iter = 1000,
  warmup = 990
)

usethis::use_data(example_flocker_model_single, overwrite = TRUE)
usethis::use_data(example_flocker_model_single_C, overwrite = TRUE)
usethis::use_data(example_flocker_model_single_fp, overwrite = TRUE)
usethis::use_data(example_flocker_model_aug, overwrite = TRUE)
usethis::use_data(example_flocker_model_multi_colex_ex, overwrite = TRUE)
usethis::use_data(example_flocker_model_multi_colex_eq, overwrite = TRUE)
usethis::use_data(example_flocker_model_multi_auto_ex, overwrite = TRUE)
usethis::use_data(example_flocker_model_multi_auto_eq, overwrite = TRUE)
usethis::use_data(example_flocker_model_multi_colex_ex_fp, overwrite = TRUE)
