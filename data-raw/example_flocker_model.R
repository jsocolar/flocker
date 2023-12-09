## code to prepare `example_flocker_model_xx` datasets

set.seed(1)

#### single-season model ####
fd <- simulate_flocker_data(n_pt = 20, n_sp = 5)
mfd <- make_flocker_data(
  fd$obs, 
  fd$unit_covs,
  fd$event_covs
)
example_flocker_model_single <- flock(
  f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
  f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | species),
  flocker_data = mfd,
  cores = 2,
  chains = 2,
  iter = 1000,
  warmup = 996,
  save_warmup = FALSE
)

usethis::use_data(example_flocker_model_single, overwrite = TRUE)
