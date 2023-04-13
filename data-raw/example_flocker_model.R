## code to prepare `example_flocker_fit` dataset goes here
set.seed(1)
fd <- make_flocker_data(
  example_flocker_data$obs, 
  example_flocker_data$unit_covs, 
  example_flocker_data$event_covs
)
example_flocker_model_single <- flock(
  f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
  f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | species),
  flocker_data = fd,
  backend = "cmdstanr",
  cores = 4,
  iter = 500,
  warmup = 480
)

usethis::use_data(example_flocker_model_single, overwrite = TRUE)
