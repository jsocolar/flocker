## code to prepare `example_flocker_model_single`
library(flocker)

set.seed(1)

#### single-season model ####
fd <- simulate_flocker_data(n_pt = 20, n_sp = 8)
mfd_single <- make_flocker_data(
  fd$obs, 
  fd$unit_covs,
  fd$event_covs
)
example_flocker_model_single <- flock(
  f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
  f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | species),
  flocker_data = mfd_single,
  cores = 2,
  chains = 2,
  iter = 1000,
  warmup = 998,
  save_warmup = FALSE
)

usethis::use_data(mfd_single, overwrite = TRUE)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  example_flocker_model_single,
  file = "inst/extdata/example_flocker_model_single.rds",
  compress = "xz"
)
