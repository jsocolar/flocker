# Regenerate from the repository root in a fresh R session:
# Rscript data-raw/example_flocker_model.R
# Load this checkout so the saved data and fit use its schema and Stan source.
pkgload::load_all(helpers = FALSE, quiet = TRUE)

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
  # Retain only four draws in total to keep the source package small.
  iter = 1000,
  warmup = 998,
  seed = 1,
  backend = "rstan",
  # Ship portable draws, not the platform-specific compiled shared library.
  stan_model_args = list(save_dso = FALSE),
  save_warmup = FALSE
)

save(mfd_single, file = "data/mfd_single.rda", compress = "xz")

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  example_flocker_model_single,
  file = "inst/extdata/example_flocker_model_single.rds",
  compress = "xz"
)
