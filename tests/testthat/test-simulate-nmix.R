test_that("N-mixture simulation preserves RNG state, covariates, and ragged visits", {
  set.seed(281)
  before <- .Random.seed
  sim <- simulate_nmix_data(20, 4, ragged_rep = TRUE, seed = 12)
  expect_identical(.Random.seed, before)
  expect_identical(sim, simulate_nmix_data(20, 4, ragged_rep = TRUE, seed = 12))
  expect_true(all(rowSums(!is.na(sim$obs)) >= 1))
  expect_true(anyNA(sim$obs))
  expect_identical(is.na(sim$obs), is.na(sim$event_covs$ec1))
  expect_true(all(sim$obs <= sim$abundance, na.rm = TRUE))
  expect_true(all(sim$obs == floor(sim$obs), na.rm = TRUE))
  expect_silent(make_nmix_data(sim$obs, sim$unit_covs, sim$event_covs, quiet = TRUE))
  covs <- list(uc1 = seq_len(5), ec1 = matrix(seq_len(15), 5, 3))
  explicit <- simulate_nmix_data(5, 3, covariates = covs)
  expect_equal(explicit$unit_covs$uc1, covs$uc1)
  expect_equal(unname(explicit$event_covs$ec1), covs$ec1)
  simulate_nmix_data(seed = NULL)
  expect_false(identical(.Random.seed, before))
  expect_identical(dim(simulate_nmix_data(1, 1, ragged_rep = TRUE)$obs), c(1L, 1L))
})

test_that("simulated counts have Poisson-binomial moments and shared abundance", {
  coefs <- c(abund_intercept = log(6), abund_slope_unit = 0,
             det_intercept = qlogis(.4), det_slope_unit = 0, det_slope_visit = 0)
  sim <- simulate_nmix_data(30000, 2, params = coefs, seed = 8)
  expect_lt(abs(mean(sim$abundance) - 6), .08)
  expect_lt(abs(mean(sim$obs) - 2.4), .05)
  expect_lt(abs(var(sim$obs[, 1]) - 2.4), .10)
  expect_lt(abs(cov(sim$obs[, 1], sim$obs[, 2]) - .96), .08)
})

test_that("simulation rejects malformed inputs", {
  expect_error(simulate_nmix_data(n_pt = 0), "n_pt")
  expect_error(simulate_nmix_data(n_rep = 1.5), "n_rep")
  expect_error(simulate_nmix_data(params = 1:5), "params")
  expect_error(simulate_nmix_data(ragged_rep = NA), "ragged_rep")
  expect_error(simulate_nmix_data(covariates = list(foo = 1)), "covariates")
  expect_error(simulate_nmix_data(covariates = list(uc1 = NA)), "uc1")
  expect_error(simulate_nmix_data(covariates = list(ec1 = 1)), "ec1")
})

test_that("uncentered abundance intercepts are supported for SBC", {
  sim <- simulate_nmix_data(10, 3)
  data <- make_nmix_data(sim$obs, sim$unit_covs, sim$event_covs, quiet = TRUE)
  sd <- nmix_standata(~ 0 + Intercept + uc1, ~ 0 + Intercept + uc1 + ec1, data)
  expect_identical(colnames(sd$X_abund), c("Intercept", "uc1"))
  expect_equal(unname(sd$X_abund[, 1]), rep(1, nrow(data$data)))
  expect_error(nmix_standata(~ 0 + Intercept + ec1, ~ 1, data), "unit covariates")
})
