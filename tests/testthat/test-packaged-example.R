test_that("packaged example data use the current schema", {
  env <- new.env()
  utils::data("mfd_single", package = "flocker", envir = env)
  data <- env$mfd_single
  expect_named(data, c("data", "metadata"))
  expect_identical(data$metadata$schema_version, 2L)
  expect_identical(data$metadata$data_role, "model_data")
  expect_silent(validate_packed_schema(data))
  expect_identical(dim(get_positions(data)), data$metadata$original_dim)
})

test_that("the compact packaged fit supports public post-processing", {
  # Avoid loading the serialized brmsfit and its optional namespaces on CRAN.
  skip_on_cran()
  skip_if_not_installed("rstan")
  fit <- readRDS(system.file("extdata", "example_flocker_model_single.rds",
                            package = "flocker", mustWork = TRUE))
  env <- new.env()
  utils::data("mfd_single", package = "flocker", envir = env)
  # brms may reorder columns and add its explicit Intercept column.
  for (name in names(env$mfd_single$data)) {
    expect_equal(fit$data[[name]], env$mfd_single$data[[name]])
  }
  expect_identical(get_flocker_metadata(fit), env$mfd_single$metadata)
  expect_equal(brms::ndraws(fit), 4)
  dims <- c(env$mfd_single$metadata$original_dim, 4L)
  fitted <- fitted_flocker(fit)
  expect_identical(dim(fitted$linpred_det), dims)
  expect_identical(dim(predict_flocker(fit)), dims)
  expect_identical(dim(predict_flocker(fit, history_condition = TRUE)), dims)
  states <- get_Z(fit)
  expect_true(all(is.na(states) | (states >= 0 & states <= 1)))
  ll <- log_lik_flocker(fit)
  expect_true(all(is.finite(ll)))
})
