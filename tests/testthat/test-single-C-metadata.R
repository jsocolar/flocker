test_that("rep-constant data use the common metadata contract", {
  obs <- rbind(c(1, 0), c(0, NA), c(0, 1))
  fd <- make_flocker_data(obs, data.frame(x = 1:3), quiet = TRUE)
  expect_named(fd, c("data", "metadata"))
  expect_identical(fd$metadata$type, "single_C")
  expect_identical(fd$metadata$schema_version, 2L)
  expect_identical(fd$metadata$data_role, "model_data")
  expect_identical(fd$metadata$unit_covs, "x")
  expect_equal(fd$data$ff_n_trial, c(2, 1, 2))
  expect_equal(fd$data$ff_n_suc, c(1, 0, 1))
  expect_no_error(flocker_standata(f_occ = ~ x, f_det = ~ 1, flocker_data = fd))
  expect_no_error(flocker_stancode(f_occ = ~ x, f_det = ~ 1,
                                  flocker_data = fd, threads = 2))
  nd <- make_flocker_data(obs, data.frame(x = 1:3), newdata_checks = TRUE, quiet = TRUE)
  expect_identical(nd$metadata$data_role, "new_data")
  expect_error(flocker_standata(f_occ = ~ x, f_det = ~ 1, flocker_data = nd),
               "new_data cannot")
  old <- structure(c(list(data = fd$data), fd$metadata), class = "flocker_data")
  expect_error(get_flocker_metadata(old), "Incompatible flocker_data structure")
  fd$metadata$schema_version <- NULL
  expect_error(flocker_standata(f_occ = ~ x, f_det = ~ 1, flocker_data = fd),
               "Incompatible rep-constant single-season")
})

test_that("rep-constant prediction uses nested metadata and rejects old fits", {
  testthat::skip_if(windows_cran)
  fit <- example_flocker_model_single_C
  expect_no_error(predict_flocker(fit, draw_ids = 1:2, mixed = TRUE,
                     allow_new_levels = TRUE, sample_new_levels = "gaussian"))
  attr(fit, "flocker_metadata") <- NULL
  expect_error(predict_flocker(fit, draw_ids = 1:2),
               "Incompatible rep-constant single-season")
})
