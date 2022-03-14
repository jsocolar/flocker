test_that("loo_flocker_onefit works correctly", {
  # read test-case flocker fit
  test_fit <- readRDS("../test_data/test_fit.rds")
  
  # check error
  expect_error(loo_compare_flocker(test_fit), "model_list must be a list of flocker_fit objects.")
  expect_error(loo_compare_flocker(list(test_fit)), "model_list must contain at least two flocker_fit objects.")

  # check model naming
  test_compare <- loo_compare_flocker(list(test_fit, test_fit), model_names = c("m1", "m2"))
  expect_identical(row.names(test_compare), c("m1", "m2"))
  
  # check test output identical
  expect_true(all(apply(test_compare, 2, diff) == 0))
})