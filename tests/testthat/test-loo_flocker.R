test_that("loo_flocker_onefit works correctly", {
  # read test-case flocker fit
  test_fit <- readRDS("../test_data/test_fit.rds")
  
  test_loo <- loo_flocker(test_fit)
  test_loo_thinned <- loo_flocker(test_fit, 5)
  test_loo_list <- loo_flocker(list(test_fit, test_fit))
  
  # check error
  expect_error(loo_flocker(), "argument \"x\" is missing, with no default")
  expect_error(loo_flocker(1), "x must be a flocker_fit object or a list of flocker_fit objects")
  expect_error(loo_flocker(list(1)), "x is a list, but x\\[\\[1\\]\\] is not a flocker_fit object.")
  expect_error(loo_flocker(list(test_fit, 1)), "x is a list, but x\\[\\[2\\]\\] is not a flocker_fit object.")
  
  # check dims
  expect_equal(attributes(test_loo)$dims, c(1000, 50))
  expect_equal(attributes(test_loo_thinned)$dims, c(200, 50))
  
  # check list output
  expect_identical(class(test_compare), "list")
  expect_identical(test_compare[[1]], test_compare[[2]])
  expect_identical(test_loo_list[[1]], test_loo)
})