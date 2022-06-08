test_that("loo_flocker_onefit works correctly", {
  # read test-case flocker fit
  test_fit <- readRDS("../test_data/test_fit.rds")
  
  # check error
  expect_error(loo_flocker_onefit(1), "x must be a flocker_fit object")
  
  test_loo <- loo_flocker_onefit(test_fit, thin = NULL)
  test_loo_alt <- loo_flocker_onefit(test_fit, thin=1)
  test_loo_thinned <- loo_flocker_onefit(test_fit, thin = 5)
  
  # check dimensions
  expect_equal(attributes(test_loo)$dims, c(1000, 50))
  expect_equal(attributes(test_loo_thinned)$dims, c(200, 50))
  
  # check thin = NULL specification returns same output 
  expect_identical(test_loo, test_loo_alt)
})