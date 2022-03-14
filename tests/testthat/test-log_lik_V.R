test_that("check log_lik_V works correctly", {
  # read test-case flocker fit
  test_fit <- readRDS("../test_data/test_fit.rds")

  # check dims
  ll_test <- log_lik_V(test_fit)
  expect_equal(dim(ll_test), c(1000, 50))
  ll_test <- log_lik_V(test_fit, draw_ids = 1:50)
  expect_equal(dim(ll_test), c(50, 50))
  
  # check classes/values
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_equal(class(as.vector(ll_test)), "numeric")
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # check rep-constant error
  attributes(test_fit)$lik_type <- "C"
  expect_error(log_lik_V(test_fit), "flocker_fit_V works only for rep-varying flocker_fits")
})