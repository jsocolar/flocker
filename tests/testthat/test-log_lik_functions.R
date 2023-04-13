test_that("check log_lik functions work correctly", {
  # check dims
  ll_test <- log_lik_flocker(example_flocker_model_single)
  expect_equal(dim(ll_test), c(80, 900))
  ll_test <- log_lik_flocker(example_flocker_model_single, draw_ids = 1:20)
  expect_equal(dim(ll_test), c(20, 900))
  
  # check classes/values
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_equal(class(as.vector(ll_test)), "numeric")
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
})
