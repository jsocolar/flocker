test_that("loo_flocker_onefit works correctly", {
  suppressWarnings(test_loo <- loo_flocker(example_flocker_model_single))
  suppressWarnings(test_loo_thinned <- loo_flocker(example_flocker_model_single, 2))
  suppressWarnings(test_loo_list <- loo_flocker(rep(list(example_flocker_model_single), 2)))
  
  # check error
  expect_error(loo_flocker(), "argument \"x\" is missing, with no default")
  expect_error(loo_flocker(1), "x must be a flocker_fit object or a list of flocker_fit objects")
  expect_error(loo_flocker(list(1)), "x is a list, but x\\[\\[1\\]\\] is not a flocker_fit object.")
  expect_error(loo_flocker(list(example_flocker_model_single, 1)), "x is a list, but x\\[\\[2\\]\\] is not a flocker_fit object.")
  
  # check dims
  expect_equal(attributes(test_loo)$dims, c(8, 100))
  expect_equal(attributes(test_loo_thinned)$dims, c(4, 100))
  
  # check list output
  expect_identical(class(test_loo_list), "list")
  expect_identical(test_loo_list[[1]], test_loo_list[[2]])
  expect_identical(test_loo_list[[1]], test_loo)
})


test_that("loo_flocker_onefit works correctly", {
  # check error
  expect_error(loo_flocker_onefit(1), "x must be a flocker_fit object")
  
  suppressWarnings(test_loo <- loo_flocker_onefit(example_flocker_model_single, thin = NULL))
  suppressWarnings(test_loo_alt <- loo_flocker_onefit(example_flocker_model_single, thin=1))
  suppressWarnings(test_loo_thinned <- loo_flocker_onefit(example_flocker_model_single, thin = 2))
  
  # check dimensions
  expect_equal(attributes(test_loo)$dims, c(8, 100))
  expect_equal(attributes(test_loo_thinned)$dims, c(4, 100))
  
  # check thin = NULL specification returns same output 
  expect_identical(test_loo, test_loo_alt)
})

test_that("loo_flocker_onefit works correctly", {
  # check error
  expect_error(loo_compare_flocker(list(example_flocker_model_single)), "model_list must contain at least two flocker_fit objects.")
  
  # check model naming
  suppressWarnings(test_compare <- loo_compare_flocker(list(example_flocker_model_single, example_flocker_model_single), 
                                      model_names = c("m1", "m2"), thin = 2))
  expect_identical(row.names(test_compare), c("m1", "m2"))
  
  # check test output identical
  expect_true(all(apply(test_compare, 2, diff) == 0))
})
