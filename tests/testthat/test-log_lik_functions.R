test_that("check log_lik functions work correctly", {
  # single (package data model - 100 units, 4 draws)
  ll_test <- log_lik_flocker(example_flocker_model_single)
  expect_equal(dim(ll_test), c(4, 160))
  ll_test_sub <- log_lik_flocker(example_flocker_model_single, draw_ids = 1:2)
  expect_equal(dim(ll_test_sub), c(2, 160))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_equal(class(as.vector(ll_test)), "numeric")
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # single_C (100 units, 8 draws)
  ll_test <- log_lik_flocker(example_flocker_model_single_C)
  expect_equal(dim(ll_test), c(16, 100))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # augmented (20 species total, 8 draws)
  ll_test <- log_lik_flocker(example_flocker_model_aug)
  expect_equal(dim(ll_test), c(16, 20))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # multiseason colex explicit (20 series, 8 draws)
  ll_test <- log_lik_flocker(example_flocker_model_multi_colex_ex)
  expect_equal(dim(ll_test), c(16, 20))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # multiseason colex equilibrium (20 series, 8 draws)
  ll_test <- log_lik_flocker(example_flocker_model_multi_colex_eq)
  expect_equal(dim(ll_test), c(16, 20))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # multiseason autologistic explicit (20 series, 8 draws)
  ll_test <- log_lik_flocker(example_flocker_model_multi_auto_ex)
  expect_equal(dim(ll_test), c(16, 20))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # multiseason autologistic equilibrium (20 series, 8 draws)
  ll_test <- log_lik_flocker(example_flocker_model_multi_auto_eq)
  expect_equal(dim(ll_test), c(16, 20))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
})

test_that("log_lik_flocker new_data argument works correctly", {
  
  # --- input validation ---
  
  # non-flocker_data should error
  expect_error(
    log_lik_flocker(example_flocker_model_single, new_data = data.frame(x = 1)),
    "must be a flocker_data object"
  )
  
  # wrong data type should error
  expect_error(
    log_lik_flocker(example_flocker_model_single, new_data = mfd_multi_colex_ex),
    "data type in new_data does not match"
  )
  
  # --- equality checks: passing training data as new_data should match default ---
  
  # single (uses shipped package data mfd_single)
  ll_default <- log_lik_flocker(example_flocker_model_single)
  ll_newdata <- log_lik_flocker(
    example_flocker_model_single,
    new_data = mfd_single,
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # single_C (uses mfd_single_C from setup.R)
  ll_default <- log_lik_flocker(example_flocker_model_single_C)
  ll_newdata <- log_lik_flocker(
    example_flocker_model_single_C,
    new_data = mfd_single_C,
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # multiseason colex explicit (uses mfd_multi_colex_ex from setup.R)
  ll_default <- log_lik_flocker(example_flocker_model_multi_colex_ex)
  ll_newdata <- log_lik_flocker(
    example_flocker_model_multi_colex_ex,
    new_data = mfd_multi_colex_ex,
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # multiseason colex equilibrium (uses mfd_multi_colex_eq from setup.R)
  ll_default <- log_lik_flocker(example_flocker_model_multi_colex_eq)
  ll_newdata <- log_lik_flocker(
    example_flocker_model_multi_colex_eq,
    new_data = mfd_multi_colex_eq,
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # multiseason autologistic explicit (uses mfd_multi_auto_ex from setup.R)
  ll_default <- log_lik_flocker(example_flocker_model_multi_auto_ex)
  ll_newdata <- log_lik_flocker(
    example_flocker_model_multi_auto_ex,
    new_data = mfd_multi_auto_ex,
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # multiseason autologistic equilibrium (uses mfd_multi_auto_eq from setup.R)
  ll_default <- log_lik_flocker(example_flocker_model_multi_auto_eq)
  ll_newdata <- log_lik_flocker(
    example_flocker_model_multi_auto_eq,
    new_data = mfd_multi_auto_eq,
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
})