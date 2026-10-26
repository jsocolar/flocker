test_that("check log_lik functions work correctly", {
  # single (package data model - 100 units, 4 draws)
  ll_test <- log_lik_flocker(fixture_fit("single"))
  expect_equal(dim(ll_test), c(4, 160))
  ll_test_sub <- log_lik_flocker(fixture_fit("single"), draw_ids = 1:2)
  expect_equal(dim(ll_test_sub), c(2, 160))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_equal(class(as.vector(ll_test)), "numeric")
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # single_C (100 units, 8 draws)
  ll_test <- log_lik_flocker(fixture_fit("single_C"))
  expect_equal(dim(ll_test), c(16, 100))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  testthat::skip_on_cran()
  
  # augmented (variable observed species plus 10 pseudospecies, 8 draws)
  ll_test <- log_lik_flocker(fixture_fit("aug"))
  n_group <- fixture_fit("aug")$data$ff_n_group[1]
  known_present <- fixture_fit("aug")$data$ff_group_Q[
    get_flocker_metadata(fixture_fit("aug"))$group_start
  ]
  n_observed <- sum(known_present)
  n_aug <- n_group - n_observed
  expect_equal(dim(ll_test)[1], 16)
  expect_equal(n_aug, 10)
  expect_equal(
    colnames(ll_test),
    c(
      paste0("observed_species", seq_len(n_observed)),
      paste0("pseudospecies", seq_len(n_aug))
    )
  )
  
  expect_equal(dim(ll_test)[2], n_group)
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # multiseason colex explicit (20 series, 8 draws)
  ll_test <- log_lik_flocker(fixture_fit("multi_colex_ex"))
  expect_equal(dim(ll_test), c(16, 20))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # multiseason colex equilibrium (20 series, 8 draws)
  ll_test <- log_lik_flocker(fixture_fit("multi_colex_eq"))
  expect_equal(dim(ll_test), c(16, 20))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # multiseason autologistic explicit (20 series, 8 draws)
  ll_test <- log_lik_flocker(fixture_fit("multi_auto_ex"))
  expect_equal(dim(ll_test), c(16, 20))
  expect_equal(class(ll_test), c("matrix", "array"))
  expect_lte(max(ll_test), 0)
  expect_false(any(is.infinite(ll_test)))
  expect_false(any(is.na(ll_test)))
  
  # multiseason autologistic equilibrium (20 series, 8 draws)
  ll_test <- log_lik_flocker(fixture_fit("multi_auto_eq"))
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
    log_lik_flocker(fixture_fit("single"), new_data = data.frame(x = 1)),
    "must be a flocker_data object"
  )
  
  # wrong data type should error
  expect_error(
    log_lik_flocker(fixture_fit("single"), new_data = fixture_data("multi_colex_ex")),
    "new_data is formatted for a different model type"
  )
  
  # --- equality checks: passing training data as new_data should match default ---
  
  # single (uses current-schema data from setup.R)
  ll_default <- log_lik_flocker(fixture_fit("single"))
  ll_newdata <- log_lik_flocker(
    fixture_fit("single"),
    new_data = fixture_data("single"),
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # single_C (uses fixture_data("single_C") from setup.R)
  ll_default <- log_lik_flocker(fixture_fit("single_C"))
  ll_newdata <- log_lik_flocker(
    fixture_fit("single_C"),
    new_data = fixture_data("single_C"),
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  testthat::skip_on_cran()

  # generic two-level single-season
  ll_default <- log_lik_flocker(fixture_fit("twolevel"))
  ll_newdata <- log_lik_flocker(
    fixture_fit("twolevel"),
    new_data = fixture_data("twolevel"),
    allow_new_levels = FALSE
  )
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)

  # augmented
  ll_default <- log_lik_flocker(fixture_fit("aug"))
  ll_newdata <- log_lik_flocker(
    fixture_fit("aug"),
    new_data = fixture_data("aug"),
    allow_new_levels = FALSE
  )
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # multiseason colex explicit (uses fixture_data("multi_colex_ex") from setup.R)
  ll_default <- log_lik_flocker(fixture_fit("multi_colex_ex"))
  ll_newdata <- log_lik_flocker(
    fixture_fit("multi_colex_ex"),
    new_data = fixture_data("multi_colex_ex"),
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # multiseason colex equilibrium (uses fixture_data("multi_colex_eq") from setup.R)
  ll_default <- log_lik_flocker(fixture_fit("multi_colex_eq"))
  ll_newdata <- log_lik_flocker(
    fixture_fit("multi_colex_eq"),
    new_data = fixture_data("multi_colex_eq"),
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # multiseason autologistic explicit (uses fixture_data("multi_auto_ex") from setup.R)
  ll_default <- log_lik_flocker(fixture_fit("multi_auto_ex"))
  ll_newdata <- log_lik_flocker(
    fixture_fit("multi_auto_ex"),
    new_data = fixture_data("multi_auto_ex"),
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
  
  # multiseason autologistic equilibrium (uses fixture_data("multi_auto_eq") from setup.R)
  ll_default <- log_lik_flocker(fixture_fit("multi_auto_eq"))
  ll_newdata <- log_lik_flocker(
    fixture_fit("multi_auto_eq"),
    new_data = fixture_data("multi_auto_eq"),
    allow_new_levels = FALSE
  )
  expect_equal(dim(ll_newdata), dim(ll_default))
  expect_equal(ll_newdata, ll_default, tolerance = 1e-10)
})

test_that("two-level log likelihood remains finite for very small probabilities", {
  n_unit <- 679L
  occ_lp <- matrix(1000, nrow = n_unit, ncol = 1)
  det_lp <- array(0, dim = c(n_unit, 2, 1))
  obs <- matrix(0, nrow = n_unit, ncol = 2)

  group_ll <- log_lik_twolevel_single_from_components(
    occ_lp, det_lp, matrix(1000, nrow = 1),
    rep(1L, n_unit), 1L, obs
  )
  expect_equal(group_ll[1, 1], n_unit * 2 * log(0.5))

  n_visit <- 1500L
  unit_obs <- matrix(0, nrow = 1, ncol = n_visit)
  unit_obs[1, 1] <- 1
  unit_obs[1, n_visit] <- NA
  unit_det_lp <- array(0, dim = c(1, n_visit, 1))
  unit_det_lp[1, n_visit, 1] <- NA
  unit_ll <- log_lik_twolevel_single_from_components(
    matrix(1000, nrow = 1),
    unit_det_lp,
    matrix(1000, nrow = 1), 1L, 1L,
    unit_obs
  )
  expect_equal(unit_ll[1, 1], (n_visit - 1) * log(0.5))
})
