expect_singleton_draw_dimension <- function(x) {
  expect_identical(tail(dim(x), 1), 1L)
}

expect_singleton_fitted_output <- function(flocker_fit) {
  fitted <- fitted_flocker(flocker_fit, draw_ids = 1)
  lapply(fitted, function(x) {
    expect_singleton_draw_dimension(x)
    expect_identical(tail(dimnames(x), 1)[[1]], "draw_1")
  })
  invisible(fitted)
}

test_that("single-season post-processing supports one draw", {
  expect_singleton_fitted_output(example_flocker_model_single2)
  expect_singleton_draw_dimension(
    get_Z(example_flocker_model_single2, draw_ids = 1)
  )
  expect_singleton_draw_dimension(
    predict_flocker(example_flocker_model_single2, draw_ids = 1)
  )
  expect_identical(
    nrow(log_lik_flocker(example_flocker_model_single2, draw_ids = 1)),
    1L
  )
})

test_that("two-level post-processing supports one draw", {
  testthat::skip_on_cran()

  expect_singleton_fitted_output(example_flocker_model_twolevel)
  Z <- get_Z(example_flocker_model_twolevel, draw_ids = 1)
  expect_named(Z, c("unit", "level2"))
  expect_singleton_draw_dimension(Z$unit)
  expect_singleton_draw_dimension(Z$level2)
  expect_equal(rownames(Z$level2), letters[1:3])
  sampled_Z <- get_Z(
    example_flocker_model_twolevel, draw_ids = 1, sample = TRUE
  )
  packed_group <- get_unit_group(example_flocker_model_twolevel$data)
  unit_order <- get_flocker_metadata(example_flocker_model_twolevel)$unit_order
  original_group <- integer(length(packed_group))
  original_group[unit_order] <- packed_group
  expect_true(all(
    sampled_Z$unit <= sampled_Z$level2[original_group, , drop = FALSE]
  ))
  expect_singleton_draw_dimension(
    predict_flocker(example_flocker_model_twolevel, draw_ids = 1)
  )
  expect_identical(
    nrow(log_lik_flocker(example_flocker_model_twolevel, draw_ids = 1)),
    1L
  )
  expect_equal(
    colnames(log_lik_flocker(example_flocker_model_twolevel, draw_ids = 1)),
    letters[1:3]
  )
})

test_that("multiseason post-processing supports one draw", {
  testthat::skip_on_cran()

  expect_singleton_fitted_output(example_flocker_model_multi_colex_ex)
  expect_singleton_draw_dimension(
    get_Z(example_flocker_model_multi_colex_ex, draw_ids = 1)
  )
  expect_singleton_draw_dimension(
    predict_flocker(example_flocker_model_multi_colex_ex, draw_ids = 1)
  )
  expect_identical(
    nrow(log_lik_flocker(example_flocker_model_multi_colex_ex, draw_ids = 1)),
    1L
  )
})
