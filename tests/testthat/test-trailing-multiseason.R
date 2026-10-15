test_that("multiseason post-processing retains globally trailing seasons", {
  testthat::skip_on_cran()

  fit <- example_flocker_model_multi_colex_ex
  data_base <- format_ragged_multi_fixture()
  data_trailing <- format_ragged_multi_fixture(global_trailing_season = TRUE)
  fixture_trailing <- make_ragged_multi_fixture(global_trailing_season = TRUE)
  draw_ids <- 1:2

  fitted_base <- fitted_flocker(
    fit,
    components = "col",
    new_data = data_base,
    unit_level = TRUE,
    draw_ids = draw_ids
  )$linpred_col
  fitted_trailing <- fitted_flocker(
    fit,
    components = "col",
    new_data = data_trailing,
    unit_level = TRUE,
    draw_ids = draw_ids
  )$linpred_col
  expect_equal(fitted_trailing[, 1:5, ], fitted_base)
  expect_equal(
    is.na(fitted_trailing),
    array(
      rep(is.na(fixture_trailing$expected_unit), length(draw_ids)),
      dim = c(3, 6, length(draw_ids))
    ),
    check.attributes = FALSE
  )

  fitted_det <- fitted_flocker(
    fit,
    components = "det",
    new_data = data_trailing,
    unit_level = FALSE,
    draw_ids = draw_ids
  )$linpred_det
  expect_equal(
    is.na(fitted_det),
    array(
      rep(is.na(fixture_trailing$obs), length(draw_ids)),
      dim = c(3, 4, 6, length(draw_ids))
    ),
    check.attributes = FALSE
  )

  for(history_condition in c(FALSE, TRUE)) {
    Z_base <- get_Z(
      fit,
      draw_ids = draw_ids,
      history_condition = history_condition,
      new_data = data_base
    )
    Z_trailing <- get_Z(
      fit,
      draw_ids = draw_ids,
      history_condition = history_condition,
      new_data = data_trailing
    )
    expect_equal(
      as.vector(Z_trailing[, 1:5, ]),
      as.vector(Z_base)
    )
    expect_identical(
      is.na(Z_trailing),
      array(
        rep(is.na(fixture_trailing$expected_unit), length(draw_ids)),
        dim = c(3, 6, length(draw_ids))
      )
    )
  }

  expected_prediction_nas <- array(
    rep(is.na(fixture_trailing$obs), length(draw_ids)),
    dim = c(3, 4, 6, length(draw_ids))
  )
  for(history_condition in c(FALSE, TRUE)) {
    predictions <- predict_flocker(
      fit,
      draw_ids = draw_ids,
      history_condition = history_condition,
      new_data = data_trailing
    )
    expect_identical(is.na(predictions), expected_prediction_nas)
  }

  expect_equal(
    log_lik_flocker(fit, draw_ids = draw_ids, new_data = data_trailing),
    log_lik_flocker(fit, draw_ids = draw_ids, new_data = data_base),
    tolerance = 1e-10
  )
})
