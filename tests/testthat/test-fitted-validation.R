test_that("N-mixture predicates distinguish data and fit objects", {
  counts <- make_nmix_data(matrix(0, 2, 2), quiet = TRUE)
  occupancy <- make_flocker_data(matrix(0, 2, 2), quiet = TRUE)
  expect_true(is_nmix_data(counts))
  expect_false(is_nmix_data(occupancy))
  expect_false(is_nmix_data(list(metadata = list(type = "nmix"))))
  expect_false(is_nmix_data(NULL))
  expect_false(is_nmix_fit(counts))
  expect_false(is_nmix_fit(structure(list(), class = c("brmsfit", "flocker_fit"))))
  expect_false(is_nmix_fit(structure(list(), class = "nmix_fit")))
  expect_true(is_nmix_fit(structure(list(), class = c("brmsfit", "nmix_fit"))))
})

test_that("fitted-value checks validate flags, intervals, data, and draw indices", {
  # The public example becomes a placeholder during R CMD check. Use the
  # existing test fixture instead, including under check-like environments.
  fit <- example_flocker_model_single2
  check <- function(...) {
    args <- list(fit = fit, new_data = NULL, response = TRUE, summarise = FALSE,
                  CI = c(0.05, 0.95), draw_ids = NULL)
    args <- utils::modifyList(args, list(...), keep.null = TRUE)
    do.call(validate_fitted_params, args)
  }
  expect_identical(check(), seq_len(brms::ndraws(fit)))
  expect_identical(check(draw_ids = c(3L, 1L)), c(3L, 1L))
  expect_silent(check(response = FALSE, summarise = TRUE))
  expect_silent(check(response = TRUE, summarise = FALSE))
  for (bad in list(NA, NULL, logical(), c(TRUE, FALSE), 1)) {
    expect_error(check(response = bad), "response must be TRUE or FALSE")
    expect_error(check(summarise = bad), "summarise must be TRUE or FALSE")
  }
  for (bad in list(numeric(), NA, 0, -1, 1.5, Inf, brms::ndraws(fit) + 1)) {
    expect_error(check(draw_ids = bad), "draw_ids must")
  }
  expect_error(check(CI = "x"), "CI must be numeric")
  expect_error(check(CI = 0.5), "CI must be of length 2")
  expect_error(check(CI = c(NA, 1)), "between zero and one")
  expect_error(check(CI = c(-0.1, 1)), "between zero and one")
  expect_error(check(CI = c(0.9, 0.1)), "increasing")
  expect_error(check(new_data = fit$data), "formatted data")
  expect_silent(check(new_data = fit$data, allow_dataframe = TRUE))
  expect_error(check(new_data = make_nmix_data(matrix(0, 2, 2), quiet = TRUE)),
               "different model type")
  expect_error(fitted_flocker(fit, response = NA), "response must")
  expect_error(fitted_flocker(fit, summarise = NA), "summarise must")
  expect_error(fitted_flocker(fit, draw_ids = 0), "draw_ids must")
  out <- fitted_flocker(fit, components = "det", draw_ids = 1)
  expect_identical(tail(dim(out$linpred_det), 1), 1L)
})
