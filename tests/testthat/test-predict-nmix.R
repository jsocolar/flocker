test_that("unconditioned predictions share abundance across visits", {
  n <- 30000L
  lps <- list(abund = matrix(log(c(4, 10)), 2, n),
              det = array(qlogis(.5), c(2, 3, n)))
  lps$det[1, 3, ] <- NA
  calls <- 0L
  local_mocked_bindings(fitted_nmix = function(...) {
    calls <<- calls + 1L
    lps
  })
  set.seed(361)
  out <- predict_nmix(NULL)
  expect_identical(calls, 1L)
  expect_identical(dim(out), c(2L, 3L, n))
  expect_true(all(is.na(out[1, 3, ])))
  expect_true(all(out[!is.na(out)] >= 0 & out[!is.na(out)] == floor(out[!is.na(out)])))
  # Poisson thinning: E[Y_j] = lambda*p_j, Cov(Y_j,Y_k) = lambda*p_j*p_k.
  expect_lt(abs(mean(out[1, 1, ]) - 2), .06)
  expect_lt(abs(mean(out[2, 2, ]) - 5), .10)
  expect_lt(abs(stats::cov(out[1, 1, ], out[1, 2, ]) - 1), .10)
  expect_lt(abs(stats::cov(out[2, 1, ], out[2, 2, ]) - 2.5), .20)
})

test_that("conditioned predictions replicate rather than copy observed counts", {
  n <- 20000L
  data <- make_nmix_data(matrix(3, 1, 1), newdata_checks = TRUE, quiet = TRUE)
  lps <- list(abund = matrix(log(4), 1, n), det = array(qlogis(.5), c(1, 1, n)))
  local_mocked_bindings(fitted_nmix = function(...) lps)
  set.seed(983)
  out <- predict_nmix(NULL, new_data = data, history_condition = TRUE)
  # N | y = 3 + Poisson(2); a replicated count has mean 2.5.
  expect_lt(abs(mean(out) - 2.5), .06)
  expect_true(any(out < 3))
  expect_true(any(out > 3))
})

test_that("prediction flags and fits are validated", {
  expect_error(predict_nmix(NULL, history_condition = NA), "history_condition")
  expect_error(predict_nmix(NULL), "returned by nmix")
})

test_that("unit and draw order are retained when repeating abundance across visits", {
  abundance <- matrix(c(2, 5, 11, 17), 2, 2)
  det <- array(Inf, c(2, 3, 2),
               dimnames = list(c("oak", "ash"), c("am", "noon", "pm"),
                               c("draw_3", "draw_1")))
  det[1, 3, ] <- NA
  local_mocked_bindings(
    fitted_nmix = function(...) list(abund = abundance, det = det),
    nmix_abundance_from_predictors = function(...) abundance
  )
  out <- predict_nmix(NULL)
  expect_identical(dimnames(out), dimnames(det))
  expect_equal(unname(out[, 1, ]), abundance)
  expect_equal(unname(out[, 2, ]), abundance)
  expect_equal(unname(out[2, 3, ]), abundance[2, ])
  expect_true(all(is.na(out[1, 3, ])))
})
