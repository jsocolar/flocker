# Independent finite-sum references use a generous bound only in tests.
nmix_abundance_reference <- function(y, lambda, p) {
  n <- 0:500
  lp <- stats::dpois(n, lambda, log = TRUE)
  for (j in seq_along(y)) lp <- lp + stats::dbinom(y[j], n, p[j], log = TRUE)
  w <- exp(lp - max(lp))
  w <- w / sum(w)
  list(n = n, probability = w, mean = sum(n * w))
}

test_that("exact conditional abundance means agree with independent sums", {
  for (y in list(c(0, 0), 3, c(2, 4, 1), c(10, 10, 8))) {
    for (lambda in c(0.1, 5, 50)) {
      p <- seq(0.1, 0.9, length.out = length(y))
      expected <- nmix_abundance_reference(y, lambda, p)
      actual <- nmix_poisson_binomial_abundance(
        y, log(lambda), matrix(qlogis(p), ncol = 1), FALSE
      )
      expect_equal(actual, expected$mean, tolerance = 1e-10)
      expect_gte(actual, max(y))
    }
  }
  # One visit: observed individuals plus an independent missed Poisson count.
  expect_equal(nmix_poisson_binomial_abundance(
    3, log(c(2, 7)), matrix(qlogis(c(.2, .8)), 1), FALSE), c(4.6, 4.4))
  expect_equal(nmix_poisson_binomial_abundance(
    c(0, 0), log(5), matrix(qlogis(c(.2, .3)), 2), FALSE), 2.8)
  expect_equal(nmix_poisson_binomial_abundance(
    c(2, 2), log(5), matrix(1000, 2), FALSE), 2)
})

test_that("conditional samples follow the shifted Poisson mixture", {
  set.seed(185)
  y <- c(2, 4, 1)
  p <- c(.2, .7, .4)
  expected <- nmix_abundance_reference(y, 8, p)
  n <- 20000L
  actual <- nmix_poisson_binomial_abundance(
    y, rep(log(8), n), matrix(qlogis(p), 3, n), TRUE
  )
  expect_true(all(actual >= max(y) & actual == floor(actual)))
  expect_equal(mean(actual), expected$mean, tolerance = .03)
  for (k in 4:10) {
    expect_lt(abs(mean(actual == k) - expected$probability[k + 1]), .015)
  }
})

test_that("get_N validates its flags and fitted model", {
  expect_error(get_N(NULL, history_condition = NA), "history_condition")
  expect_error(get_N(NULL, sample = c(TRUE, FALSE)), "sample")
  expect_error(get_N(NULL), "returned by nmix")
})
