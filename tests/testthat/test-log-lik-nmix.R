test_that("exact log likelihoods agree with independent abundance sums", {
  for (y in list(c(0, 0), 3, c(2, 4, 1), c(10, 10, 8))) {
    lambda <- c(.1, 5, 50)
    p <- seq(.1, .9, length.out = length(y))
    reference <- vapply(lambda, function(lam) {
      n <- 0:500
      terms <- dpois(n, lam, log = TRUE)
      for (j in seq_along(y)) terms <- terms + dbinom(y[j], n, p[j], log = TRUE)
      matrixStats::logSumExp(terms)
    }, numeric(1))
    expect_equal(nmix_poisson_binomial_log_lik(y, log(lambda),
      matrix(qlogis(p), length(y), 3)), reference, tolerance = 1e-10)
  }
  # Single-visit Poisson thinning, including detection near both boundaries.
  for (eta in c(-1000, -30, 0, 30, 1000)) {
    for (y in c(0, 3)) {
      expected <- if (y == 0) -exp(log(5) + log_inv_logit(eta)) else
        y * (log(5) + log_inv_logit(eta)) - exp(log(5) + log_inv_logit(eta)) - lgamma(y + 1)
      expect_equal(nmix_poisson_binomial_log_lik(y, log(5), matrix(eta)), expected)
    }
  }
})

test_that("likelihood reconstruction preserves histories, dimensions, and names", {
  d <- make_nmix_data(rbind(a = c(am = 2, pm = NA), b = c(am = 0, pm = 1)), quiet = TRUE)
  abund <- matrix(log(c(4, 7, 5, 8)), 2, 2,
                   dimnames = list(c("a", "b"), c("draw_3", "draw_1")))
  det <- array(qlogis(.3), c(2, 2, 2))
  det[1, 2, ] <- NA
  local_mocked_bindings(fitted_nmix = function(...) list(abund = abund, det = det))
  actual <- log_lik_nmix(NULL, new_data = d)
  expect_identical(dimnames(actual), rev(dimnames(abund)))
  expect_equal(unname(actual[, 1]), dpois(2, c(4, 5) * .3, log = TRUE))
  expect_equal(unname(actual[, 2]), nmix_poisson_binomial_log_lik(c(0, 1),
    log(c(7, 8)), matrix(qlogis(.3), 2, 2)))
})

test_that("LOO rejects invalid models and thinning", {
  for (thin in list(0, -1, 1.5, NA, Inf, c(1, 2), "2")) {
    expect_error(loo_nmix(NULL, thin), "thin must")
  }
  expect_error(loo_nmix(NULL), "nmix fit")
  expect_error(loo_nmix(list()), "nmix fit")
  expect_error(loo_compare_nmix(list(NULL, NULL)), "at least two nmix fits")
  expect_error(log_lik_nmix(NULL), "returned by nmix")
})

test_that("LOO thinning preserves chain blocks and uses likelihood-scale ESS", {
  fit <- structure(list(), class = c("brmsfit", "nmix_fit"))
  ll <- cbind(a = -1000 + sin(1:40), b = -1500 + cos(1:40))
  local_mocked_bindings(niterations = function(...) 20L,
                        nchains = function(...) 2L, .package = "brms")
  selected <- NULL
  local_mocked_bindings(log_lik_nmix = function(x, draw_ids) {
    selected <<- draw_ids
    ll[draw_ids, , drop = FALSE]
  })
  actual <- suppressWarnings(loo_nmix(fit, thin = 2))
  expect_equal(selected, c(seq(1, 20, 2), seq(21, 40, 2)))
  subset <- ll[selected, , drop = FALSE]
  reff <- loo::relative_eff(exp(sweep(subset, 2, apply(subset, 2, max), "-")),
                            chain_id = rep(1:2, each = 10))
  expected <- suppressWarnings(loo::loo(subset, r_eff = reff))
  expect_equal(actual$estimates, expected$estimates)
  expect_identical(rownames(actual$pointwise), c("a", "b"))
})
