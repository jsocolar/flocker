# Independent reference: a generously bounded direct sum. The bound is only a
# test oracle; the implemented likelihood has no abundance truncation.
n_mixture_poisson_binomial_direct <- function(y, eta, log_lambda) {
  lambda <- exp(log_lambda)
  N <- seq.int(max(y), max(max(y) + 300, ceiling(lambda + 40 * sqrt(lambda) + 300)))
  terms <- dpois(N, lambda, log = TRUE)
  for (t in seq_along(y)) terms <- terms + dbinom(y[t], N, plogis(eta[t]), log = TRUE)
  matrixStats::logSumExp(terms)
}

test_that("falling-factorial coefficients satisfy independent polynomial identities", {
  expect_equal(exp(n_mixture_binomial_log_coefficients(c(2, 3))), c(6, 6, 1))
  expect_equal(exp(n_mixture_binomial_log_coefficients(c(2, 2))), c(2, 4, 1))
  expect_identical(n_mixture_binomial_log_coefficients(c(0, 0)), 0)
  expect_identical(n_mixture_binomial_log_coefficients(5), 0)
  for (y in list(c(3, 5, 2, 0), c(2, 2, 2), c(0, 1), c(4, 0, 4))) {
    lc <- n_mixture_binomial_log_coefficients(y)
    for (k in 0:12) {
      degree <- seq_along(lc) - 1
      actual <- sum(exp(lc) * factorial(degree) * choose(k, degree))
      other <- y[-which.max(y)]
      expected <- prod(vapply(other, function(n) prod(max(y) + k - seq_len(n) + 1), numeric(1)))
      expect_equal(actual, expected, tolerance = 1e-10)
    }
    expect_equal(n_mixture_binomial_log_coefficients(rev(y)), lc)
  }
  # Ordinary coefficients overflow here; the stored logs remain finite.
  expect_true(all(is.finite(n_mixture_binomial_log_coefficients(c(200, 200, 200)))))
  for (bad in list(numeric(), NA_real_, -1, 0.5, Inf, "1")) {
    expect_error(n_mixture_binomial_log_coefficients(bad), "Counts")
  }
})

test_that("auxiliary data reach a fixed Stan likelihood independently of visit rows", {
  build <- function(counts) {
    auxiliary <- n_mixture_binomial_auxiliary_data(counts)
    data <- data.frame(y = unlist(counts, use.names = FALSE))
    # Hold priors fixed: brms otherwise embeds a response-dependent default
    # intercept prior, independently of our auxiliary-data declarations.
    args <- list(formula = brms::bf(y ~ 1, abund ~ 1), data = data,
      family = n_mixture_poisson_binomial(), stanvars = n_mixture_poisson_binomial_stanvars(auxiliary),
      prior = c(brms::set_prior("normal(0, 2)", class = "Intercept"),
                brms::set_prior("normal(0, 2)", class = "Intercept", dpar = "abund")))
    list(auxiliary = auxiliary, code = do.call(brms::make_stancode, args),
         data = do.call(brms::make_standata, args))
  }
  a <- build(list(c(2, 3), 0, c(0, 1, 2)))
  b <- build(list(0))
  expect_identical(a$code, b$code)
  expect_true(rstan::stanc(model_code = a$code)$status)
  expect_identical(a$auxiliary$ff_unit_offsets, c(1L, 3L, 4L, 7L))
  expect_identical(a$auxiliary$ff_coef_offsets, c(1L, 4L, 5L, 7L))
  for (name in names(a$auxiliary)) expect_equal(as.vector(a$data[[name]]), a$auxiliary[[name]])
  expect_equal(a$auxiliary$ff_coef_power, c(0, 1, 2, 0, 0, 1))
  expect_identical(dim(b$data$ff_coef_power), 1L)
  expect_equal(as.vector(b$data$ff_coef_power), 0)
  expect_identical(dim(b$data$ff_log_coef), 1L)
  expect_identical(dim(b$data$ff_max_count), 1L)
  expect_match(a$code, "n_mixture_poisson_binomial_lpmf(Y | mu, abund, ff_unit_offsets, ff_coef_offsets, ff_max_count, ff_log_factorial, ff_log_coef, ff_coef_power)", fixed = TRUE)
  expect_equal(nrow(a$data$X), 6)
  expect_error(n_mixture_binomial_auxiliary_data(list(numeric())), "Counts")
})

test_that("compiled exact likelihood and gradients agree with independent calculations", {
  skip_on_cran()
  skip_on_os("windows")
  source <- paste(readLines(system.file("stan", "n_mixture_poisson_binomial.stan",
    package = "flocker", mustWork = TRUE)), collapse = "\n")
  code <- paste0("functions {", source, "}
    data { int<lower=1> T; array[T] int y;
           int<lower=1> ff_nunit; int<lower=1> ff_ncoef;
           array[ff_nunit+1] int ff_unit_offsets;
           array[ff_nunit+1] int ff_coef_offsets;
           array[ff_nunit] int ff_max_count;
           vector[ff_nunit] ff_log_factorial; vector[ff_ncoef] ff_log_coef; vector[ff_ncoef] ff_coef_power; }
    parameters { real log_lambda; vector[T] eta; }
    model { target += n_mixture_poisson_binomial_lpmf(y | eta, rep_vector(log_lambda, T),
      ff_unit_offsets, ff_coef_offsets, ff_max_count, ff_log_factorial, ff_log_coef, ff_coef_power); }")
  model <- rstan::stan_model(model_code = code)
  cases <- list(
    list(y = c(3, 5, 2, 0), eta = c(-1, 0.5, 1, -2), log_lambda = log(8)),
    list(y = c(2, 2), eta = c(0.3, -0.8), log_lambda = log(3)),
    list(y = c(0, 0), eta = c(-3, 1), log_lambda = log(4)),
    list(y = 4, eta = -0.5, log_lambda = log(7)),
    list(y = c(0, 1), eta = c(-20, 15), log_lambda = -8),
    list(y = c(100, 90), eta = c(1, 0.7), log_lambda = log(140)))
  evaluate <- function(y, eta, log_lambda) {
    counts <- if (is.list(y)) y else list(y)
    auxiliary <- n_mixture_binomial_auxiliary_data(counts)
    for (name in setdiff(names(auxiliary), c("ff_nunit", "ff_ncoef"))) {
      auxiliary[[name]] <- as.array(auxiliary[[name]])
    }
    y <- unlist(counts, use.names = FALSE)
    fit <- rstan::sampling(model, data = c(list(T = length(y), y = as.array(as.integer(y))), auxiliary),
      algorithm = "Fixed_param", chains = 1, iter = 1, warmup = 0,
      init = list(list(log_lambda = log_lambda, eta = as.array(eta))), refresh = 0)
    pars <- c(log_lambda, eta)
    list(lp = rstan::log_prob(fit, pars), grad = rstan::grad_log_prob(fit, pars))
  }
  for (case in cases) {
    actual <- do.call(evaluate, case)
    expect_equal(actual$lp, do.call(n_mixture_poisson_binomial_direct, case), tolerance = 1e-8)
    pars <- c(case$log_lambda, case$eta)
    reference <- function(x) n_mixture_poisson_binomial_direct(case$y, x[-1], x[1])
    gradient <- vapply(seq_along(pars), function(i) {
      up <- down <- pars
      up[i] <- up[i] + 1e-5
      down[i] <- down[i] - 1e-5
      (reference(up) - reference(down)) / 2e-5
    }, numeric(1))
    expect_equal(as.numeric(actual$grad), gradient, tolerance = 1e-5)
  }
  # These predictor pairs have lambda*p approximately 1, although lambda
  # alone overflows. A one-visit model is exactly a thinned Poisson.
  for (y in c(0, 3)) {
    extreme <- evaluate(y, -1000, 1000)
    expect_equal(extreme$lp, dpois(y, 1, log = TRUE), tolerance = 1e-10)
    expect_true(all(is.finite(extreme$grad)))
  }
  expect_equal(evaluate(3, 1000, log(4))$lp, dpois(3, 4, log = TRUE), tolerance = 1e-10)
  counts <- list(c(2, 3), 0, c(0, 1, 2))
  eta <- c(-1, 0.3, 1, -0.5, 0, 0.7)
  expected <- n_mixture_poisson_binomial_direct(counts[[1]], eta[1:2], log(5)) +
    n_mixture_poisson_binomial_direct(counts[[2]], eta[3], log(5)) +
    n_mixture_poisson_binomial_direct(counts[[3]], eta[4:6], log(5))
  expect_equal(evaluate(counts, eta, log(5))$lp, expected, tolerance = 1e-10)
})
