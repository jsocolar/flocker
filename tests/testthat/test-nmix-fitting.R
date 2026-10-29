nmix_fitting_data <- function() {
  obs <- rbind(c(2, 0, 1), c(0, NA, NA), c(3, 2, NA), c(1, 1, 0))
  make_nmix_data(obs, data.frame(x = c(-1, 0, 1, 2),
                                  g = factor(c("a", "a", "b", "b"))),
                list(w = matrix(seq_len(12) / 10, 4)), quiet = TRUE)
}

test_that("N-mixture formulas and auxiliary data reach the likelihood", {
  d <- nmix_fitting_data()
  sd <- nmix_standata(~ x, ~ x + w, d)
  expect_equal(as.vector(sd$Y), d$data$ff_y)
  expect_equal(unname(sd$X[, "w"]), d$data$w)
  expect_equal(unname(sd$X_abund[, "x"]), d$data$x)
  for (name in names(d$auxiliary)) {
    expect_equal(as.vector(sd[[name]]), d$auxiliary[[name]])
  }
  priors <- c(brms::set_prior("normal(0, 1)", class = "Intercept"),
              brms::set_prior("normal(2, 1)", class = "Intercept", dpar = "abund"))
  code <- nmix_stancode(~ x, ~ x + w, d, prior = priors)
  expect_match(code, "n_mixture_poisson_binomial_lpmf(Y | mu, abund,", fixed = TRUE)
  expect_match(code, "normal_lpdf(Intercept | 0, 1)", fixed = TRUE)
  expect_match(code, "normal_lpdf(Intercept_abund | 2, 1)", fixed = TRUE)
  expect_true(rstan::stanc(model_code = code)$status)

  # Formula-local transformations must remain resolvable.
  det <- local({ twice <- function(x) 2 * x; ~ twice(w) })
  expect_equal(unname(nmix_standata(~ x, det, d)$X[, 2]), 2 * d$data$w)
  covariance <- diag(2)
  dimnames(covariance) <- list(levels(d$data$g), levels(d$data$g))
  correlated <- nmix_standata(~ x + (1 | gr(g, cov = A)), ~ w,
                              d, data2 = list(A = covariance))
  expect_true(any(grepl("^Lcov", names(correlated))))
})

test_that("N-mixture fitting rejects incompatible data and formulas", {
  d <- nmix_fitting_data()
  expect_error(n_mix(~ 1, ~ 1, make_flocker_data(matrix(0, 2, 2), quiet = TRUE)),
               "requires N-mixture data")
  bad <- d
  bad$metadata$schema_version <- 1L
  expect_error(n_mix(~ 1, ~ 1, bad), "Incompatible")
  bad <- d
  bad$metadata$data_role <- "new_data"
  expect_error(n_mix(~ 1, ~ 1, bad), "new_data cannot")
  bad <- d
  bad$auxiliary <- NULL
  expect_error(n_mix(~ 1, ~ 1, bad), "auxiliary")
  bad <- d
  bad$auxiliary$ff_unit_offsets[2] <- 1L
  expect_error(n_mix(~ 1, ~ 1, bad), "offsets")
  expect_error(n_mix(~ w, ~ 1, d), "unit covariates")
  expect_error(n_mix(~ 1, ~ ff_Q, d), "reserved")
  expect_error(n_mix(~ 1, ~ ., d), "explicitly")
  expect_error(n_mix(~ 1, y ~ 1, d), "one-sided")
  expect_error(n_mix(~ 1, brms::bf(det ~ 1, abund ~ 1), d), "one-sided")
  expect_error(n_mix(~ 1, ~ 1, d, threads = 2), "Threading")
  expect_error(n_mix(~ ar(time = x), ~ 1, d), "Autoregressive")
  expect_error(n_mix(~ 1, ~ arma(time = x), d), "Autoregressive")
  expect_error(n_mix(~ 1, ~ 1, d, autocor = ~ ar(time = x)), "Autoregressive")
  expect_error(n_mix(~ 1, ~ 1, d, family = poisson()), "internally")
  expect_error(n_mix(~ 1, ~ 1, d, prior = "normal(0, 1)"), "brmsprior")
  expect_error(suppressWarnings(nmix_standata(~ 1, ~ log(x), d)),
               "packed visits|infinite|NA|finite")
})

test_that("Poisson-binomial N-mixture models fit through the public interface", {
  skip_on_cran()
  skip_on_os("windows")
  set.seed(414)
  abundance <- rpois(60, 5)
  obs <- matrix(rbinom(180, rep(abundance, 3), 0.6), 60, 3)
  d <- make_nmix_data(obs, quiet = TRUE)
  fit <- n_mix(~ 1, ~ 1, d, chains = 1, iter = 1200, warmup = 400,
                seed = 415, refresh = 0, cores = 1,
                prior = c(brms::set_prior("normal(0, 1.5)", class = "Intercept"),
                          brms::set_prior("normal(1, 2)", class = "Intercept",
                                          dpar = "abund")))
  expect_s3_class(fit, "brmsfit")
  expect_s3_class(fit, "n_mix_fit")
  expect_identical(get_flocker_metadata(fit), d$metadata)
  expect_identical(attr(fit, "flocker_metadata"), d$metadata)
  expect_equal(fit$data$ff_y, d$data$ff_y)
  draws <- as.data.frame(fit)
  expect_true(all(is.finite(draws$b_Intercept)))
  expect_true(all(is.finite(draws$b_abund_Intercept)))
  expect_lt(abs(mean(draws$b_Intercept) - qlogis(0.6)), 0.7)
  expect_lt(abs(mean(draws$b_abund_Intercept) - log(5)), 0.6)
  # Check the fitted-value adapter against native brms predictors using the
  # same small fit, avoiding a second compilation or a packaged fit fixture.
  ids <- c(3L, 1L)
  fitted <- fitted_nmix(fit, draw_ids = ids)
  expect_identical(dim(fitted$abund), c(60L, 2L))
  expect_identical(dim(fitted$det), c(60L, 3L, 2L))
  expected <- t(brms::posterior_linpred(fit, dpar = "abund", draw_ids = ids))
  expect_equal(unname(fitted$abund), exp(expected[d$metadata$unit_start, , drop = FALSE]))
  expected <- t(brms::posterior_linpred(fit, dpar = "mu", draw_ids = ids))
  expect_equal(as.vector(fitted$det[, , 1]), stats::plogis(expected[d$metadata$original_to_packed, 1]))
  expect_identical(names(fitted_nmix(fit, components = "det", draw_ids = 1)), "det")
  nd <- make_nmix_data(matrix(c(1, NA, 2, 0), 2, byrow = TRUE,
                              dimnames = list(c("a", "b"), c("am", "pm"))),
                       newdata_checks = TRUE, quiet = TRUE)
  new <- fitted_nmix(fit, new_data = nd, draw_ids = 1, response = FALSE)
  expect_identical(dim(new$abund), c(2L, 1L))
  expect_identical(dim(new$det), c(2L, 2L, 1L))
  expect_equal(unname(new$abund[, 1]), rep(draws$b_abund_Intercept[1], 2))
  expect_true(is.na(new$det[1, 2, 1]))
  expect_identical(rownames(new$abund), c("a", "b"))
  expect_error(fitted_nmix(fit, draw_ids = 0), "draw_ids")
  expect_error(fitted_nmix(fit, new_data = d$data), "formatted")
  bad <- nd
  bad$metadata$schema_version <- 1L
  expect_error(fitted_nmix(fit, new_data = bad), "Incompatible")
  expect_error(fitted_nmix(fit, new_data = make_flocker_data(matrix(0, 2, 2), quiet = TRUE)),
               "different model type")
})

test_that("detection defaults are independent of counts and user priors override them", {
  d <- nmix_fitting_data()
  large <- make_nmix_data(matrix(20, 4, 3),
                          data.frame(g = factor(c("a", "a", "b", "b"))), quiet = TRUE)
  p <- get_nmix_prior(~ 1, ~ 1 + (1 | g), d)
  q <- get_nmix_prior(~ 1, ~ 1 + (1 | g), large)
  expect_identical(p, q)
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == ""], "student_t(3, 0, 2.5)")
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == "abund"], "")
  expect_true(all(p$prior[p$class == "sd" & p$group == ""] == "student_t(3, 0, 2.5)"))
  expect_match(nmix_stancode(~ 1, ~ 1 + (1 | g), large),
               "student_t_lpdf(sd_1 | 3, 0, 2.5)", fixed = TRUE)
  custom <- c(brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "mu"),
              brms::set_prior("normal(1, 1)", class = "Intercept", dpar = "abund"))
  p <- get_nmix_prior(~ 1, ~ 1, d, prior = custom)
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == ""], "normal(0, 1)")
  expect_equal(p$prior[p$class == "Intercept" & p$dpar == "abund"], "normal(1, 1)")
  expect_silent(nmix_stancode(~ 0 + x, ~ 0 + x, d))
  spline <- get_nmix_prior(~ 1, ~ s(x, k = 3), d)
  expect_equal(spline$prior[spline$class == "sds" & spline$coef == ""],
                "student_t(3, 0, 2.5)")
  # With intercept defaults fixed, data sizes/counts do not alter Stan source.
  expect_identical(nmix_stancode(~ 1, ~ 1, d), nmix_stancode(~ 1, ~ 1, large))
})
