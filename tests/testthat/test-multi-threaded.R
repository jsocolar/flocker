threaded_multi_args <- function(threads = 1) {
  list(f_occ = ~ uc1, f_det = ~ ec1, f_col = ~ uc1, f_ex = ~ uc1,
       flocker_data = format_ragged_multi_fixture(), multiseason = "colex",
       multi_init = "explicit", threads_per_chain = threads)
}

test_that("threaded colex code keeps the vectorized family and data layout", {
  args <- threaded_multi_args(2)
  code <- do.call(flocker_stancode, args)
  expect_true(rstan::stanc(model_code = code)$status)
  expect_match(code, "target += occupancy_multi_colex_threaded_lpmf(Y", fixed = TRUE)
  expect_match(code, "reduce_sum(occupancy_multi_colex_partial_sum", fixed = TRUE)
  expect_false(grepl("partial_log_lik", code, fixed = TRUE))
  sd <- do.call(flocker_standata, args)
  args$threads_per_chain <- 1
  expect_identical(sd, do.call(flocker_standata, args))
  expect_false(grepl("reduce_sum", do.call(flocker_stancode, args), fixed = TRUE))
  serial_code <- do.call(flocker_stancode, args)
  withr::local_options(list(brms.threads = 2))
  expect_identical(do.call(flocker_stancode, args), serial_code)
  expect_error(do.call(flocker_stancode, c(args, list(threads = 2))),
               "Use threads_per_chain")
  for (bad in list(NULL, NA, 0, -1, 1.5, Inf, TRUE, c(1, 2))) {
    args["threads_per_chain"] <- list(bad)
    expect_error(do.call(flocker_stancode, args), "threads_per_chain must")
  }
  args$threads_per_chain <- 2
  args$f_occ <- NULL
  args$multi_init <- "equilibrium"
  expect_error(do.call(flocker_stancode, args), "Threading is only supported")
  args$multi_init <- "explicit"
  args$f_occ <- ~ uc1
  args$multiseason <- "autologistic"
  args$f_ex <- NULL
  args$f_auto <- ~ uc1
  expect_error(do.call(flocker_stancode, args), "Threading is only supported")
})

# One compiled model exercises both implementations on identical parameters.
# This checks autodiff as well as values, including ragged/unvisited seasons.
test_that("threaded and serial likelihood values and gradients agree", {
  skip_on_cran()
  skip_on_os("windows")
  code <- paste0("functions {", make_occupancy_multi_lpmf(TRUE), "}
    data {
      int<lower=1> N;
      array[N] int y;
      array[N] int n_series;
      array[N] int series_starts;
      array[N] int n_unit;
      array[N] int unit_starts;
      array[N] int Q;
      int parallel;
    }
    parameters { vector[N] mu; vector[N] occ; vector[N] colo; vector[N] ex; }
    model {
      if (parallel == 1)
        target += occupancy_multi_colex_threaded_lpmf(y | mu, occ, colo, ex,
          n_series, series_starts, n_unit, unit_starts, Q);
      else
        target += occupancy_multi_colex_lpmf(y | mu, occ, colo, ex,
          n_series, series_starts, n_unit, unit_starts, Q);
    }")
  old <- rstan::rstan_options("threads_per_chain")
  on.exit(rstan::rstan_options(threads_per_chain = old))
  rstan::rstan_options(threads_per_chain = 2L)
  model <- rstan::stan_model(model_code = code)
  for (fd in list(format_ragged_multi_fixture(), make_flocker_data(
      array(c(0, 1), c(1, 1, 2)), event_covs = list(w = array(0, c(1, 1, 2))),
      type = "multi", quiet = TRUE, newdata_checks = TRUE))) {
    d <- fd$data
    data <- list(N = nrow(d), y = d$ff_y, n_series = d$ff_n_series,
                 series_starts = d$ff_series_start, n_unit = d$ff_n_unit,
                 unit_starts = d$ff_unit_start, Q = d$ff_Q, parallel = 0L)
    serial <- rstan::sampling(model, data = data, chains = 0)
    data$parallel <- 1L
    threaded <- rstan::sampling(model, data = data, chains = 0)
    for (shift in c(-3, 0, 3)) {
      pars <- seq(-2, 2, length.out = 4 * nrow(d)) + shift
      expect_equal(rstan::log_prob(threaded, pars), rstan::log_prob(serial, pars),
                   tolerance = 1e-10)
      expect_equal(rstan::grad_log_prob(threaded, pars),
                   rstan::grad_log_prob(serial, pars), tolerance = 1e-10)
    }
  }
})

test_that("threaded fits retain serial multiseason post-processing", {
  skip_on_cran()
  skip_on_os("windows")
  args <- threaded_multi_args(2)
  args$chains <- 1
  args$cores <- 1
  args$iter <- 80
  args$warmup <- 40
  args$seed <- 197
  args$refresh <- 0
  args$backend <- "rstan"
  args$prior <- c(brms::set_prior("normal(0, 1)", class = "Intercept"),
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "occ"),
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "colo"),
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "ex"),
    brms::set_prior("normal(0, 1)", class = "b"),
    brms::set_prior("normal(0, 1)", class = "b", dpar = "occ"),
    brms::set_prior("normal(0, 1)", class = "b", dpar = "colo"),
    brms::set_prior("normal(0, 1)", class = "b", dpar = "ex"))
  fit <- suppressWarnings(do.call(flock, args))
  expect_identical(fit$family$name, "occupancy_multi_colex_threaded")
  expect_identical(type_flocker_fit(fit), "multi_colex")
  expect_true(fit$threads$force)
  expect_equal(fit$threads$threads, 2)
  serial <- fit
  serial$family <- occupancy_multi_colex()
  for (new_data in list(NULL, args$flocker_data)) {
    expect_equal(get_positions(fit), get_positions(serial))
    expect_equal(fitted_flocker(fit, new_data = new_data, draw_ids = 1:2),
                 fitted_flocker(serial, new_data = new_data, draw_ids = 1:2))
    expect_equal(get_Z(fit, new_data = new_data, draw_ids = 1:2),
                 get_Z(serial, new_data = new_data, draw_ids = 1:2))
    expect_equal(log_lik_flocker(fit, new_data = new_data, draw_ids = 1:2),
                 log_lik_flocker(serial, new_data = new_data, draw_ids = 1:2))
    set.seed(198)
    predicted <- predict_flocker(fit, new_data = new_data, draw_ids = 1:2)
    set.seed(198)
    expect_equal(predicted, predict_flocker(serial, new_data = new_data, draw_ids = 1:2))
  }
  expect_s3_class(suppressWarnings(loo_flocker(fit)), "loo")
})
