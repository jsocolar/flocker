test_that("complete-block slices preserve occupancy likelihoods and gradients", {
  skip_on_cran()
  skip_on_os("windows")
  source <- paste(make_occupancy_twolevel_single_lpmf(),
                  make_occupancy_multi_lpmf(),
                  # The augmented wrapper delegates to the two-level function.
                  paste(readLines(system.file("stan", "occupancy_augmented.stan",
                    package = "flocker")), collapse = "\n"), sep = "\n")
  calls <- c(
    "occupancy_single_lpmf(y[a:b] | mu[a:b], occ[a:b], unit_size[a:b], Q[a:b])",
    "occupancy_multi_colex_lpmf(y[a:b] | mu[a:b], occ[a:b], colo[a:b], ex[a:b], outer_size[a:b], unit_size[a:b], Q[a:b])",
    "occupancy_multi_colex_eq_lpmf(y[a:b] | mu[a:b], colo[a:b], ex[a:b], outer_size[a:b], unit_size[a:b], Q[a:b])",
    "occupancy_multi_autologistic_lpmf(y[a:b] | mu[a:b], occ[a:b], colo[a:b], ex[a:b], outer_size[a:b], unit_size[a:b], Q[a:b])",
    "occupancy_multi_autologistic_eq_lpmf(y[a:b] | mu[a:b], colo[a:b], ex[a:b], outer_size[a:b], unit_size[a:b], Q[a:b])",
    "occupancy_twolevel_single_lpmf(y[a:b] | mu[a:b], occ[a:b], Omega[a:b], outer_size[a:b], unit_size[a:b], Q[a:b], group_Q[a:b])",
    "occupancy_augmented_lpmf(y[a:b] | mu[a:b], occ[a:b], Omega[a:b], outer_size[a:b], unit_size[a:b], Q[a:b], group_Q[a:b])"
  )
  code <- paste0("functions {", source, "}", '
    data {
      int<lower=1> n; int<lower=1,upper=7> family;
      array[n] int y; array[n] int unit_size; array[n] int outer_size;
      array[n] int Q; array[n] int group_Q;
      int<lower=1> n_slice; array[n_slice] int ends; vector[n] x;
    }
    parameters { vector[5] beta; }
    model {
      vector[n] mu = rep_vector(beta[1], n) + x;
      vector[n] occ = rep_vector(beta[2], n) + 0.2*x;
      vector[n] colo = rep_vector(beta[3], n) - 0.3*x;
      vector[n] ex = rep_vector(beta[4], n) + 0.4*x;
      vector[n] Omega = rep_vector(beta[5], n) - 0.1*x;
      int a = 1;
      for (slice in 1:n_slice) {
        int b = ends[slice];',
        paste0(ifelse(seq_along(calls) == 1, "if", "else if"),
          " (family == ", seq_along(calls), ") target += ", calls, ";", collapse = "\n"),
      'a = b + 1;
      }
    }')
  model <- rstan::stan_model(model_code = code)
  single <- make_flocker_data(rbind(c(1, 0, 0), c(0, NA, NA), c(0, 0, NA), c(1, 1, NA)),
    event_covs = list(w = matrix(0, 4, 3)), quiet = TRUE)
  multi <- format_ragged_multi_fixture()
  two <- make_flocker_data(rbind(c(1, 0), c(0, NA), c(0, 0), c(0, NA), c(1, 1)),
    data.frame(g = factor(c("a", "a", "b", "b", "c"))), type = "twolevel_single",
    level2_group = "g", quiet = TRUE)
  augmented <- make_flocker_data(array(rep(c(1, 0, 0, 0, 0, 0), 2), c(3, 2, 2)),
    type = "augmented", n_aug = 2, quiet = TRUE)
  for (family in seq_along(calls)) {
    fd <- if (family == 1) single else if (family <= 5) multi else if (family == 6) two else augmented
    d <- fd$data
    outer <- if (family == 1) d$ff_unit_size else if (family <= 5) d$ff_series_size else d$ff_group_size
    n <- nrow(d)
    dat <- list(n = n, family = family, y = d$ff_y, unit_size = d$ff_unit_size,
      outer_size = outer, Q = d$ff_Q,
      group_Q = if (family >= 6) d$ff_group_Q else rep(0L, n),
      x = seq(-.8, .8, length.out = n), n_slice = 1L, ends = as.array(n))
    full <- rstan::sampling(model, data = dat, chains = 0)
    starts <- which(outer > 0)
    ends <- starts + outer[starts] - 1L
    # Every possible two-slice partition, plus one slice per complete block.
    partitions <- c(lapply(ends[-length(ends)], function(end) c(end, n)), list(ends))
    for (partition in partitions) {
      dat$n_slice <- length(partition)
      dat$ends <- as.array(partition)
      sliced <- rstan::sampling(model, data = dat, chains = 0)
      for (beta in list(c(-.5, .3, -.8, .4, .2), c(2, -2, 1, -1, .7))) {
        expect_equal(rstan::log_prob(sliced, beta), rstan::log_prob(full, beta), tolerance = 1e-10)
        expect_equal(rstan::grad_log_prob(sliced, beta), rstan::grad_log_prob(full, beta), tolerance = 1e-10)
      }
    }
  }
})
