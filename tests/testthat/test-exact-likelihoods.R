single_rv_lpmf_from_flocker_data <- function(fd, psi, theta) {
  d <- fd$data
  n_unit <- d$ff_n_unit[1]
  rep_cols <- grep("^ff_rep_index", names(d), value = TRUE)
  ll <- 0

  for (i in seq_len(n_unit)) {
    indices <- as.integer(d[i, rep_cols])
    indices <- indices[indices != -99]
    y_i <- d$ff_y[indices]
    theta_i <- theta[indices]

    if (d$ff_Q[i] == 1) {
      p_i <- psi[i] * prod(theta_i^y_i * (1 - theta_i)^(1 - y_i))
    } else {
      p_i <- (1 - psi[i]) + psi[i] * prod(1 - theta_i)
    }
    ll <- ll + log(p_i)
  }

  ll
}

emission_prob_exact <- function(z, y, theta) {
  keep <- !is.na(y)
  y <- y[keep]
  theta <- theta[keep]

  if (z == 0) {
    return(as.numeric(!any(y == 1)))
  }

  prod(theta^y * (1 - theta)^(1 - y))
}

dynamic_lik_brute_force <- function(init, colo, ex, obs, det) {
  n_year <- dim(obs)[3]
  histories <- expand.grid(rep(list(0:1), n_year))
  prob <- 0

  for (i in seq_len(nrow(histories))) {
    z <- as.integer(histories[i, ])
    p_z <- if (z[1] == 1) init else 1 - init

    if (n_year > 1) {
      for (t in 2:n_year) {
        if (z[t - 1] == 0) {
          p_z <- p_z * if (z[t] == 1) colo[t] else 1 - colo[t]
        } else {
          p_z <- p_z * if (z[t] == 0) ex[t] else 1 - ex[t]
        }
      }
    }

    p_y <- 1
    for (t in seq_len(n_year)) {
      p_y <- p_y * emission_prob_exact(z[t], obs[1, , t], det[1, , t])
    }

    prob <- prob + p_z * p_y
  }

  prob
}

drop_draw <- function(x, draw = 1) {
  array(x[, , , draw], dim = dim(x)[1:3])
}

augmented_lpmf_from_flocker_data <- function(fd, psi, theta, Omega) {
  d <- fd$data
  n_unit <- d$ff_n_unit[1]
  n_sp <- d$ff_n_sp[1]
  rep_cols <- grep("^ff_rep_index", names(d), value = TRUE)
  ll <- 0

  for (sp in seq_len(n_sp)) {
    p_y_given_available <- 1
    unit_rows <- which(d$ff_species[seq_len(n_unit)] == sp)

    for (i in unit_rows) {
      indices <- as.integer(d[i, rep_cols])
      indices <- indices[indices != -99]
      y_i <- d$ff_y[indices]
      theta_i <- theta[indices]

      if (d$ff_Q[i] == 1) {
        p_i <- psi[i] * prod(theta_i^y_i * (1 - theta_i)^(1 - y_i))
      } else {
        p_i <- (1 - psi[i]) + psi[i] * prod(1 - theta_i)
      }

      p_y_given_available <- p_y_given_available * p_i
    }

    if (d$ff_superQ[sp] == 1) {
      ll <- ll + log(Omega * p_y_given_available)
    } else {
      ll <- ll + log((1 - Omega) + Omega * p_y_given_available)
    }
  }

  ll
}

test_that("single-season rep-varying likelihood matches exact marginal likelihood", {
  obs <- matrix(
    c(
      0, 0, NA,
      1, 0, 0
    ),
    nrow = 2,
    byrow = TRUE
  )
  theta <- matrix(
    c(
      0.2, 0.5, NA,
      0.4, 0.6, 0.7
    ),
    nrow = 2,
    byrow = TRUE
  )
  psi <- c(0.3, 0.8)

  fd <- make_flocker_data(obs, event_covs = list(theta = theta), quiet = TRUE)
  gp <- get_positions(fd)

  expected <- log((1 - psi[1]) + psi[1] * (1 - theta[1, 1]) * (1 - theta[1, 2])) +
    log(psi[2] * theta[2, 1] * (1 - theta[2, 2]) * (1 - theta[2, 3]))

  expect_equal(fd$data$ff_y[gp], obs, check.attributes = FALSE)
  expect_equal(fd$data$theta[gp], theta, check.attributes = FALSE)
  expect_equal(
    single_rv_lpmf_from_flocker_data(fd, psi, fd$data$theta),
    expected,
    tolerance = 1e-12
  )
})

test_that("multiseason colex explicit likelihood matches brute force", {
  obs <- array(c(0, 0, 1, 0, 0, 0), dim = c(1, 2, 3))
  det <- array(c(0.2, 0.5, 0.4, 0.6, 0.3, 0.7), dim = c(1, 2, 3, 1))
  init <- matrix(0.35, nrow = 1)
  colo <- array(c(0.11, 0.25, 0.40), dim = c(1, 3, 1))
  ex <- array(c(0.13, 0.30, 0.20), dim = c(1, 3, 1))

  expected <- log(dynamic_lik_brute_force(init[1, 1], colo[1, , 1], ex[1, , 1], obs, drop_draw(det)))
  actual <- log_lik_dynamic(init, colo, ex, obs, det)[1, 1]

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("multiseason colex equilibrium likelihood matches brute force", {
  obs <- array(c(0, 0, 1, 0, 0, 0), dim = c(1, 2, 3))
  det <- array(c(0.2, 0.5, 0.4, 0.6, 0.3, 0.7), dim = c(1, 2, 3, 1))
  colo <- array(c(0.11, 0.25, 0.40), dim = c(1, 3, 1))
  ex <- array(c(0.13, 0.30, 0.20), dim = c(1, 3, 1))
  init <- matrix(colo[1, 1, 1] / (colo[1, 1, 1] + ex[1, 1, 1]), nrow = 1)

  expected <- log(dynamic_lik_brute_force(init[1, 1], colo[1, , 1], ex[1, , 1], obs, drop_draw(det)))
  actual <- log_lik_dynamic(init, colo, ex, obs, det)[1, 1]

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("multiseason autologistic explicit likelihood matches brute force", {
  obs <- array(c(0, 0, 1, 0, 0, 0), dim = c(1, 2, 3))
  det <- array(c(0.2, 0.5, 0.4, 0.6, 0.3, 0.7), dim = c(1, 2, 3, 1))
  init <- matrix(0.35, nrow = 1)
  colo <- array(c(0.11, 0.25, 0.40), dim = c(1, 3, 1))
  auto <- array(c(0.5, 1.2, 0.8), dim = c(1, 3, 1))
  ex <- 1 - boot::inv.logit(boot::logit(colo) + auto)

  expected <- log(dynamic_lik_brute_force(init[1, 1], colo[1, , 1], ex[1, , 1], obs, drop_draw(det)))
  actual <- log_lik_dynamic(init, colo, ex, obs, det)[1, 1]

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("multiseason autologistic equilibrium likelihood matches brute force", {
  obs <- array(c(0, 0, 1, 0, 0, 0), dim = c(1, 2, 3))
  det <- array(c(0.2, 0.5, 0.4, 0.6, 0.3, 0.7), dim = c(1, 2, 3, 1))
  colo <- array(c(0.11, 0.25, 0.40), dim = c(1, 3, 1))
  auto <- array(c(0.5, 1.2, 0.8), dim = c(1, 3, 1))
  ex <- 1 - boot::inv.logit(boot::logit(colo) + auto)
  init <- matrix(colo[1, 1, 1] / (colo[1, 1, 1] + ex[1, 1, 1]), nrow = 1)

  expected <- log(dynamic_lik_brute_force(init[1, 1], colo[1, , 1], ex[1, , 1], obs, drop_draw(det)))
  actual <- log_lik_dynamic(init, colo, ex, obs, det)[1, 1]

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("augmented likelihood matches exact marginal likelihood", {
  obs <- array(NA, dim = c(2, 2, 2))
  obs[, , 1] <- matrix(
    c(
      0, 0,
      1, 0
    ),
    nrow = 2,
    byrow = TRUE
  )
  obs[, , 2] <- matrix(
    c(
      0, 1,
      0, 0
    ),
    nrow = 2,
    byrow = TRUE
  )
  theta <- matrix(
    c(
      0.2, 0.5,
      0.4, 0.6
    ),
    nrow = 2,
    byrow = TRUE
  )
  psi_by_site_species <- matrix(
    c(
      0.3, 0.7, 0.2,
      0.8, 0.4, 0.6
    ),
    nrow = 2,
    ncol = 3,
    byrow = TRUE
  )
  Omega <- 0.65

  fd <- make_flocker_data(
    obs,
    event_covs = list(theta = theta),
    type = "augmented",
    n_aug = 1,
    quiet = TRUE
  )
  gp <- get_positions(fd)
  psi <- psi_by_site_species[as.matrix(gp[, 1, ])]

  observed_species_lik <-
    ((1 - psi_by_site_species[1, 1]) +
       psi_by_site_species[1, 1] * (1 - theta[1, 1]) * (1 - theta[1, 2])) *
    (psi_by_site_species[2, 1] * theta[2, 1] * (1 - theta[2, 2]))
  observed_species_lik2 <-
    (psi_by_site_species[1, 2] * (1 - theta[1, 1]) * theta[1, 2]) *
    ((1 - psi_by_site_species[2, 2]) +
       psi_by_site_species[2, 2] * (1 - theta[2, 1]) * (1 - theta[2, 2]))
  augmented_species_lik <-
    ((1 - psi_by_site_species[1, 3]) +
       psi_by_site_species[1, 3] * (1 - theta[1, 1]) * (1 - theta[1, 2])) *
    ((1 - psi_by_site_species[2, 3]) +
       psi_by_site_species[2, 3] * (1 - theta[2, 1]) * (1 - theta[2, 2]))
  expected <- log(Omega * observed_species_lik) +
    log(Omega * observed_species_lik2) +
    log((1 - Omega) + Omega * augmented_species_lik)

  expect_equal(fd$data$theta[gp[, , 1]], theta, check.attributes = FALSE)
  expect_equal(fd$data$theta[gp[, , 2]], theta, check.attributes = FALSE)
  expect_equal(fd$data$theta[gp[, , 3]], theta, check.attributes = FALSE)
  expect_equal(
    augmented_lpmf_from_flocker_data(fd, psi, fd$data$theta, Omega),
    expected,
    tolerance = 1e-12
  )
})
