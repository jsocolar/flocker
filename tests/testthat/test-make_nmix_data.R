test_that("N-mixture formatting retains visits and round-trips counts and covariates", {
  obs <- rbind(a = c(2, 3, 0), b = c(0, 1, NA), c = c(5, NA, NA))
  colnames(obs) <- c("morning", "noon", "evening")
  unit <- data.frame(x = c(10, 20, 30))
  w <- matrix(seq_len(9), 3)
  w[is.na(obs)] <- NA
  fd <- make_nmix_data(obs, unit, list(w = w), quiet = TRUE)
  expect_named(fd, c("data", "metadata", "auxiliary"))
  expect_identical(fd$data$ff_y, c(2L, 3L, 0L, 0L, 1L, 5L))
  expect_equal(fd$data$x, c(10, 10, 10, 20, 20, 30))
  expect_equal(fd$data$w, c(1, 4, 7, 2, 5, 3))
  expect_identical(fd$metadata$type, "nmix")
  expect_identical(fd$metadata$data_role, "model_data")
  expect_identical(fd$metadata$schema_version, 2L)
  expect_equal(matrix(fd$data$ff_y[get_positions(fd)], nrow(obs),
                      dimnames = dimnames(obs)), obs)
  expect_equal(unname(get_positions(fd, unit_level = TRUE)), c(1L, 4L, 6L))
  expect_identical(make_nmix_data(as.data.frame(obs), unit, list(w = as.data.frame(w)),
                                quiet = TRUE), fd)
  # The same formatter branch produces the same row order and real covariates.
  occupancy <- make_flocker_data(1L * (obs > 0), unit, list(w = w), quiet = TRUE)
  expect_identical(fd$data[names(fd$data) != "ff_y"],
                   occupancy$data[names(occupancy$data) != "ff_y"])
  expect_identical(get_positions(fd), get_positions(occupancy))
  expect_identical(fd$auxiliary,
    nmix_binomial_auxiliary_data(list(c(2L, 3L, 0L), c(0L, 1L), 5L)))
  expect_identical(fd$auxiliary$ff_unit_offsets, c(1L, 4L, 6L, 7L))
  expect_identical(fd$auxiliary$ff_coef_offsets, c(1L, 4L, 5L, 6L))
  expect_false(any(grepl("coef", names(fd$data))))

  # Formatter output connects directly to the existing Stan-data declarations.
  sd <- brms::make_standata(brms::bf(ff_y ~ 1, abund ~ 1), data = fd$data,
    family = nmix_poisson_binomial(),
    stanvars = nmix_poisson_binomial_stanvars(fd$auxiliary))
  for (name in names(fd$auxiliary)) expect_equal(as.vector(sd[[name]]), fd$auxiliary[[name]])
  expect_equal(as.vector(sd$Y), fd$data$ff_y)
})

test_that("N-mixtures never collapse visits, including binary and all-zero responses", {
  for (obs in list(matrix(0, 2, 2), rbind(c(0, 1), c(1, 0)), rbind(c(3, 2), c(0, 4)))) {
    fd <- make_nmix_data(obs, quiet = TRUE)
    expect_identical(fd$metadata$type, "nmix")
    expect_identical(fd$data$ff_y, as.integer(t(obs)))
    expect_identical(fd$data$ff_Q, rep(as.integer(rowSums(obs) > 0), each = ncol(obs)))
    expect_equal(nrow(fd$data), length(obs))
    expect_false("ff_n_suc" %in% names(fd$data))
  }
  obs <- rbind(c(0, 1), c(1, 0))
  expect_identical(make_flocker_data(obs, quiet = TRUE)$metadata$type, "single_C")
  expect_error(make_flocker_data(obs + 2, quiet = TRUE), "0, 1, or NA")
  expect_error(flocker_standata(~ 1, ~ 1, make_nmix_data(obs, quiet = TRUE)),
               "requires occupancy data")
})

test_that("count-specific validation and shared structural checks run before packing", {
  for (value in list(-1, 0.5, Inf, -Inf, NaN, "2", .Machine$integer.max)) {
    obs <- matrix(0, 2, 2)
    obs[1, 1] <- value
    expect_error(make_nmix_data(obs, quiet = TRUE), "nonnegative integer")
  }
  expect_error(make_nmix_data(rbind(c(0, NA), c(NA, NA)), quiet = TRUE), "first column")
  expect_error(make_nmix_data(rbind(c(2, NA, 1), c(1, 2, 0)), quiet = TRUE), "non-trailing")
  expect_error(make_nmix_data(matrix(1, 2, 2), list(x = 1), quiet = TRUE))
  expect_error(make_nmix_data(matrix(1, 2, 2), event_covs = list(w = matrix(NA, 2, 2)),
                             quiet = TRUE), "event covariate contains missing")
  expect_error(make_nmix_data(array(1, c(2, 2, 2)), quiet = TRUE), "two dimensions")
  expect_error(make_flocker_data_(matrix(1, 2, 2), type = "multi", model = "nmix"),
               "single-season, one-level")
  expect_error(make_nmix_data(matrix(c(2, 0), 2), quiet = TRUE), "at least two columns")
  fd <- make_nmix_data(matrix(c(2, 0), 2), quiet = TRUE, newdata_checks = TRUE)
  expect_identical(fd$metadata$data_role, "new_data")
  expect_equal(fd$data$ff_y, c(2, 0))
  expect_message(make_nmix_data(matrix(0, 2, 2)), "N-mixture")
})
