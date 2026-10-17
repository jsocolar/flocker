test_that("single-season packing preserves identities and real covariates", {
  obs <- rbind(c(0, 1, 0), c(0, NA, NA), c(1, 0, NA))
  dimnames(obs) <- list(c("ash", "birch", "cedar"), c("am", "noon", "pm"))
  w <- matrix(seq_len(9), 3, dimnames = dimnames(obs))
  unit <- data.frame(x = c(10, 20, 30), group = factor(c("a", "b", "a")))
  fd <- make_flocker_data(obs, unit, list(w = w), quiet = TRUE)
  expect_named(fd, c("data", "metadata"))
  expect_named(fd$data, c("ff_y", "x", "group", "w", "ff_n_unit", "ff_unit_start"))
  expect_equal(fd$data$ff_y, c(0L, 1L, 0L, 0L, 1L, 0L))
  expect_equal(fd$data$ff_n_unit, c(3L, rep(-99L, 5)))
  expect_equal(fd$data$ff_unit_start, c(1L, 4L, 5L, rep(-99L, 3)))
  expect_equal(fd$data$x, c(10, 10, 10, 20, 30, 30))
  expect_equal(as.character(fd$data$group), c("a", "a", "a", "b", "a", "a"))
  m <- get_flocker_metadata(fd)
  expect_identical(m$schema_version, 2L)
  expect_identical(m$data_role, "model_data")
  expect_equal(m$unit_start, c(1L, 4L, 5L))
  expect_equal(m$event_unit, c(1L, 1L, 1L, 2L, 3L, 3L))
  restored <- array(NA_real_, m$original_dim, m$original_dimnames)
  restored[m$packed_to_original] <- fd$data$ff_y
  expect_equal(restored, obs)
  expect_equal(fd$data$ff_y[m$original_to_packed[!is.na(obs)]], obs[!is.na(obs)])
  expect_equal(fd$data$w, w[m$packed_to_original])
  from_frame <- make_flocker_data(as.data.frame(obs), unit, list(w = as.data.frame(w)), quiet = TRUE)
  expect_identical(from_frame$data, fd$data)

  sd <- flocker_standata(~ x, ~ x + w, fd)
  expect_equal(as.vector(sd$Y), fd$data$ff_y)
  expect_equal(as.vector(sd$vint1), fd$data$ff_n_unit)
  expect_equal(as.vector(sd$vint2), fd$data$ff_unit_start)
  expect_equal(unname(sd$X_occ[, "x"]), fd$data$x)
  code <- flocker_stancode(~ x, ~ x + w, fd)
  expect_match(code, make_occupancy_single_lpmf(), fixed = TRUE)
  expect_true(rstan::stanc(model_code = code)$status)
  # The public combined-formula route uses the same packed input.
  combined <- flocker_standata(f_det = brms::bf(det ~ x + w, occ ~ x), flocker_data = fd)
  expect_equal(combined, sd)
  expect_error(flocker_stancode(~ w, ~ w, fd), "occupancy formula")
  expect_error(flocker_stancode(~ log(w), ~ w, fd), "occupancy formula")
  expect_error(flocker_stancode(f_det = brms::bf(det ~ w, occ ~ log(w)),
                               flocker_data = fd), "occupancy formula")
})

test_that("single-season fitting rejects incompatible packing and roles early", {
  obs <- rbind(c(1, 0), c(0, NA))
  fd <- make_flocker_data(obs, event_covs = list(w = matrix(1:4, 2)), quiet = TRUE)
  for (tag in list(NULL, 1L, 999L)) {
    bad <- fd
    bad$metadata$schema_version <- tag
    expect_error(flock(~ 1, ~ w, bad), "Incompatible single-season packing schema")
  }
  old <- structure(list(data = fd$data, type = "single", n_rep = 2L),
                   class = c("list", "flocker_data"))
  expect_error(flock(~ 1, ~ w, old), "Pre-2.0")
  prediction <- make_flocker_data(obs, event_covs = list(w = matrix(1:4, 2)),
                                  newdata_checks = TRUE, quiet = TRUE)
  expect_identical(prediction$metadata$data_role, "new_data")
  expect_error(flock(~ 1, ~ w, prediction), "new_data cannot")
  expect_error(flock(~ 1, ~ w, fd, threads = 2), "Threading is not supported")
  bad <- fd
  bad$data$ff_n_unit[1] <- 99L
  expect_error(flock(~ 1, ~ w, bad), "Invalid single-season event blocks")
  bad <- fd
  bad$data$ff_unit_start[2] <- 1L
  expect_error(flock(~ 1, ~ w, bad), "Invalid single-season event blocks")
  fd$data$ff_unit_start[1] <- 0L
  expect_error(flock(~ 1, ~ w, fd), "Invalid single-season event blocks")
})

test_that("single-season generated programs are independent of maximum visits", {
  codes <- lapply(c(2L, 5L), function(n) {
    fd <- make_flocker_data(matrix(0L, 2, n),
                            event_covs = list(w = matrix(seq_len(2 * n), 2)),
                            quiet = TRUE)
    flocker_stancode(~ 1, ~ w, fd)
  })
  expect_identical(codes[[1]], codes[[2]])
  fd <- make_flocker_data(matrix(c(1, 0), 1),
                          event_covs = list(w = matrix(c(2, 3), 1)), quiet = TRUE)
  expect_identical(fd$metadata$original_dim, c(1L, 2L))
  expect_equal(fd$data$ff_n_unit, c(1L, -99L))
  expect_equal(fd$data$ff_unit_start, c(1L, -99L))
})
