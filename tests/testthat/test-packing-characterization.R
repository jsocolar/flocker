# Characterize ecological values rather than the layout of the 1.1 index
# columns. Only this adapter should need replacing for the v2 schema.
recover_packed_values <- function(fd, variable, unit_level = FALSE) {
  positions <- get_positions(fd, unit_level = unit_level)
  values <- fd$data[[variable]][positions]
  if (is.null(dim(positions))) return(unname(values))
  array(values, dim = dim(positions))
}

test_that("ragged two-level packing keeps all ecological levels aligned", {
  obs <- rbind(c(0, 1, 0), c(0, NA, NA), c(1, 0, NA),
               c(0, 0, 1), c(0, NA, NA))
  dimnames(obs) <- list(paste0("unit-", 1:5), c("early", "middle", "late"))
  unit_covs <- data.frame(
    group = factor(c("birch", "ash", "birch", "cedar", "ash"),
                   levels = c("cedar", "ash", "birch")),
    unit_marker = c(51, 12, 83, 24, 95)
  )
  # Deliberately different input-row, factor-level, and unit encounter orders.
  level2_covs <- data.frame(
    group = factor(c("birch", "cedar", "ash"),
                   levels = c("cedar", "ash", "birch")),
    group_marker = c(300, 100, 200)
  )
  event_marker <- matrix(seq_len(15), 5)
  event_marker[is.na(obs)] <- NA
  fd <- make_flocker_data(
    obs, unit_covs, list(event_marker = event_marker),
    type = "twolevel_single", level2_group = "group",
    level2_covs = level2_covs, quiet = TRUE
  )

  expect_equal(recover_packed_values(fd, "ff_y"), unname(obs))
  expect_equal(recover_packed_values(fd, "event_marker"), event_marker)
  expect_equal(recover_packed_values(fd, "unit_marker", TRUE),
               c(51, 12, 83, 24, 95))
  expect_equal(recover_packed_values(fd, "group_marker", TRUE),
               c(300, 200, 300, 100, 200))
  expected_unit_events <- matrix(rep(unit_covs$unit_marker, 3), 5)
  expected_unit_events[is.na(obs)] <- NA
  expect_equal(recover_packed_values(fd, "unit_marker"), expected_unit_events)
  expect_identical(dim(get_positions(fd)), c(5L, 3L))
  expect_identical(get_flocker_metadata(fd)$level2_group_names, c("cedar", "ash", "birch"))
  expect_identical(get_flocker_metadata(fd)$unit_covs,
                   c("group", "unit_marker"))
  expect_identical(get_flocker_metadata(fd)$event_covs, "event_marker")
})

test_that("aggregated single-season packing preserves ragged sufficient statistics", {
  obs <- rbind(c(1, 0, 1), c(0, NA, NA), c(1, 0, NA))
  fd <- make_flocker_data(
    obs, unit_covs = data.frame(unit_marker = c(17, 31, 59)), quiet = TRUE
  )
  expect_identical(fd$type, "single_C")
  expect_equal(fd$data$ff_n_trial, c(3, 1, 2))
  expect_equal(fd$data$ff_n_suc, c(2, 0, 1))
  expect_equal(recover_packed_values(fd, "unit_marker", TRUE), c(17, 31, 59))
  # The 1.1 aggregated path repeats unit rows even at unsampled visits.
  # This is characterization, not the v2 prediction-availability contract.
  expect_equal(recover_packed_values(fd, "unit_marker"),
               matrix(rep(c(17, 31, 59), 3), 3))
})

test_that("single-unit and single-visit event maps retain matrix dimensions", {
  for (obs in list(matrix(c(1, 0), 1), matrix(c(1, 0), 2))) {
    marker <- matrix(seq_along(obs), nrow(obs))
    fd <- make_flocker_data(
      obs, event_covs = list(event_marker = marker),
      newdata_checks = TRUE, quiet = TRUE
    )
    expect_identical(dim(get_positions(fd)), dim(obs))
    expect_equal(recover_packed_values(fd, "ff_y"), obs)
    expect_equal(recover_packed_values(fd, "event_marker"), marker)
    expect_length(get_positions(fd, unit_level = TRUE), nrow(obs))
  }
})

test_that("augmented species names align with reconstructed site and event values", {
  obs <- array(c(1, 0, 0, NA, 0, 1, 0, NA), c(2, 2, 2),
               dimnames = list(c("north", "south"), c("am", "pm"),
                               c("wren", "thrush")))
  event_marker <- matrix(c(11, 22, 33, NA), 2)
  fd <- make_flocker_data(
    obs, unit_covs = data.frame(site_marker = c(71, 92)),
    event_covs = list(event_marker = event_marker),
    type = "augmented", n_aug = 1, quiet = TRUE
  )
  expected_obs <- array(NA_real_, c(2, 2, 3))
  expected_obs[, , 1:2] <- obs
  expected_obs[, , 3] <- matrix(c(0, 0, 0, NA), 2)
  expect_identical(get_flocker_metadata(fd)$level2_group_names, c("wren", "thrush", "pseudospecies1"))
  expect_equal(recover_packed_values(fd, "ff_y"), expected_obs)
  expect_equal(recover_packed_values(fd, "event_marker"),
               array(rep(event_marker, 3), c(2, 2, 3)))
  expect_equal(recover_packed_values(fd, "site_marker", TRUE),
               matrix(rep(c(71, 92), 3), 2))
  expect_identical(dim(get_positions(fd)), c(2L, 2L, 3L))
  expect_identical(dim(get_positions(fd, TRUE)), c(2L, 3L))
})

test_that("dynamic reconstruction distinguishes internal gaps from trailing units", {
  fixture <- make_ragged_multi_fixture(global_trailing_season = TRUE)
  # Give trailing units distinctive, complete covariates: they still must not
  # enter fitting rows. Internal gaps, in contrast, retain their unit predictor.
  fixture$unit_covs[[5]]$uc1[2] <- 999
  fixture$unit_covs[[6]]$uc1[] <- 1001:1003
  fd <- suppressWarnings(make_flocker_data(
    fixture$obs, fixture$unit_covs, list(event_marker = fixture$event),
    type = "multi", quiet = TRUE
  ))
  units <- recover_packed_values(fd, "uc1", TRUE)
  expect_equal(units, fixture$expected_unit)
  expect_equal(units[1, 2], fixture$unit_covs[[2]]$uc1[1])
  expect_true(all(is.na(recover_packed_values(fd, "ff_y")[1, , 2])))
  expect_false(any(fd$data$uc1 %in% c(999, 1001:1003)))
  expect_equal(recover_packed_values(fd, "ff_y"), fixture$obs)
  expect_equal(recover_packed_values(fd, "event_marker"), fixture$event)
  expect_identical(dim(get_positions(fd)), c(3L, 4L, 6L))
  expect_identical(dim(get_positions(fd, TRUE)), c(3L, 6L))
})
