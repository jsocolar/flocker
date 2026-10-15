make_ragged_multi_fixture <- function(global_trailing_season = FALSE) {
  n_year <- if(global_trailing_season) 6 else 5
  obs <- array(NA_real_, dim = c(3, 4, n_year))
  obs[1, , 1] <- c(1, 0, 1, 0)
  obs[2, 1:3, 1] <- c(0, 1, 0)
  obs[3, 1:2, 1] <- c(1, 0)
  obs[2, 1, 2] <- 1
  obs[3, , 2] <- c(0, 1, 0, 1)
  obs[1, 1:2, 3] <- c(1, 0)
  obs[3, 1:3, 3] <- c(0, 0, 1)
  obs[1, 1, 4] <- 0
  obs[2, , 4] <- c(1, 0, 1, 0)
  obs[1, 1:3, 5] <- c(0, 1, 0)
  obs[3, 1, 5] <- 1

  event <- array(seq_along(obs) / 100 - 0.4, dim = dim(obs))
  event[is.na(obs)] <- NA
  unit_values <- matrix(
    seq_len(3 * n_year) / 20 - 0.5,
    nrow = 3
  )
  unit_covs <- lapply(seq_len(n_year), function(year) {
    data.frame(uc1 = unit_values[, year])
  })
  expected_unit <- unit_values
  expected_unit[2, 5] <- NA
  if(global_trailing_season) {
    expected_unit[, 6] <- NA
  }

  list(
    obs = obs,
    event = event,
    unit_covs = unit_covs,
    expected_unit = expected_unit
  )
}


format_ragged_multi_fixture <- function(global_trailing_season = FALSE) {
  fixture <- make_ragged_multi_fixture(global_trailing_season)
  suppressWarnings(make_flocker_data(
    fixture$obs,
    unit_covs = fixture$unit_covs,
    event_covs = list(ec1 = fixture$event),
    type = "multi",
    quiet = TRUE
  ))
}
