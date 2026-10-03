test_that("make_flocker_data works correctly", {
  example_data <- simulate_flocker_data()
  obs <- example_data$obs
  unit_covs <- example_data$unit_covs
  event_covs <- example_data$event_covs
  
  fd <- make_flocker_data(obs, unit_covs, event_covs, quiet = TRUE)
  expect_equal(fd$type, "single")
  expect_equal(class(fd), c("list", "flocker_data"))
  expect_equal(names(fd), c("data", "n_rep", "type", "unit_covs", "event_covs"))
  expect_true(
    all(fd$data[1:1500, c("ff_rep_index1", "ff_rep_index2", 
                          "ff_rep_index3", "ff_rep_index4")] == 
          matrix(1:6000, ncol = 4))
    )  
  expect_equal(names(fd$data), c("ff_y", "uc1", "species", "ec1",
                                 "ff_n_unit", "ff_n_rep", "ff_Q", "ff_unit", "ff_rep_index1", 
                                 "ff_rep_index2", "ff_rep_index3", "ff_rep_index4"))
  expect_true(all(fd$data$y %in% c(0,1,-99)))
  
  event_covs[[1]] <- matrix(1:6000, ncol = 4)
  fd <- make_flocker_data(obs, unit_covs, event_covs, quiet = TRUE)
  expect_equal(fd$data$ec1, 1:6000)
  
  
  obs[2,4] <- NA
  fd <- make_flocker_data(obs, unit_covs, event_covs, quiet = TRUE)
  expect_equal(
    sum(fd$data[1:1500, c("ff_rep_index1", "ff_rep_index2", "ff_rep_index3")] - 
          matrix(1:4500, ncol = 3)), 0)
  expect_equal(fd$data$ff_rep_index4[1:1500], c(4501, -99, 4502:5999))
  expect_true(all(fd$data$ff_y %in% c(0,1)))
  expect_equal(fd$data$ec1[1:4500], 1:4500)
  expect_equal(fd$data$ec1[4501:5999], c(4501, 4503:6000))
  
  obs[,4] <- NA
  expect_error(fd <- make_flocker_data(obs, quiet = TRUE),
                 "The final column of obs contains only NAs.")
  
  obs <- array(obs, dim = c(nrow(obs), ncol(obs), 1))
  expect_error(fd <- make_flocker_data(obs, quiet = TRUE), 
               "in a single-season model, obs must have exactly two dimensions")
  obs <- rep(1, 3000)
  expect_error(fd <- make_flocker_data(obs, quiet = TRUE), 
               "in a single-season model, obs must have exactly two dimensions")
  obs <- matrix(rep(1, 3000), ncol=1)
  expect_error(fd <- make_flocker_data(obs, quiet = TRUE), 
               "obs must contain at least two columns.")
  obs <- example_data$obs
  obs[1,1] <- 2
  expect_error(fd <- make_flocker_data(obs, quiet = TRUE))
  obs[1,1] <- NA
  expect_error(fd <- make_flocker_data(obs, quiet = TRUE), 
               "obs has NAs in its first column")
  obs[1,1] <- 1
  obs[1,2] <- NA
  expect_error(fd <- make_flocker_data(obs, quiet = TRUE), 
               "Some rows of obs have non-trailing NAs")

  # Still need to add checks for the rest of the error messages, and for proper error messages using rep-constant data
})

test_that("make_flocker_data handles two-level single-season data", {
  obs <- matrix(c(
    1, 0, 0,
    0, 0, NA,
    0, 1, 0,
    0, 0, 0
  ), nrow = 4, byrow = TRUE)
  unit_covs <- data.frame(
    group = factor(
      c("b", "a", "b", "a"),
      levels = c("a", "b", "unused")
    ),
    unit_x = 1:4
  )
  event_covs <- list(event_x = matrix(seq_len(12), nrow = 4))
  level2_covs <- data.frame(
    group = ordered(c("b", "a"), levels = c("a", "b")),
    group_x = c(20, 10)
  )
  
  fd <- make_flocker_data(
    obs, unit_covs, event_covs, type = "twolevel_single",
    level2_covs = level2_covs, level2_group = "group", quiet = TRUE
  )
  
  expect_equal(fd$type, "twolevel_single")
  expect_equal(fd$level2_group, "group")
  expect_equal(fd$max_unit_group, 2)
  expect_equal(fd$level2_covs, "group_x")
  expect_equal(
    as.character(fd$data$group[seq_len(fd$data$ff_n_group[1])]),
    c("a", "b")
  )
  expect_equal(fd$data$group_x[seq_len(fd$data$ff_n_group[1])], c(10, 20))
  expect_equal(fd$data$ff_group_known_present[seq_len(fd$data$ff_n_group[1])], c(0, 1))
  expect_equal(fd$data$ff_n_unit_group[seq_len(fd$data$ff_n_group[1])], c(2, 2))
  expect_true(all(c("ff_group_index1", "ff_group_index2") %in% names(fd$data)))
  expect_equal(fd$unit_order, c(2, 1, 3, 4))
  expect_equal(get_unit_group(fd$data), c(1, 2, 2, 1))
  expect_false("unit_group" %in% names(fd))
  expect_false(any(c("ff_orig_unit", "ff_group") %in% names(fd$data)))
  expect_equal(
    fd$data$group_x[seq_len(fd$data$ff_n_unit[1])],
    c(10, 20, 20, 10)
  )

  fd_minimal <- make_flocker_data(
    obs, unit_covs, event_covs, type = "twolevel_single",
    level2_group = "group", quiet = TRUE
  )
  expect_equal(fd_minimal$data$ff_n_group[1], 2)
  expect_equal(fd_minimal$level2_covs, character(0))
  expect_equal(
    as.character(fd_minimal$data$group[seq_len(fd_minimal$data$ff_n_group[1])]),
    c("a", "b")
  )
  
  unit_covs_bad <- unit_covs
  unit_covs_bad$group <- as.character(unit_covs_bad$group)
  expect_error(
    make_flocker_data(
      obs, unit_covs_bad, event_covs, type = "twolevel_single",
      level2_covs = level2_covs, level2_group = "group", quiet = TRUE
    ),
    "level2_group must identify a factor column"
  )
  
  level2_covs_bad <- transform(level2_covs, unit_x = c(1, 2))
  expect_error(
    make_flocker_data(
      obs, unit_covs, event_covs, type = "twolevel_single",
      level2_covs = level2_covs_bad, level2_group = "group", quiet = TRUE
    ),
    "may only share the level2_group column"
  )

  level2_covs_empty <- data.frame(
    group = factor(c("c", "b", "a"), levels = c("a", "b", "c")),
    group_x = c(30, 20, 10)
  )
  expect_error(
    make_flocker_data(
      obs, unit_covs, event_covs, type = "twolevel_single",
      level2_covs = level2_covs_empty, level2_group = "group", quiet = TRUE
    ),
    "groups with no corresponding units: c"
  )
})

test_that("augmented data reject original species without detections", {
  obs <- array(0, dim = c(2, 2, 2))
  obs[1, 1, 2] <- 1
  expect_error(
    make_flocker_data(
      obs, type = "augmented", n_aug = 1, quiet = TRUE
    ),
    "must contain only species with at least one detection"
  )
})

test_that("augmented internal identifiers use reserved ff_ names", {
  obs <- array(0, dim = c(2, 2, 2))
  obs[1, 1, 1] <- 1
  obs[2, 1, 2] <- 1
  site_covs <- data.frame(species = c(10, 20), site_id = c(100, 200))

  fd <- make_flocker_data(
    obs,
    unit_covs = site_covs,
    type = "augmented",
    n_aug = 1,
    quiet = TRUE
  )
  gp <- get_positions(fd)
  expect_true(all(c("species", "site_id", "ff_species") %in% names(fd$data)))
  expect_false(any(c("ff_orig_unit", "ff_group", "ff_site") %in% names(fd$data)))
  expect_false(any(c("ff_n_sp", "ff_superQ", "ff_known_present") %in% names(fd$data)))
  expect_equal(
    fd$data$ff_group_known_present[seq_len(fd$data$ff_n_group[1])],
    c(1, 1, 0)
  )
  for (sp in seq_len(dim(gp)[3])) {
    expect_equal(
      matrix(as.integer(fd$data$ff_species[gp[, , sp]]), nrow = 2),
      matrix(sp, nrow = 2, ncol = 2)
    )
    expect_equal(
      matrix(as.numeric(fd$data$species[gp[, , sp]]), nrow = 2),
      matrix(site_covs$species, nrow = 2, ncol = 2)
    )
  }
  unit_rows <- seq_len(fd$data$ff_n_unit[1])
  orig_unit <- fd$unit_order
  expect_equal(fd$unit_site, rep(1:2, 3)[orig_unit])
  expect_equal(get_unit_group(fd$data), rep(1:3, each = 2)[orig_unit])
  expect_equal(
    as.integer(fd$data$ff_species[unit_rows]),
    rep(1:3, each = 2)[orig_unit]
  )

  event_species <- matrix(seq_len(4), nrow = 2)
  fd_event <- make_flocker_data(
    obs,
    event_covs = list(species = event_species),
    type = "augmented",
    n_aug = 1,
    quiet = TRUE
  )
  gp_event <- get_positions(fd_event)
  for (sp in seq_len(dim(gp_event)[3])) {
    expect_equal(
      matrix(as.numeric(fd_event$data$species[gp_event[, , sp]]), nrow = 2),
      event_species
    )
  }

  expect_error(
    make_flocker_data(
      obs,
      unit_covs = data.frame(ff_species = 1:2),
      type = "augmented",
      n_aug = 1,
      quiet = TRUE
    ),
    "reserved string"
  )
  expect_error(
    make_flocker_data(
      obs,
      event_covs = list(ff_site = matrix(seq_len(4), nrow = 2)),
      type = "augmented",
      n_aug = 1,
      quiet = TRUE
    ),
    "reserved string"
  )
})
