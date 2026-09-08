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
    group = factor(c("b", "a", "b", "a"), levels = c("a", "b")),
    unit_x = 1:4
  )
  event_covs <- list(event_x = matrix(seq_len(12), nrow = 4))
  group_covs <- data.frame(
    group = factor(c("b", "a"), levels = c("a", "b")),
    group_x = c(20, 10),
    known = c(FALSE, FALSE)
  )
  
  fd <- make_flocker_data(
    obs, unit_covs, event_covs, type = "twolevel_single",
    group_covs = group_covs, top_level = "group",
    known_present = "known", quiet = TRUE
  )
  
  expect_equal(fd$type, "twolevel_single")
  expect_equal(fd$top_level, "group")
  expect_equal(fd$max_unit_group, 2)
  expect_equal(fd$group_covs, names(group_covs))
  expect_equal(fd$data$group[seq_len(fd$data$ff_n_group[1])], factor(c("a", "b"), levels = c("a", "b")))
  expect_equal(fd$data$group_x[seq_len(fd$data$ff_n_group[1])], c(10, 20))
  expect_equal(fd$data$ff_group_known_present[seq_len(fd$data$ff_n_group[1])], c(0, 1))
  expect_equal(fd$data$ff_n_unit_group[seq_len(fd$data$ff_n_group[1])], c(2, 2))
  expect_true(all(c("ff_group_index1", "ff_group_index2") %in% names(fd$data)))
  expect_equal(
    fd$data$ff_orig_unit[seq_len(fd$data$ff_n_unit[1])],
    c(2, 1, 3, 4)
  )
  
  unit_covs_bad <- unit_covs
  unit_covs_bad$group <- as.character(unit_covs_bad$group)
  expect_error(
    make_flocker_data(
      obs, unit_covs_bad, event_covs, type = "twolevel_single",
      group_covs = group_covs, top_level = "group", quiet = TRUE
    ),
    "top_level must identify a factor column"
  )
  
  group_covs_bad <- transform(group_covs, unit_x = c(1, 2))
  expect_error(
    make_flocker_data(
      obs, unit_covs, event_covs, type = "twolevel_single",
      group_covs = group_covs_bad, top_level = "group", quiet = TRUE
    ),
    "may only share the top_level column"
  )
})
