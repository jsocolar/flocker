test_that("make_flocker_data works correctly", {
  example_data <- example_flocker_data
  obs <- example_data$obs
  unit_covs <- example_data$unit_covs
  event_covs <- example_data$event_covs
  
  fd <- make_flocker_data(obs, unit_covs, event_covs)
  expect_equal(fd$type, "single")
  expect_equal(class(fd), c("list", "flocker_data"))
  expect_equal(names(fd), c("data", "n_rep", "type", "unit_covs", "event_covs"))
  expect_equal(
    sum(fd$data[1:900, c("ff_rep_index1", "ff_rep_index2", 
                          "ff_rep_index3", "ff_rep_index4")] - 
          matrix(1:3600, ncol = 4)), 0)  
  expect_equal(names(fd$data), c("ff_y", "uc1", "uc2", "grp", "species", "ec1", "ec2",
                                 "ff_n_unit", "ff_n_rep", "ff_Q", "ff_unit", "ff_rep_index1", 
                                 "ff_rep_index2", "ff_rep_index3", "ff_rep_index4"))
  expect_true(all(fd$data$y %in% c(0,1,-99)))
  
  event_covs[[1]] <- matrix(1:3600, ncol = 4)
  fd <- make_flocker_data(obs, unit_covs, event_covs)
  expect_equal(fd$data$ec1, 1:3600)
  
  
  obs[2,4] <- NA
  fd <- make_flocker_data(obs, unit_covs, event_covs)
  expect_equal(
    sum(fd$data[1:900, c("ff_rep_index1", "ff_rep_index2", "ff_rep_index3")] - 
          matrix(1:2700, ncol = 3)), 0)
  expect_equal(fd$data$ff_rep_index4[1:900], c(2701, -99, 2702:3599))
  expect_true(all(fd$data$ff_y %in% c(0,1)))
  expect_equal(fd$data$ec1[1:2700], 1:2700)
  expect_equal(fd$data$ec1[2701:3599], c(2701, 2703:3600))
  
  obs[,4] <- NA
  expect_error(fd <- make_flocker_data(obs),
                 "The final column of obs contains only NAs.")
  
  obs <- array(obs, dim = c(nrow(obs), ncol(obs), 1))
  expect_error(fd <- make_flocker_data(obs), 
               "in a single-season model, obs must have exactly two dimensions")
  obs <- rep(1, 3000)
  expect_error(fd <- make_flocker_data(obs), 
               "in a single-season model, obs must have exactly two dimensions")
  obs <- matrix(1:3000, ncol=1)
  expect_error(fd <- make_flocker_data(obs), 
               "obs must contain at least two columns.")
  obs <- example_data$obs
  obs[1,1] <- 2
  expect_error(fd <- make_flocker_data(obs))
  obs[1,1] <- NA
  expect_error(fd <- make_flocker_data(obs), 
               "obs has NAs in its first column")
  obs[1,1] <- 1
  obs[1,2] <- NA
  expect_error(fd <- make_flocker_data(obs), 
               "Some rows of obs have non-trailing NAs")

  # Still need to add checks for the rest of the error messages, and for proper error messages using rep-constant data
})
