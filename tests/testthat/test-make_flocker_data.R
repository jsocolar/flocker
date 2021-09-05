test_that("make_flocker_data works correctly", {
  example_data <- example_flocker_data()
  obs <- example_data$obs
  unit_covs <- example_data$unit_covs
  event_covs <- example_data$event_covs
  
  fd <- make_flocker_data(obs, unit_covs, event_covs)
  expect_equal(fd$type, "V")
  expect_equal(class(fd), c("list", "flocker_data"))
  expect_equal(names(fd), c("data", "max_rep", "type"))
  expect_equal(
    sum(fd$data[1:3000, c("rep_index1", "rep_index2", 
                          "rep_index3", "rep_index4")] - 
          matrix(1:12000, ncol = 4)), 0)  
  expect_equal(names(fd$data), c("y", "uc1", "uc2", "grp", "species", "ec1", "ec2",
                                 "n_unit", "n_rep", "Q", "unit", "rep_index1", 
                                 "rep_index2", "rep_index3", "rep_index4"))
  expect_true(all(fd$data$occ %in% c(0,1,-99)))
  
  event_covs[[1]] <- matrix(1:12000, ncol = 4)
  fd <- make_flocker_data(obs, unit_covs, event_covs)
  expect_equal(fd$data$ec1, 1:12000)
  
  
  obs[2,4] <- NA
  fd <- make_flocker_data(obs, unit_covs, event_covs)
  expect_equal(fd$type, "V")
  expect_equal(class(fd), c("list", "flocker_data"))
  expect_equal(names(fd), c("data", "max_rep", "type"))
  expect_equal(
    sum(fd$data[1:3000, c("rep_index1", "rep_index2", "rep_index3")] - 
          matrix(1:9000, ncol = 3)), 0)
  expect_equal(fd$data$rep_index4[1:3000], c(9001, -99, 9002:11999))
  expect_equal(names(fd$data), c("y", "uc1", "uc2", "grp", "species", "ec1", "ec2",
                                 "n_unit", "n_rep", "Q", "unit", "rep_index1", 
                                 "rep_index2", "rep_index3", "rep_index4"))
  expect_true(all(fd$flocker_data$occ %in% c(0,1,-99)))
  expect_equal(fd$data$ec1[1:9000], 1:9000)
  expect_equal(fd$data$ec1[9001:11999], c(9001, 9003:12000))
  
  obs[,4] <- NA
  expect_error(fd <- make_flocker_data(obs),
                 "The final column of obs contains only NAs.")
  
  obs <- array(obs, dim = c(nrow(obs), ncol(obs), 1))
  expect_error(fd <- make_flocker_data(obs), 
               "obs must have exactly two dimensions.")
  obs <- rep(1, 3000)
  expect_error(fd <- make_flocker_data(obs), 
               "obs must have exactly two dimensions.")
  obs <- matrix(1:3000, ncol=1)
  expect_error(fd <- make_flocker_data(obs), 
               "obs must contain at least two columns.")
  obs <- example_data$obs
  obs[1,1] <- 2
  expect_error(fd <- make_flocker_data(obs), 
               "obs contains values other than 0, 1, NA.")
  obs[1,1] <- NA
  expect_error(fd <- make_flocker_data(obs), 
               "obs has NAs in its first column.")
  obs[1,1] <- 1
  obs[1,2] <- NA
  expect_error(fd <- make_flocker_data(obs), 
               "Some rows of obs have non-trailing NAs")

  # Still need to add checks for the rest of the error messages, and for proper error messages using rep-constant data
})
