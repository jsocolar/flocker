test_that("make_flocker_data works correctly", {
  example_data <- example_flocker_data()
  obs <- example_data$obs
  site_covs <- example_data$site_covs
  visit_covs <- example_data$visit_covs
  
  fd <- make_flocker_data(obs, site_covs, visit_covs)
  expect_equal(fd$.type, "V")
  expect_equal(class(fd), c("list", "flocker_data"))
  expect_equal(names(fd), c("flocker_data", ".nsite", ".nvisit", ".indices", ".type"))
  expect_equal(nrow(fd$flocker_data), sum(fd$.nvisit))
  expect_equal(fd$.indices, matrix(1:12000, ncol = 4))
  expect_equal(names(fd$flocker_data), c("det", "vc1", "vc2", "occ", "sc1", "sc2", "grp", 
                                         "species", ".site", ".visit", "occupancy_subset"))
  expect_identical(fd$flocker_data$occupancy_subset, c(rep(1,3000), rep(0, 9000)))
  expect_true(all(fd$flocker_data$occ %in% c(0,1,-99)))
  
  visit_covs[[1]] <- matrix(1:12000, ncol = 4)
  fd <- make_flocker_data(obs, site_covs, visit_covs)
  expect_equal(fd$flocker_data$vc1, 1:12000)
  
  
  obs[2,4] <- NA
  fd <- make_flocker_data(obs, site_covs, visit_covs)
  expect_equal(fd$.type, "V")
  expect_equal(class(fd), c("list", "flocker_data"))
  expect_equal(names(fd), c("flocker_data", ".nsite", ".nvisit", ".indices", ".type"))
  expect_equal(nrow(fd$flocker_data), sum(fd$.nvisit))
  expect_equal(fd$.indices[,1:3], matrix(1:9000, ncol = 3))
  expect_equal(fd$.indices[,4], c(9001, -99, 9002:11999))
  expect_equal(names(fd$flocker_data), c("det", "vc1", "vc2", "occ", "sc1", "sc2", "grp", 
                                         "species", ".site", ".visit", "occupancy_subset"))
  expect_identical(fd$flocker_data$occupancy_subset, c(rep(1,3000), rep(0, 8999)))
  expect_true(all(fd$flocker_data$occ %in% c(0,1,-99)))
  expect_equal(fd$flocker_data$vc1[1:9000], 1:9000)
  expect_equal(fd$flocker_data$vc1[9001:11999], c(9001, 9003:12000))
  
  obs[,4] <- NA
  expect_warning(fd <- make_flocker_data(obs),
                 "The final column of obs contains only NAs")
  
  obs <- array(obs, dim = c(nrow(obs), ncol(obs), 1))
  expect_error(fd <- make_flocker_data(obs), 
               "obs must have exactly two dimensions")
  obs <- rep(1, 3000)
  expect_error(fd <- make_flocker_data(obs), 
               "obs must have exactly two dimensions")
  obs <- matrix(1:3000, ncol=1)
  expect_error(fd <- make_flocker_data(obs), 
               "obs must contain at least two columns")
  obs <- example_data$obs
  obs[1,1] <- 2
  expect_error(fd <- make_flocker_data(obs), 
               "obs contains values other than 0, 1, NA")
  obs[1,1] <- NA
  expect_error(fd <- make_flocker_data(obs), 
               "obs has NAs in its first column")
  obs[1,1] <- 1
  obs[1,2] <- NA
  expect_error(fd <- make_flocker_data(obs), 
               "Some rows of obs have non-trailing NAs")

  # Still need to add checks for the rest of the error messages.
})