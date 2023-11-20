sd <- simulate_flocker_data()
fd_single <- make_flocker_data(sd$obs, sd$unit_covs, sd$event_covs)
fd_single_C <- make_flocker_data(sd$obs, sd$unit_covs)

sd <- simulate_flocker_data(augmented = TRUE)
fd_augmented <- make_flocker_data(sd$obs, sd$unit_covs, sd$event_covs, type = "augmented", n_aug = 10)

sd <- simulate_flocker_data(n_season = 3, multiseason = "colex", multi_init = "explicit")
fd_multi <- make_flocker_data(sd$obs, sd$unit_covs, sd$event_covs, type = "multi")


test_that("flocker_stancode works as expected", {
  # This is effectively a set of tests for `flock` itself, but can be run
  # without actually doing model fitting in Stan
  
  f_occ <- ~ uc1
  f_det <- ~ uc1 + ec1
  f_col <- NULL
  f_ex <- NULL
  f_auto <- NULL
  flocker_data <- fd_single
  data2 <- NULL
  multiseason <- NULL
  multi_init <- NULL
  augmented <- FALSE
  threads <- NULL
  
  expect_type(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                  f_col, f_ex, multi_init, f_auto, augmented, threads),
              "character")
  
  
  f_occ <- ~ uc1 + ec1
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  
  f_occ <- y ~ uc1
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  
  flocker_data <- fd_single_C
  f_occ <- ~ uc1
  f_det <- ~ uc1
  
  expect_type(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads),
              "character")
  
  
  flocker_data <- fd_augmented
  f_det <- ~ uc1 + ec1
  augmented <- TRUE
  
  expect_type(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads),
              "character")
  
  
  flocker_data <- fd_multi
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  augmented <- FALSE
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  multiseason <- "colex"
  multi_init <- "explicit"
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_col <- ~ uc1
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_ex <- ~ uc1
  
  expect_type(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads),
              "character")
  
  multiseason <- "autologistic"
  multi_init <- "equilibrium"
  
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_auto <- ~ uc1
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_auto <- NULL
  f_occ <- NULL
  f_ex <- NULL
  expect_error(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_auto <- ~ uc1
  expect_type(flocker_stancode(f_occ, f_det, flocker_data, data2, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads),
              "character")

})

