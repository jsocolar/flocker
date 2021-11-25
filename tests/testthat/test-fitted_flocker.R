test_that("fitted_flocker works correctly", {
  # read test-case flocker fit
  test_fit <- readRDS("tests/test_data/test_fit.rds")

  # new data
  newdat <- expand.grid(uc1 = rnorm(10), 
                        uc2 = rnorm(10), 
                        species = unique(test_fit$data$species))
  
  # check default return
  test_out <- fitted_flocker(test_fit, type = "occupancy")
  expect_equal(names(test_out), c("estimate", "Q5", "Q95"))
  # dimensions
  expect_equal(nrow(test_out), nrow(test_fit$data))
  
  # check varying CIs 
  expect_error(fitted_flocker(test_fit, type = "occupancy", CI = c(0, 2)), 
               "CI cannot have bounds <0 or >1")
  expect_error(fitted_flocker(test_fit, type = "occupancy", CI = c(-1, 1)), 
               "CI cannot have bounds <0 or >1")
  expect_error(fitted_flocker(test_fit, type = "occupancy", CI = "x"), 
               "CI should be numeric length 2")
  expect_error(fitted_flocker(test_fit, type = "occupancy", CI = 1), 
               "CI should be numeric length 2")
  test_out <- fitted_flocker(test_fit, type = "occupancy", CI = c(.1, .8))
  expect_equal(names(test_out), c("estimate", "Q10", "Q80"))
  
  # check new_data 
  test_out <- fitted_flocker(test_fit, type = "occupancy", new_data = newdat)
  expect_equal(names(test_out), c("estimate", "Q5", "Q95"))
  expect_equal(nrow(test_out), nrow(newdat))
  
  test_out2 <- fitted_flocker(test_fit, type = "occupancy", new_data = test_fit$data)
  test_out3 <- fitted_flocker(test_fit, type = "occupancy")
  expect_identical(test_out2, test_out3)
  
  # check error for missing predictor column
  expect_error(fitted_flocker(test_fit, type = "occupancy", new_data = newdat[,-1]), 
               "The following variables can neither be found in 'data' nor in 'data2':\n'uc1'")
  
  # check non-summary version: 
  test_out <- fitted_flocker(test_fit, type = "occupancy", summarise = F)
  expect_equal(colnames(test_out), paste0("iter_", 1:brms::ndraws(test_fit)))
  expect_equal(dim(test_out), c(nrow(test_fit$data), brms::ndraws(test_fit)))
  
  # check new levels
  newdat2 <- expand.grid(uc1 = rnorm(10),
                        uc2 = rnorm(10),
                        species = "A")
  expect_error(fitted_flocker(test_fit, type = "occupancy", new_data = newdat2),
               "Levels 'A' of grouping factor 'species' cannot be found in the fitted model")

  expect_equal(names(fitted_flocker(test_fit, type = "occupancy", new_data = newdat2, allow_new_levels = T)),
               c("estimate", "Q5", "Q95"))

  # to add: detection and both checks
})
