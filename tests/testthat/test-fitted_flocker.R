test_that("fitted_flocker works correctly", {
  # new data
  newdat <- expand.grid(uc1 = rnorm(10), 
                        uc2 = rnorm(10), 
                        ec1 = rnorm(10),
                        species = unique(example_flocker_model_single$data$species))
  
  # check default CI return
  test_out <- fitted_flocker(example_flocker_model_single, summarise = TRUE)
  expect_equal(names(test_out), c("linpred_occ", "linpred_det"))
  
  expect_equal(dimnames(test_out$linpred_occ)[[3]], c("mean", "Q5", "Q95"))
  
  # check varying CIs 
  
  expect_error(fitted_flocker(example_flocker_model_single, summarise = TRUE, CI = c(0, 2)), 
               "CI must be between zero and one inclusive")
  expect_error(fitted_flocker(example_flocker_model_single, summarise = TRUE, CI = c(-1, 1)), 
               "CI must be between zero and one inclusive")
  expect_error(fitted_flocker(example_flocker_model_single, summarise = TRUE, CI = "x"), 
               "CI must be numeric")
  expect_error(fitted_flocker(example_flocker_model_single, summarise = TRUE, CI = 1), 
               "CI must be of length 2")
  test_out <- fitted_flocker(example_flocker_model_single, summarise = TRUE, CI = c(.1, .8))
  expect_equal(dimnames(test_out$linpred_occ)[[3]], c("mean", "Q10", "Q80"))
  
  test_out <- fitted_flocker(example_flocker_model_single, summarise = TRUE, CI = c(.111, .88))
  expect_equal(dimnames(test_out$linpred_occ)[[3]], c("mean", "Q11.1", "Q88"))
  
  # check new_data 
  test_out <- fitted_flocker(example_flocker_model_single, summarise = TRUE, new_data = newdat)
  expect_equal(names(test_out), c("linpred_occ", "linpred_det"))
  expect_equal(names(test_out$linpred_occ), c("mean", "Q5", "Q95"))
  expect_equal(names(test_out$linpred_det), c("mean", "Q5", "Q95"))
  
  expect_equal(nrow(test_out$linpred_occ), nrow(newdat))
  expect_equal(nrow(test_out$linpred_det), nrow(newdat))
  
  new_flocker_dat <- list(data = example_flocker_model_single$data, type = "single")
  class(new_flocker_dat) <- "flocker_data"
  
  test_out2 <- fitted_flocker(example_flocker_model_single, new_data = new_flocker_dat)
  test_out3 <- fitted_flocker(example_flocker_model_single)
  expect_identical(test_out2, test_out3)
  
  # check error for missing predictor column
  expect_error(fitted_flocker(example_flocker_model_single, new_data = newdat[,-1]))
  
  # check non-summary version: 
  test_out <- fitted_flocker(example_flocker_model_single)
  expect_equal(dimnames(test_out$linpred_occ)[[3]], paste0("draw_", 1:brms::ndraws(example_flocker_model_single)))
  expect_equal(
    dim(test_out$linpred_occ), 
    c(
      example_flocker_model_single$data$ff_n_unit[1], 
      max(example_flocker_model_single$data$ff_n_rep),
      brms::ndraws(example_flocker_model_single)))
  
  # check new levels
  newdat2 <- newdat
  levels(newdat2$species) <- c(levels(newdat2$species), "Morphnus_guianensis")
  newdat2$species[1] <- "Morphnus_guianensis"
  expect_error(fitted_flocker(example_flocker_model_single, new_data = newdat2),
               "Levels 'Morphnus_guianensis' of grouping factor 'species' cannot be found in the fitted model")

  expect_equal(names(fitted_flocker(example_flocker_model_single, new_data = newdat2, allow_new_levels = T)),
               c("linpred_occ", "linpred_det"))
})
