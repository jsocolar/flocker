test_that("fitted_flocker works correctly", {
  # new data
  newdat <- expand.grid(uc1 = rnorm(10), 
                        uc2 = rnorm(10), 
                        ec1 = rnorm(10),
                        species = unique(fixture_fit("single")$data$species))
  
  # check default CI return
  test_out <- fitted_flocker(fixture_fit("single"), summarise = TRUE)
  expect_equal(names(test_out), c("linpred_occ", "linpred_det"))
  
  expect_equal(dimnames(test_out$linpred_occ)[[3]], c("mean", "Q5", "Q95"))
  
  # check varying CIs 
  
  expect_error(fitted_flocker(fixture_fit("single"), summarise = TRUE, CI = c(0, 2)),
               "CI must be between zero and one inclusive")
  expect_error(fitted_flocker(fixture_fit("single"), summarise = TRUE, CI = c(-1, 1)),
               "CI must be between zero and one inclusive")
  expect_error(fitted_flocker(fixture_fit("single"), summarise = TRUE, CI = "x"),
               "CI must be numeric")
  expect_error(fitted_flocker(fixture_fit("single"), summarise = TRUE, CI = 1),
               "CI must be of length 2")
  test_out <- fitted_flocker(fixture_fit("single"), summarise = TRUE, CI = c(.1, .8))
  expect_equal(dimnames(test_out$linpred_occ)[[3]], c("mean", "Q10", "Q80"))
  
  test_out <- fitted_flocker(fixture_fit("single"), summarise = TRUE, CI = c(.111, .88))
  expect_equal(dimnames(test_out$linpred_occ)[[3]], c("mean", "Q11.1", "Q88"))
  
  # check new_data 
  test_out <- fitted_flocker(fixture_fit("single"), summarise = TRUE, new_data = newdat)
  expect_equal(names(test_out), c("linpred_occ", "linpred_det"))
  expect_equal(names(test_out$linpred_occ), c("mean", "Q5", "Q95"))
  expect_equal(names(test_out$linpred_det), c("mean", "Q5", "Q95"))
  
  expect_equal(nrow(test_out$linpred_occ), nrow(newdat))
  expect_equal(nrow(test_out$linpred_det), nrow(newdat))
  
  new_flocker_dat <- list(
    data = fixture_fit("single")$data,
    metadata = get_flocker_metadata(fixture_fit("single"))
  )
  class(new_flocker_dat) <- "flocker_data"
  
  test_out2 <- fitted_flocker(fixture_fit("single"), new_data = new_flocker_dat)
  test_out3 <- fitted_flocker(fixture_fit("single"))
  expect_identical(test_out2, test_out3)
  
  # check error for missing predictor column
  expect_error(fitted_flocker(fixture_fit("single"), new_data = newdat[,-1]))
  
  # check non-summary version: 
  test_out <- fitted_flocker(fixture_fit("single"))
  expect_equal(dimnames(test_out$linpred_occ)[[3]], paste0("draw_", 1:brms::ndraws(fixture_fit("single"))))
  expect_equal(
    dim(test_out$linpred_occ), 
    c(
      get_flocker_metadata(fixture_fit("single"))$original_dim,
      brms::ndraws(fixture_fit("single"))))
  
  # check new levels
  newdat2 <- newdat
  levels(newdat2$species) <- c(levels(newdat2$species), "Morphnus_guianensis")
  newdat2$species[1] <- "Morphnus_guianensis"
  expect_error(fitted_flocker(fixture_fit("single"), new_data = newdat2),
               "Levels 'Morphnus_guianensis' of grouping factor 'species' cannot be found in the fitted model")

  expect_equal(names(fitted_flocker(fixture_fit("single"), new_data = newdat2, allow_new_levels = T)),
               c("linpred_occ", "linpred_det"))
})

test_that("Omega follows the standard fitted_flocker return shape", {
  testthat::skip_on_cran()

  observation_level <- fitted_flocker(
    fixture_fit("aug"),
    components = c("occ", "Omega"),
    draw_ids = 1:2
  )
  expect_identical(
    dim(observation_level$linpred_Omega),
    dim(observation_level$linpred_occ)
  )
  expect_identical(attr(observation_level, "unit_level"), FALSE)

  unit_level <- fitted_flocker(
    fixture_fit("aug"),
    components = c("occ", "Omega"),
    draw_ids = 1:2,
    unit_level = TRUE
  )
  expect_identical(
    dim(unit_level$linpred_Omega),
    dim(unit_level$linpred_occ)
  )
  expect_identical(attr(unit_level, "unit_level"), TRUE)

  one_draw <- fitted_flocker(
    fixture_fit("aug"),
    components = c("occ", "Omega"),
    draw_ids = 1
  )
  expect_identical(
    dim(one_draw$linpred_Omega),
    dim(one_draw$linpred_occ)
  )
})
