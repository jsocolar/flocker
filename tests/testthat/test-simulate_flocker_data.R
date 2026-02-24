expect_named <- function(x, nm) {
  testthat::expect_true(setequal(names(x), nm))
}

expect_array_dim <- function(x, d) {
  testthat::expect_true(is.array(x))
  testthat::expect_equal(dim(x), d)
}

test_that("single-season rep_constant shapes are correct", {
  fd <- simulate_flocker_data(n_pt = 7, n_sp = 3, n_rep = 4, rep_constant = TRUE, seed = 1)
  
  expect_true(is.matrix(fd$obs))
  expect_equal(dim(fd$obs), c(7*3, 4))
  
  expect_true(is.data.frame(fd$unit_covs))
  expect_named(fd$unit_covs, c("uc1","species"))
  expect_equal(nrow(fd$unit_covs), 7*3)
  
  expect_true(is.null(fd$event_covs))
  expect_true(is.list(fd$params))
  expect_true(is.list(fd$covariates))
})

test_that("single-season rep_varying includes event covariates with correct shapes", {
  fd <- simulate_flocker_data(n_pt = 5, n_sp = 2, n_rep = 3, rep_constant = FALSE, seed = 2)
  
  expect_true(is.matrix(fd$obs))
  expect_equal(dim(fd$obs), c(5*2, 3))
  
  expect_true(is.list(fd$event_covs))
  expect_true("ec1" %in% names(fd$event_covs))
  expect_true(is.matrix(fd$event_covs$ec1))
  expect_equal(dim(fd$event_covs$ec1), c(5*2, 3))
})

test_that("seed produces deterministic output", {
  a <- simulate_flocker_data(n_pt=6, n_sp=2, n_rep=3, seed=123)
  b <- simulate_flocker_data(n_pt=6, n_sp=2, n_rep=3, seed=123)
  
  expect_identical(a$obs, b$obs)
  expect_identical(a$unit_covs, b$unit_covs)
  expect_identical(a$event_covs, b$event_covs)
  expect_identical(a$params$coefs, b$params$coefs)
})

test_that("non-NULL seed does not advance global RNG", {
  set.seed(999)
  x1 <- runif(1)
  
  simulate_flocker_data(n_pt=3, n_sp=1, n_rep=2, seed=1)
  
  x2 <- runif(1)
  
  set.seed(999)
  y1 <- runif(1)
  y2 <- runif(1)
  
  expect_identical(x1, y1)
  expect_identical(x2, y2)
})

test_that("NULL seed advances global RNG", {
  set.seed(999)
  x1 <- runif(1)
  
  simulate_flocker_data(n_pt=3, n_sp=1, n_rep=2, seed=NULL)
  
  x2 <- runif(1)
  
  set.seed(999)
  y1 <- runif(1)
  y2 <- runif(1)
  
  expect_identical(x1, y1)
  expect_false(identical(x2, y2))
})

test_that("multiseason/multi_init constraints are enforced", {
  expect_error(simulate_flocker_data(n_season=1, multiseason="colex", multi_init=NULL))
  expect_error(simulate_flocker_data(n_season=1, multiseason=NULL, multi_init="explicit"))
  
  expect_error(simulate_flocker_data(n_season=2, multiseason=NULL, multi_init="explicit"))
  expect_error(simulate_flocker_data(n_season=2, multiseason="colex", multi_init=NULL))
})

test_that("augmented constraints are enforced", {
  expect_error(simulate_flocker_data(n_season=2, augmented=TRUE, multiseason="colex", multi_init="explicit"))
})

test_that("params rejects unknown keys", {
  expect_error(simulate_flocker_data(params=list(zzz=1), seed=1))
})

test_that("params$coefs validation works", {
  expect_error(simulate_flocker_data(n_sp=2, params=list(coefs=matrix(1,2,2)), seed=1))
  
  bad_rows <- data.frame(det_intercept=1, det_slope_unit=1, occ_intercept=1, occ_slope_unit=1)
  expect_error(simulate_flocker_data(n_sp=2, params=list(coefs=bad_rows), seed=1))
  
  bad_col <- data.frame(det_intercept=1, det_slope_unit="a", occ_intercept=1, occ_slope_unit=1)
  expect_error(simulate_flocker_data(n_sp=1, params=list(coefs=bad_col), seed=1))
  
  bad_names <- data.frame(a=1, b=1, c=1, d=1)
  expect_error(simulate_flocker_data(n_sp=1, params=list(coefs=bad_names), seed=1))
})

test_that("cannot supply coefs with coef_means or Sigma", {
  cf <- data.frame(det_intercept=0, det_slope_unit=0, occ_intercept=0, occ_slope_unit=0)
  expect_error(simulate_flocker_data(n_sp=1, params=list(coefs=cf, coef_means=rep(0,4)), seed=1))
  expect_error(simulate_flocker_data(n_sp=1, params=list(coefs=cf, Sigma=diag(4)), seed=1))
})

test_that("covariates validation works", {
  expect_error(simulate_flocker_data(covariates=1, seed=1))
  expect_error(simulate_flocker_data(covariates=list(), seed=1))
  expect_error(simulate_flocker_data(n_pt=3, covariates=list(uc1=1:2), seed=1))
  
  # ec1 length check only when rep_constant=FALSE
  expect_error(simulate_flocker_data(n_pt=3, n_rep=2, n_season=1, rep_constant=FALSE,
                                     covariates=list(uc1=rnorm(3), ec1=rnorm(3)), seed=1))
})


test_that("multiseason output shapes are correct", {
  fd <- simulate_flocker_data(n_pt=4, n_sp=2, n_rep=3, n_season=5,
                              multiseason="colex", multi_init="explicit",
                              rep_constant=FALSE, seed=1)
  
  expect_true(is.array(fd$obs))
  expect_equal(dim(fd$obs), c(4*2, 3, 5))
  
  expect_true(is.list(fd$unit_covs))
  expect_length(fd$unit_covs, 5)
  expect_true(all(vapply(fd$unit_covs, is.data.frame, logical(1))))
  expect_equal(nrow(fd$unit_covs[[1]]), 8)
  
  expect_true(is.list(fd$event_covs))
  expect_true(is.array(fd$event_covs$ec1))
  expect_equal(dim(fd$event_covs$ec1), c(8, 3, 5))
})

test_that("ragged_rep introduces missing visits", {
  fd <- simulate_flocker_data(n_pt=10, n_sp=1, n_rep=6, ragged_rep=TRUE,
                              rep_constant=FALSE, seed=1)
  expect_true(any(is.na(fd$obs)))
  expect_true(any(is.na(fd$event_covs$ec1)))
})

test_that("missing_seasons introduces missing seasons", {
  fd <- simulate_flocker_data(n_pt=10, n_sp=1, n_rep=3, n_season=4,
                              multiseason="colex", multi_init="explicit",
                              rep_constant=FALSE, missing_seasons=TRUE, seed=1)
  # Should be some NA in obs
  expect_true(any(is.na(fd$obs)))
})

test_that("augmented output is trimmed and has expected shapes", {
  fd <- simulate_flocker_data(n_pt=30, n_sp=20, n_rep=4, augmented=TRUE,
                              rep_constant=FALSE, seed=2)
  
  expect_true(is.array(fd$obs))
  expect_equal(dim(fd$obs)[1:2], c(30, 4))
  expect_true(dim(fd$obs)[3] <= 20)      # trimmed, maybe equal but often less
  expect_true(dim(fd$obs)[3] >= 1)
  
  expect_true(is.data.frame(fd$unit_covs))
  expect_named(fd$unit_covs, c("uc1"))
  expect_equal(nrow(fd$unit_covs), 30)
  
  expect_true(is.list(fd$event_covs))
  expect_true(is.matrix(fd$event_covs$ec1))
  expect_equal(dim(fd$event_covs$ec1), c(30, 4))
})

test_that("coef set is correct: single-season rep-constant", {
  fd <- simulate_flocker_data(n_sp=2, rep_constant=TRUE, augmented=FALSE, seed=1)
  expect_true(setequal(colnames(fd$params$coefs),
                       c("det_intercept","det_slope_unit","occ_intercept","occ_slope_unit")))
})

test_that("coef set is correct: multiseason colex equilibrium", {
  fd <- simulate_flocker_data(n_sp=2, n_season=3, multiseason="colex", multi_init="equilibrium",
                              rep_constant=FALSE, seed=1)
  expect_true(setequal(colnames(fd$params$coefs),
                       c("det_intercept","det_slope_unit","det_slope_visit",
                         "col_intercept","col_slope_unit","ex_intercept","ex_slope_unit")))
})

test_that("coef set is correct: multiseason autologistic explicit", {
  fd <- simulate_flocker_data(n_sp=2, n_season=3, multiseason="autologistic", multi_init="explicit",
                              rep_constant=FALSE, seed=1)
  expect_true(setequal(colnames(fd$params$coefs),
                       c("det_intercept","det_slope_unit","det_slope_visit",
                         "occ_intercept","occ_slope_unit",
                         "col_intercept","col_slope_unit",
                         "auto_intercept","auto_slope_unit")))
})

test_that("coef set is correct: multiseason autologistic equilibrium", {
  fd <- simulate_flocker_data(n_sp=2, n_season=3, multiseason="autologistic", multi_init="equilibrium",
                              rep_constant=FALSE, seed=1)
  expect_true(setequal(colnames(fd$params$coefs),
                       c("det_intercept","det_slope_unit","det_slope_visit",
                         "col_intercept","col_slope_unit",
                         "auto_intercept","auto_slope_unit")))
})

test_that("sfd validates covariates$uc1", {
  expect_error(flocker:::sfd(n_rep=2,n_pt=3,n_sp=1,n_season=1,
                             multiseason=NULL,multi_init=NULL,augmented=FALSE,rep_constant=TRUE,
                             params=list(),
                             covariates=list(uc1=1:2),
                             ragged_rep=FALSE,missing_seasons=FALSE))
})

test_that("simulate_flocker_data runs across key configurations", {
  configs <- list(
    list(n_season=1, rep_constant=TRUE, augmented=FALSE),
    list(n_season=1, rep_constant=FALSE, augmented=FALSE),
    list(n_season=1, rep_constant=FALSE, augmented=TRUE),
    list(n_season=3, rep_constant=FALSE, multiseason="colex", multi_init="explicit"),
    list(n_season=3, rep_constant=FALSE, multiseason="colex", multi_init="equilibrium"),
    list(n_season=3, rep_constant=FALSE, multiseason="autologistic", multi_init="explicit"),
    list(n_season=3, rep_constant=FALSE, multiseason="autologistic", multi_init="equilibrium")
  )
  
  for (cfg in configs) {
    fd <- do.call(simulate_flocker_data, c(list(n_pt=5,n_sp=3,n_rep=2,seed=1), cfg))
    expect_true(is.list(fd))
    expect_true("obs" %in% names(fd))
  }
})