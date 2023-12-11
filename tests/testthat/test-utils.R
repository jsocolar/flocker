test_that("log_inv_logit handles scalar input", {
  expect_equal(log_inv_logit(0), log(0.5))
  expect_equal(log_inv_logit(10), log(1 / (1 + exp(-10))))
})

test_that("log_inv_logit handles vector input", {
  input <- c(-3, 0, 3)
  expected_output <- log(1 / (1 + exp(-input)))
  expect_equal(log_inv_logit(input), expected_output)
})

test_that("log_inv_logit handles extreme values", {
  expect_equal(log_inv_logit(100), log(1 / (1 + exp(-100))))
  expect_equal(log_inv_logit(-100), log(1 / (1 + exp(100))))
  expect_false(identical(log_inv_logit(100), log(1 / (1 + exp(-100)))))
  expect_true(is.finite(log_inv_logit(-10000)))
})

test_that("log_inv_logit throws error for non-numeric input", {
  expect_error(log_inv_logit("abc"), "x must be numeric")
  expect_error(log_inv_logit(TRUE), "x must be numeric")
})


test_that("log1m_inv_logit returns correct output", {
  r <- rnorm(100)
  expect_identical(log_inv_logit(r), log1m_inv_logit(-r))
})

test_that("array utils work properly", {
  testmat <- matrix(c(NA, 1:3), nrow = 2)
  testdf <- as.data.frame(testmat)
  testarray <- array(1:8, dim = c(2,2,2))
  
  # expand_matrix
  expect_error(expand_matrix(1:4))
  expect_identical(expand_matrix(testmat), c(NA, 1:3))
  expect_identical(expand_matrix(testdf), c(NA, 1:3))
  
  # expand_array_3D
  expect_error(expand_array_3D(testarray[,,1]))
  expect_equal(expand_array_3D(testarray), matrix(c(1,2,5,6,3,4,7,8), ncol = 2))
  
  # nslice
  expect_error(nslice(testmat))
  expect_equal(nslice(testarray), 2)
  
  # stack_matrix
  expect_error(stack_matrix(testarray, 2))
  expect_equal(stack_matrix(testmat, 2), matrix(c(NA, 1, NA, 1, 2, 3, 2, 3), ncol = 2))
  expect_equivalent(as.matrix(stack_matrix(testdf, 2)), matrix(c(NA, 1, NA, 1, 2, 3, 2, 3), ncol = 2))
  
  # new_matrix
  m <- matrix(1:4, nrow = 2, ncol = 2)
  new_m <- new_matrix(m)
  expect_true(is.matrix(new_m))
  expect_equal(dim(new_m), dim(m))
  expect_true(all(is.na(new_m)))
  
  data_vector <- c(5, 6, 7, 8)
  new_m <- new_matrix(m, data = data_vector)
  expect_equal(new_m, matrix(data_vector, nrow = 2, ncol = 2))
  
  new_m <- new_matrix(m, data = data_vector, byrow = TRUE)
  expect_equal(new_m, matrix(data_vector, nrow = 2, ncol = 2, byrow = TRUE))
  
  m <- matrix(1:9, nrow = 3, ncol = 3)
  data_vector <- 1:2 # insufficient length
  expect_warning(new_m <- new_matrix(m, data = data_vector))
  expected_matrix <- matrix(rep(data_vector, length.out = 9), nrow = 3, ncol = 3)
  expect_equal(new_m, expected_matrix)
  
  m <- matrix(1:4, nrow = 2, ncol = 2)
  data_vector <- 1:5 # excess length
  expect_warning(new_matrix(m, data = data_vector))
  
  data_vector <- c("a", "b", "c", "d")
  new_m <- new_matrix(m, data = data_vector)
  expect_equal(new_m, matrix(data_vector, nrow = 2, ncol = 2))
  
  m <- matrix(numeric(0), nrow = 0, ncol = 0)
  new_m <- new_matrix(m)
  expect_equal(dim(new_m), c(0, 0))

  expect_error(new_matrix(10))
  expect_error(new_matrix(NULL))
  expect_error(new_matrix(NA))
  
  m <- matrix(1:4, nrow = 2, ncol = 2, dimnames = list(c("r1", "r2"), c("c1", "c2")))
  new_m <- new_matrix(m)
  expect_false(identical(dimnames(new_m), dimnames(m)))
  expect_true(is.null(dimnames(new_m)))
  
  # new array
  m <- array(1:8, dim = c(2, 2, 2))
  new_a <- new_array(m)
  expect_true(is.array(new_a))
  expect_equal(dim(new_a), dim(m))
  expect_true(all(is.na(new_a)))
  
  data_vector <- c(9, 10, 11, 12, 13, 14, 15, 16)
  new_a <- new_array(m, data = data_vector)
  expect_equal(new_a, array(data_vector, dim = c(2, 2, 2)))
  
  m <- array(1:27, dim = c(3, 3, 3))
  data_vector <- 1:5 # insufficient length
  new_a <- new_array(m, data = data_vector)
  expected_array <- array(rep(data_vector, length.out = 27), dim = c(3, 3, 3))
  expect_equal(new_a, expected_array)
  
  non_array_input <- 10
  expect_error(new_array(non_array_input))
  
  expect_error(new_array(NULL))
  expect_error(new_array(NA))
})

test_that("bookkeeping works properly", {
  # flocker_col_names
  expect_true(all(grepl("^ff_", flocker_col_names(3, 4))))
  expect_equal(7, length(flocker_col_names(3, 4)) - length(flocker_col_names()))
  
  # flocker_reserved
  expect_true(all(grepl(flocker_reserved()[1], flocker_col_names())))
  expect_true(all(grepl(flocker_reserved()[2], paste0(".", c(".", "foo", 1:2)))))
  
  # flocker_model_types
  expect_true(all(grepl("^single|^augmented|^multi", flocker_model_types())))
  
  # flocker_data_input_types
  expect_true(
    all(
      grepl(
        paste0("^", 
               paste(flocker_data_input_types(), collapse = "|^")
               ),
        flocker_model_types()
        )
      )
    )
  for(i in seq_along(flocker_data_input_types())){
    expect_false(
      all(
        grepl(
          paste0("^", 
                 paste(flocker_data_input_types()[-i], collapse = "|^")
          ),
          flocker_model_types()
        )
      )
    )
  }
  
  # flocker_data_output_types
  expect_true(all(grepl("^single|^augmented|^multi", flocker_data_output_types())))
})

test_that("fdtl function returns expected dataframe", {
  # Call the fdtl function
  result <- fdtl()
  
  # Check if the result is a dataframe
  expect_is(result, "data.frame")
  
  # Check if the result has the correct column names
  expect_named(result, c("model_type", "data_output_type", "data_input_type"))
  
  # Check if the result has the correct number of rows (assuming 10 model types)
  expect_equal(nrow(result), 7)
  
  # Check if the result has the correct number of columns
  expect_equal(ncol(result), 3)
  
  # Check if the model_type column contains the expected values
  expect_identical(result$model_type, flocker_model_types())
  
  # Check if the data_output_type and data_input_type columns contain the expected values
  expected_data_input_types <- c(
    "single", "single", "augmented", "multi", "multi", "multi", "multi"
  )
  
  expected_data_output_types <- c(
    "single", "single_C", "augmented", "multi", "multi", "multi", "multi"
  )
  
  
  expect_identical(result$data_output_type, expected_data_output_types)
  expect_identical(result$data_input_type, expected_data_input_types)
})

test_that("is_flocker_fit works", {
  expect_true(is_flocker_fit(example_flocker_model_single))
  expect_false(is_flocker_fit("foo"))
  expect_false(is_flocker_fit(NULL))
  expect_false(is_flocker_fit(list(f = example_flocker_model_single)))
})

test_that("type_flocker_fit function returns expected string", {
  expect_identical(type_flocker_fit(example_flocker_model_single), "single")
  
  # Create a dummy flocker_fit object
  dummy_flocker_fit <- structure(
    list(),
    class = "flocker_fit",
    data_type = "multi",
    multiseason = "colex",
    multi_init = "equilibrium"
  )
  
  # Call the type_flocker_fit function
  result <- type_flocker_fit(dummy_flocker_fit)
  
  # Check if the result is a character string
  expect_is(result, "character")
  
  # Check if the result has the correct value
  expected_value <- "multi_colex_eq"
  expect_identical(result, expected_value)
  
  # Check if the function throws an error for non-flocker_fit objects
  non_flocker_fit <- list()
  expect_error(type_flocker_fit(non_flocker_fit), "x must be a flocker_fit object")
  
  # Check if the function throws an error for objects with missing or altered attributes
  corrupted_flocker_fit <- structure(
    list(),
    class = "flocker_fit",
    data_type = "MT1"
  )
  expect_error(type_flocker_fit(corrupted_flocker_fit), "the attributes of the flocker_fit object have been altered or corrupted")
})


test_that("get_positions works properly", {
  
  # single-season rep-varying
  sd <- simulate_flocker_data()
  fd <- make_flocker_data(
    sd$obs, sd$unit_covs, sd$event_covs, 
    type = "single")
  ps <- get_positions(fd)
  expect_true(
    all.equal(
      new_array(sd$event_covs$ec1, fd$data$ec1[ps]),
      sd$event_covs$ec1,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  expect_true(
    all.equal(
      fd$data$uc1[ps2],
      sd$unit_covs$uc1,
      check.attributes = FALSE
    )
  )
  
  # single-season rep-varying with missingness
  sd <- simulate_flocker_data(ragged_rep = TRUE)
  fd <- make_flocker_data(
    sd$obs, sd$unit_covs, sd$event_covs, 
    type = "single")
  ps <- get_positions(fd)
  expect_true(
    all.equal(
      new_array(sd$event_covs$ec1, fd$data$ec1[ps]),
      sd$event_covs$ec1,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  expect_true(
    all.equal(
      fd$data$uc1[ps2],
      sd$unit_covs$uc1,
      check.attributes = FALSE
    )
  )
  
  # single_season rep-constant
  sd <- simulate_flocker_data(rep_constant = TRUE)
  fd <- make_flocker_data(
    sd$obs, sd$unit_covs, 
    type = "single")
  ps <- get_positions(fd)
  expect_true(
    all.equal(
      fd$data$uc1[ps[,1]],
      sd$unit_covs$uc1,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  expect_true(
    all.equal(
      fd$data$uc1[ps2],
      sd$unit_covs$uc1,
      check.attributes = FALSE
    )
  )
  
  # single-season rep-constant with missingness
  sd <- simulate_flocker_data(rep_constant = TRUE, ragged_rep = TRUE)
  fd <- make_flocker_data(
    sd$obs, sd$unit_covs, 
    type = "single")
  ps <- get_positions(fd)
  expect_true(
    all.equal(
      fd$data$uc1[ps[,1]],
      sd$unit_covs$uc1,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  expect_true(
    all.equal(
      fd$data$uc1[ps2],
      sd$unit_covs$uc1,
      check.attributes = FALSE
    )
  )
  
  # augmented
  sd <- simulate_flocker_data(augmented = TRUE)
  fd <- make_flocker_data(
    sd$obs, sd$unit_covs, sd$event_covs, 
    type = "augmented", n_aug = 1)
  ps <- get_positions(fd)
  expect_true(
    all.equal(
      new_array(sd$event_covs$ec1, fd$data$ec1[ps]),
      sd$event_covs$ec1,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  expect_true(
    all.equal(
      fd$data$ec1[ps2],
      rep(sd$event_covs$ec1[,1], dim(sd$obs)[3]+1),
      check.attributes = FALSE
    )
  )
  
  # augmented with missingness
  sd <- simulate_flocker_data(augmented = TRUE, ragged_rep = TRUE)
  fd <- make_flocker_data(
    sd$obs, sd$unit_covs, sd$event_covs, 
    type = "augmented", n_aug = 1)
  ps <- get_positions(fd)
  expect_true(
    all.equal(
      new_array(sd$event_covs$ec1, fd$data$ec1[ps]),
      sd$event_covs$ec1,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  expect_true(
    all.equal(
      fd$data$ec1[ps2],
      rep(sd$event_covs$ec1[,1], dim(sd$obs)[3]+1),
      check.attributes = FALSE
    )
  )
  
  
  # multiseason
  sd <- simulate_flocker_data(
    n_pt = 10,
    n_sp = 1,
    n_season = 8,
    multiseason = "colex", 
    multi_init = "explicit"
  )
  fd <- make_flocker_data(
    sd$obs, sd$unit_covs, sd$event_covs, 
    type = "multi")
  ps <- get_positions(fd)
  expect_true(
    all.equal(
      new_array(sd$event_covs$ec1, fd$data$ec1[ps]),
      sd$event_covs$ec1,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  unit_covs_all <- sd$unit_covs[[1]]$uc1
  for(i in 2:8){
    unit_covs_all <- c(unit_covs_all, sd$unit_covs[[i]]$uc1)
  }
  expect_true(
    all.equal(
      fd$data$uc1[ps2],
      unit_covs_all,
      check.attributes = FALSE
    )
  )  
  
  
  # multiseason with missingness
  sd <- simulate_flocker_data(
    n_pt = 10,
    n_sp = 1,
    n_season = 8,
    multiseason = "colex", 
    multi_init = "explicit",
    ragged_rep = TRUE,
    missing_seasons = TRUE
  )
  suppressWarnings({
    fd <- make_flocker_data(
      sd$obs, sd$unit_covs, sd$event_covs, 
      type = "multi")
  })
  ps <- get_positions(fd)
  expect_true(
    all.equal(
      new_array(sd$event_covs$ec1, fd$data$ec1[ps]),
      sd$event_covs$ec1,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  # can't use unit covs here because they aren't NA in all relevant locations
  temp <- new_array(sd$obs[,1,], fd$data$ff_y[ps2])
  temp[temp == -99] <- NA
  expect_true(
    all.equal(
      temp,
      sd$obs[,1,],
      check.attributes = FALSE
    )
  )
  
  # multiseason with missingness for the whole first season
  sd <- simulate_flocker_data(
    n_pt = 10,
    n_sp = 1,
    n_season = 8,
    multiseason = "colex", 
    multi_init = "explicit",
    ragged_rep = TRUE,
    missing_seasons = TRUE
  )
  sd$obs[,,1] <- NA
  
  suppressWarnings({
    fd <- make_flocker_data(
      sd$obs, sd$unit_covs, sd$event_covs, 
      type = "multi")
  })

  ps <- get_positions(fd)
  expect_true(
    all.equal(
      new_array(sd$obs, fd$data$ff_y[ps]),
      sd$obs,
      check.attributes = FALSE
    )
  )
  ps2 <- get_positions(fd, unit_level = TRUE)
  # can't use unit covs here because they aren't NA in all relevant locations
  temp <- new_array(sd$obs[,1,], fd$data$ff_y[ps2])
  temp[temp == -99] <- NA
  expect_true(
    all.equal(
      temp,
      sd$obs[,1,],
      check.attributes = FALSE
    )
  )
})


test_that("emission_likelihood function returns expected output", {
  # Test cases for state 0
  obs1 <- matrix(c(0, 0, 0, 0, NA), nrow = 1)
  det1 <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7), nrow = 1)
  expected_output1 <- 1
  
  obs2 <- rbind(obs1, c(0, 1, 0, 1, 1))
  det2 <- rbind(det1, c(0.3, 0.4, 0.5, 0.6, 0.7))
  expected_output2 <- c(1,0)
  
  # Test cases for state 1
  obs3 <- matrix(c(1, 1, 0, 1, NA), nrow = 1)
  det3 <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7), nrow = 1)
  expected_output3 <- 0.3 * 0.4 * (1 - 0.5) * 0.6
  
  obs4 <- matrix(c(0, 0, 1, 1, 0), nrow = 1)
  det4 <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7), nrow = 1)
  expected_output4 <- (1 - 0.3) * (1 - 0.4) * 0.5 * 0.6 * (1 - 0.7)
  
  # Test the function with the test cases
  result1 <- emission_likelihood(0, obs1, det1)
  expect_identical(result1, expected_output1)
  
  result2 <- emission_likelihood(0, obs2, det2)
  expect_identical(result2, expected_output2)
  
  result3 <- emission_likelihood(1, obs3, det3)
  expect_identical(result3, expected_output3)
  
  result4 <- emission_likelihood(1, obs4, det4)
  expect_equal(result4, expected_output4)
  
  # Test the function with invalid inputs
  obs_invalid1 <- matrix(c(0, 1, -1, 1, 0), nrow = 1)
  obs_invalid2 <- matrix(c(0, 1, 2, 1, 0), nrow = 1)
  det_invalid1 <- matrix(c(0.5, 0.5, -0.5, 0.5, 0.5), nrow = 1)
  det_invalid2 <- matrix(c(0.5, 0.5, 1.5, 0.5, 0.5), nrow = 1)
  det_invalid3 <- matrix(c(NA, .5, .5, .5, .5), nrow = 1)
  
  expect_error(emission_likelihood(0, obs_invalid1, det1), "all\\(obs")
  expect_error(emission_likelihood(0, obs_invalid2, det1), "all\\(obs")
  expect_error(emission_likelihood(0, obs1, det_invalid1), "all\\(det")
  expect_error(emission_likelihood(0, obs1, det_invalid2), "all\\(det")
  expect_error(emission_likelihood(0, obs1, det_invalid3))
  expect_error(emission_likelihood(1, obs1, det_invalid3))
})


test_that("Z_from_emission returns correct values for valid inputs", {
  el0 <- c(0.1, 0.2, 0.3)
  el1 <- c(0.7, 0.8, 0.9)
  psi_unconditional <- c(0.4, 0.5, 0.6)
  expected_output <- psi_unconditional * el1 / 
    (psi_unconditional * el1 + (1 - psi_unconditional) * el0)
  expect_equal(Z_from_emission(el0, el1, psi_unconditional), expected_output)
})

test_that("Z_from_emission handles zeros and ones correctly", {
  el0 <- c(1, 0, 0)
  el1 <- c(0, 1, 1)
  psi_unconditional <- c(0, 1, 0.5)
  expected_output <- c(0, 1, 1) # When el0 is 0 and el1 is 1, the output should be 1
  expect_equal(Z_from_emission(el0, el1, psi_unconditional), expected_output)
})

test_that("Z_from_emission raises error with vectors of different lengths", {
  el0 <- c(0.1, 0.2)
  el1 <- c(0.7, 0.8, 0.9)
  psi_unconditional <- c(0.4, 0.5)
  expect_error(Z_from_emission(el0, el1, psi_unconditional))
})

test_that("Z_from_emission raises error with NA values in inputs", {
  el0 <- c(0.1, NA, 0.3)
  el1 <- c(0.7, 0.8, 0.9)
  psi_unconditional <- c(0.4, 0.5, 0.6)
  expect_error(Z_from_emission(el0, el1, psi_unconditional))
})

test_that("Z_from_emission raises error with negative values in inputs", {
  el0 <- c(0.1, -0.2, 0.3)
  el1 <- c(0.7, 0.8, 0.9)
  psi_unconditional <- c(0.4, 0.5, 0.6)
  expect_error(Z_from_emission(el0, el1, psi_unconditional))
})

test_that("Z_from_emission raises error with values greater than one in inputs", {
  el0 <- c(0.1, 0.2, 0.3)
  el1 <- c(1.1, 0.8, 0.9) # 1.1 is greater than 1
  psi_unconditional <- c(0.4, 0.5, 2.0) # 2.0 is greater than 1
  expect_error(Z_from_emission(el0, el1, psi_unconditional))
})

test_that("Z_from_emission handles equal emission likelihoods correctly", {
  el0 <- c(0.5, 0.5, 0.5)
  el1 <- c(0.5, 0.5, 0.5)
  psi_unconditional <- c(0.4, 0.5, 0.6)
  expected_output <- psi_unconditional / (psi_unconditional + (1 - psi_unconditional)) # Simplified formula for equal el0 and el1
  expect_equal(Z_from_emission(el0, el1, psi_unconditional), expected_output)
})

test_that("Z_from_emission handles scalar inputs correctly", {
  el0 <- 0.2
  el1 <- 0.8
  psi_unconditional <- 0.5
  expected_output <- 0.5 * 0.8 / (0.5 * 0.8 + (1 - 0.5) * 0.2)
  expect_equal(Z_from_emission(el0, el1, psi_unconditional), expected_output)
})

sd <- simulate_flocker_data()
fd_single <- make_flocker_data(sd$obs, sd$unit_covs, sd$event_covs)
fd_single_C <- make_flocker_data(sd$obs, sd$unit_covs)

sd <- simulate_flocker_data(augmented = TRUE)
fd_augmented <- make_flocker_data(sd$obs, sd$unit_covs, sd$event_covs, type = "augmented", n_aug = 10)

sd <- simulate_flocker_data(n_season = 3, multiseason = "colex", multi_init = "explicit")
fd_multi <- make_flocker_data(sd$obs, sd$unit_covs, sd$event_covs, type = "multi")

test_that("validate_flock_params works as expected", {
  f_occ <- ~ uc1
  f_det <- ~ uc1 + ec1
  f_col <- NULL
  f_ex <- NULL
  f_auto <- NULL
  flocker_data <- fd_single
  multiseason <- NULL
  multi_init <- NULL
  augmented <- FALSE
  threads <- NULL
  
  expect_silent(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads))
  
  
  f_occ <- ~ uc1 + ec1
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads))
  
  f_occ <- y ~ uc1
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads))
  
  flocker_data <- fd_single_C
  f_occ <- ~ uc1
  f_det <- ~ uc1
  
  expect_silent(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  
  
  flocker_data <- fd_augmented
  f_det <- ~ uc1 + ec1
  augmented <- TRUE
  
  expect_silent(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads))
  
  
  flocker_data <- fd_multi
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads))
  augmented <- FALSE
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  multiseason <- "colex"
  multi_init <- "explicit"
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_col <- ~ uc1
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_ex <- ~ uc1
  
  expect_silent(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                      f_col, f_ex, multi_init, f_auto, augmented, threads))
  
  multiseason <- "autologistic"
  multi_init <- "equilibrium"
  
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_auto <- ~ uc1
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_auto <- NULL
  f_occ <- NULL
  f_ex <- NULL
  expect_error(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
  f_auto <- ~ uc1
  expect_silent(validate_flock_params(f_occ, f_det, flocker_data, multiseason, 
                                     f_col, f_ex, multi_init, f_auto, augmented, threads))
})

test_that("formula_error works", {
  result <- formula_error("x")
  expect_is(result, "character")
  expect_identical(result, "Formula error: x formula has incorrect syntax.")
})

test_that("is_formula and is_flocker_formula work", {
  expect_true(is_formula(stats::formula(a ~ b)))
  expect_true(is_formula(stats::formula(~ b)))
  expect_true(is_formula(stats::formula(~1)))
  expect_true(is_formula(stats::formula(~ a ~ b)))
  
  expect_false(is_formula(1))
  expect_false(is_formula(list(stats::formula(~1), stats::formula(a ~ b))))
  
  expect_false(is_flocker_formula(stats::formula(a ~ b)))
  expect_false(is_flocker_formula(stats::formula(~ a ~ b)))
  expect_true(is_flocker_formula(stats::formula(~ b)))
  expect_true(is_flocker_formula(stats::formula(~ b + (1 || c))))
})

test_that("is_flocker_data works", {
  sfd <- simulate_flocker_data()
  expect_false(is_flocker_data(sfd))
  
  expect_true(
    make_flocker_data(sfd$obs, sfd$unit_covs) |>
      is_flocker_data()
  )
})  


test_that("is_named_list returns TRUE for a named list with unique names", {
  named_list <- list(a = 1, b = 2, c = 3)
  expect_true(is_named_list(named_list))
})

test_that("is_named_list returns FALSE for an unnamed list", {
  unnamed_list <- list(1, 2, 3)
  expect_false(is_named_list(unnamed_list))
})

test_that("is_named_list returns FALSE for a partially named list", {
  partially_named_list <- list(a = 1, 2, c = 3)
  expect_false(is_named_list(partially_named_list))
})

test_that("is_named_list returns FALSE for a named list with duplicate names", {
  duplicate_named_list <- list(a = 1, b = 2, a = 3)
  expect_false(is_named_list(duplicate_named_list))
})

test_that("is_named_list returns FALSE for an empty list", {
  empty_list <- list()
  expect_false(is_named_list(empty_list))
})

test_that("is_named_list returns FALSE for a non-list object", {
  non_list_object <- 42
  expect_false(is_named_list(non_list_object))
})

test_that("is_one_logical returns TRUE for a single logical value (TRUE)", {
  single_logical_true <- TRUE
  expect_true(is_one_logical(single_logical_true))
})

test_that("is_one_logical returns TRUE for a single logical value (FALSE)", {
  single_logical_false <- FALSE
  expect_true(is_one_logical(single_logical_false))
})

test_that("is_one_logical returns FALSE for a numeric value", {
  numeric_value <- 42
  expect_false(is_one_logical(numeric_value))
})

test_that("is_one_logical returns FALSE for a character value", {
  character_value <- "test"
  expect_false(is_one_logical(character_value))
})

test_that("is_one_logical returns FALSE for a vector of logical values", {
  logical_vector <- c(TRUE, FALSE)
  expect_false(is_one_logical(logical_vector))
})

test_that("is_one_logical returns FALSE for a NULL value", {
  null_value <- NULL
  expect_false(is_one_logical(null_value))
})

test_that("is_one_pos_int returns TRUE for a single positive integer greater than m", {
  single_positive_int <- 5
  m <- 3
  expect_true(is_one_pos_int(single_positive_int, m))
})

test_that("is_one_pos_int returns FALSE for a single positive integer equal to m", {
  single_positive_int <- 5
  m <- 5
  expect_false(is_one_pos_int(single_positive_int, m))
})

test_that("is_one_pos_int returns FALSE for a single positive integer less than m", {
  single_positive_int <- 3
  m <- 5
  expect_false(is_one_pos_int(single_positive_int, m))
})

test_that("is_one_pos_int returns FALSE for a single negative integer", {
  negative_int <- -5
  m <- 0
  expect_false(is_one_pos_int(negative_int, m))
})

test_that("is_one_pos_int returns FALSE for a single non-integer numeric value", {
  non_integer_numeric <- 2.5
  m <- 0
  expect_false(is_one_pos_int(non_integer_numeric, m))
})

test_that("is_one_pos_int returns FALSE for a character value", {
  character_value <- "test"
  m <- 0
  expect_false(is_one_pos_int(character_value, m))
})

test_that("is_one_pos_int returns FALSE for a vector of integers", {
  integer_vector <- c(2, 5, 8)
  m <- 0
  expect_false(is_one_pos_int(integer_vector, m))
})

test_that("is_one_pos_int returns FALSE for a NULL value", {
  null_value <- NULL
  m <- 0
  expect_false(is_one_pos_int(null_value, m))
})

test_that("shared_elements returns correct shared elements for two non-empty vectors", {
  vec1 <- c(1, 2, 3, 4, 5, 6)
  vec2 <- c(4, 5, 6, 7, 8, 9)
  expected_output <- c(4, 5, 6)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements returns an empty vector when there are no shared elements", {
  vec1 <- c(1, 2, 3)
  vec2 <- c(4, 5, 6)
  expected_output <- integer(0)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements handles vectors with duplicate elements correctly", {
  vec1 <- c(1, 1, 2, 2, 3, 3)
  vec2 <- c(2, 2, 3, 3, 4, 4)
  expected_output <- c(2, 3)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements returns an empty vector when both input vectors are empty", {
  vec1 <- integer(0)
  vec2 <- integer(0)
  expected_output <- integer(0)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements returns an empty vector when one input vector is empty", {
  vec1 <- c(1, 2, 3)
  vec2 <- integer(0)
  expected_output <- integer(0)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements works with character vectors", {
  vec1 <- c("apple", "banana", "cherry")
  vec2 <- c("banana", "cherry", "date")
  expected_output <- c("banana", "cherry")
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements works with mixed-type vectors", {
  vec1 <- c(1, "apple", 2, "banana")
  vec2 <- c("banana", 2, "date", 3)
  expected_output <- c(2, "banana")
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("max_position_not_na works correctly", {
  # Test with no NAs and no -99 values
  expect_equal(max_position_not_na(c(1, 2, 3, 4, 5)), 5)
  
  # Test with NAs only
  expect_equal(max_position_not_na(c(NA, NA, NA)), 0)
  
  # Test with -99 values treated as NAs
  expect_equal(max_position_not_na(c(1, 2, -99, 4, 5), treat_m99_NA = TRUE), 5)
  expect_equal(max_position_not_na(c(1, 2, -99, 4, -99), treat_m99_NA = TRUE), 4)
  
  # Test with NAs and -99 values at the end of the vector
  expect_equal(max_position_not_na(c(1, 2, 3, NA, -99), treat_m99_NA = TRUE), 3)
  expect_equal(max_position_not_na(c(1, 2, 3, NA, -99), treat_m99_NA = FALSE), 5)
  
  # Test with NAs and -99 values at the beginning of the vector
  expect_equal(max_position_not_na(c(NA, -99, 1, 2, 3), treat_m99_NA = TRUE), 5)
  expect_equal(max_position_not_na(c(NA, -99, 1, 2, 3), treat_m99_NA = FALSE), 5)
  
  # Test with non-numeric vector
  expect_equal(max_position_not_na(c("A", "B", "C", "D", "E")), 5)
  expect_equal(max_position_not_na(c("A", "B", "C", "D", "NA")), 5)
  expect_equal(max_position_not_na(c("A", "B", "C", "D", "-99"), treat_m99_NA = TRUE), 4)
})


test_that("remove_rownames works correctly", {
  # Test with a matrix
  m1 <- matrix(1:9, nrow = 3, ncol = 3, dimnames = list(c("r1", "r2", "r3"), c("c1", "c2", "c3")))
  m1_expected <- matrix(1:9, nrow = 3, ncol = 3, dimnames = list(c(), c("c1", "c2", "c3")))
  expect_equal(remove_rownames(m1), m1_expected)
  
  # Test with a data.frame
  df1 <- data.frame(a = 1:3, b = 4:6, row.names = c("r1", "r2", "r3"))
  df1_expected <- data.frame(a = 1:3, b = 4:6)
  expect_equal(remove_rownames(df1), df1_expected)
  
  # Test with a tibble
  if (requireNamespace("tibble", quietly = TRUE)) {
    tb1 <- tibble::tibble(a = 1:3, b = 4:6)
    suppressWarnings(rownames(tb1) <- c("r1", "r2", "r3"))
    tb1_expected <- tibble::tibble(a = 1:3, b = 4:6)
    expect_equal(remove_rownames(tb1), tb1_expected)
  }
})

test_that("rbinom2 works correctly", {
  r1 <- withr::with_seed(seed = 1, code = stats::rbinom(10, rep(c(1,2), 5), runif(10)))
  r2 <- withr::with_seed(seed = 1, code = rbinom2(10, rep(c(1,2), 5), runif(10)))
  expect_identical(r2, r1)
  r3 <- withr::with_seed(seed = 1, code = rbinom2(11, c(rep(c(1,2), 5), 1), c(runif(10), NA)))
  expect_identical(r3[1:10], r1)
  expect_true(is.na(r3[11]))
})
