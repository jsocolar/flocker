test_that("make_occupancy_single_lpmf works correctly", {
  # Test with max_rep = 2
  stan_code_2 <- make_occupancy_single_lpmf(2)
  expect_match(stan_code_2, "array[] int vint4", fixed = TRUE)
  expect_match(stan_code_2, "array[] int vint5", fixed = TRUE)
  expect_no_match(stan_code_2, "array[] int vint6", fixed = TRUE)
  
  expect_match(stan_code_2, "index_array[,1] = vint4", fixed = TRUE)
  expect_match(stan_code_2, "index_array[,2] = vint5", fixed = TRUE)
  expect_no_match(stan_code_2, "index_array[,3", fixed = TRUE)
  
  # Test with max_rep = 4
  stan_code_4 <- make_occupancy_single_lpmf(4)
  expect_match(stan_code_4, "array[] int vint4", fixed = TRUE)
  expect_match(stan_code_4, "array[] int vint5", fixed = TRUE)
  expect_match(stan_code_4, "array[] int vint6", fixed = TRUE)
  expect_match(stan_code_4, "array[] int vint7", fixed = TRUE)
  expect_no_match(stan_code_4, "array[] int vint8", fixed = TRUE)
  
  expect_match(stan_code_4, "index_array[,1] = vint4", fixed = TRUE)
  expect_match(stan_code_4, "index_array[,2] = vint5", fixed = TRUE)
  expect_match(stan_code_4, "index_array[,3] = vint6", fixed = TRUE)
  expect_match(stan_code_4, "index_array[,4] = vint7", fixed = TRUE)
  expect_no_match(stan_code_4, "index_array[,5", fixed = TRUE)
  
  # Test with max_rep = 1 (invalid input)
  expect_error(make_occupancy_single_lpmf(1), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_single_lpmf(-2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_single_lpmf(3.5), "max_rep must be an integer greater than 1")
})


test_that("make_occupancy_augmented_lpmf works correctly", {
  # Test with max_rep = 2
  stan_code_2 <- make_occupancy_augmented_lpmf(2)
  expect_match(stan_code_2, "array[] int vint7", fixed = TRUE)
  expect_match(stan_code_2, "array[] int vint8", fixed = TRUE)
  expect_no_match(stan_code_2, "array[] int vint9", fixed = TRUE)
  
  expect_match(stan_code_2, "index_array[,1] = vint7", fixed = TRUE)
  expect_match(stan_code_2, "index_array[,2] = vint8", fixed = TRUE)
  expect_no_match(stan_code_2, "index_array[,3", fixed = TRUE)
  
  # Test with max_rep = 4
  stan_code_4 <- make_occupancy_augmented_lpmf(4)
  expect_match(stan_code_4, "array[] int vint7", fixed = TRUE)
  expect_match(stan_code_4, "array[] int vint8", fixed = TRUE)
  expect_match(stan_code_4, "array[] int vint9", fixed = TRUE)
  expect_match(stan_code_4, "array[] int vint10", fixed = TRUE)
  expect_no_match(stan_code_4, "array[] int vint11", fixed = TRUE)
  
  expect_match(stan_code_4, "index_array[,1] = vint7", fixed = TRUE)
  expect_match(stan_code_4, "index_array[,2] = vint8", fixed = TRUE)
  expect_match(stan_code_4, "index_array[,3] = vint9", fixed = TRUE)
  expect_match(stan_code_4, "index_array[,4] = vint10", fixed = TRUE)
  expect_no_match(stan_code_4, "index_array[,5", fixed = TRUE)
  
  # Test with max_rep = 1 (invalid input)
  expect_error(make_occupancy_augmented_lpmf(1), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_augmented_lpmf(-2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_augmented_lpmf(3.5), "max_rep must be an integer greater than 1")
})


test_that("no-argument functions have correct return types", {
  expect_is(make_occupancy_single_C_lpmf(), "character")
  expect_is(make_emission_1(), "character")
  expect_is(make_emission_1_C(), "character")
  expect_is(make_emission_1_fp(), "character")
  expect_is(make_emission_0_fp(), "character")
  expect_is(make_colex_likelihoods(), "character")
  expect_is(make_colex_fp_likelihoods(), "character")
  expect_is(make_forward_colex(), "character")
  expect_is(make_forward_colex_fp(), "character")
})
