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
  
  # Test with max_rep = "foo" (invalid input)
  expect_error(make_occupancy_single_lpmf("foo"), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = c(2, 3) (invalid input)
  expect_error(make_occupancy_single_lpmf(c(2, 3)), "max_rep must be an integer greater than 1")
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
  
  # Test with max_rep = "foo" (invalid input)
  expect_error(make_occupancy_augmented_lpmf("foo"), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = c(2, 3) (invalid input)
  expect_error(make_occupancy_augmented_lpmf(c(2, 3)), "max_rep must be an integer greater than 1")
  
})


test_that("no-argument functions have correct return types", {
  expect_is(make_occupancy_single_C_lpmf(), "character")
  expect_is(make_emission_1(), "character")
  expect_is(make_emission_1_C(), "character")
  expect_is(make_colex_likelihoods(), "character")
  expect_is(make_forward_colex(), "character")
})

test_that("make_occupancy_multi_colex_lpmf works correctly", {
  # Test with max_rep = 2, max_year = 3
  stan_code_23 <- make_occupancy_multi_colex_lpmf(2, 3)
  expect_match(stan_code_23, "array[] int vint6", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint7", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint8", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint9", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint10", fixed = TRUE)
  expect_no_match(stan_code_23, "array[] int vint11", fixed = TRUE)

  expect_match(stan_code_23, "unit_index_array[,1] = vint6", fixed = TRUE)
  expect_match(stan_code_23, "unit_index_array[,2] = vint7", fixed = TRUE)
  expect_match(stan_code_23, "unit_index_array[,3] = vint8", fixed = TRUE)
  expect_no_match(stan_code_23, "unit_index_array[,4]", fixed = TRUE)
  
  expect_match(stan_code_23, "array[vint2[1], 2] int visit_index_array;", fixed = TRUE)
  expect_no_match(stan_code_23, "array[vint2[1], 3] int visit_index_array;", fixed = TRUE)
  
  expect_match(stan_code_23, "visit_index_array[,1] = vint9", fixed = TRUE)
  expect_match(stan_code_23, "visit_index_array[,2] = vint10", fixed = TRUE)
  expect_no_match(stan_code_23, "visit_index_array[,3] = vint11", fixed = TRUE)
  expect_no_match(stan_code_23, "visit_index_array[,1] = vint10", fixed = TRUE)
  
  # Test with max_rep = 1 (invalid input)
  expect_error(make_occupancy_multi_colex_lpmf(1, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = "foo" (invalid input)
  expect_error(make_occupancy_multi_colex_lpmf("foo", 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = c(2, 3) (invalid input)
  expect_error(make_occupancy_multi_colex_lpmf(c(2, 3), 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_multi_colex_lpmf(-2, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_multi_colex_lpmf(3.5, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_year = 1 (invalid input)
  expect_error(make_occupancy_multi_colex_lpmf(2, 1), "max_year must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_multi_colex_lpmf(2, -2), "max_year must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_multi_colex_lpmf(2, 3.5), "max_year must be an integer greater than 1")
})

test_that("make_occupancy_multi_colex_eq_lpmf works correctly", {
  # Test with max_rep = 2, max_year = 3
  stan_code_23 <- make_occupancy_multi_colex_eq_lpmf(2, 3)
  expect_match(stan_code_23, "array[] int vint6", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint7", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint8", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint9", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint10", fixed = TRUE)
  expect_no_match(stan_code_23, "array[] int vint11", fixed = TRUE)
  
  expect_match(stan_code_23, "unit_index_array[,1] = vint6", fixed = TRUE)
  expect_match(stan_code_23, "unit_index_array[,2] = vint7", fixed = TRUE)
  expect_match(stan_code_23, "unit_index_array[,3] = vint8", fixed = TRUE)
  expect_no_match(stan_code_23, "unit_index_array[,4]", fixed = TRUE)
  
  expect_match(stan_code_23, "array[vint2[1], 2] int visit_index_array;", fixed = TRUE)
  expect_no_match(stan_code_23, "array[vint2[1], 3] int visit_index_array;", fixed = TRUE)
  
  expect_match(stan_code_23, "visit_index_array[,1] = vint9", fixed = TRUE)
  expect_match(stan_code_23, "visit_index_array[,2] = vint10", fixed = TRUE)
  expect_no_match(stan_code_23, "visit_index_array[,3] = vint11", fixed = TRUE)
  expect_no_match(stan_code_23, "visit_index_array[,1] = vint10", fixed = TRUE)
  
  # Test with max_rep = 1 (invalid input)
  expect_error(make_occupancy_multi_colex_eq_lpmf(1, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = "foo" (invalid input)
  expect_error(make_occupancy_multi_colex_eq_lpmf("foo", 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = c(2, 3) (invalid input)
  expect_error(make_occupancy_multi_colex_eq_lpmf(c(2, 3), 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_multi_colex_eq_lpmf(-2, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_multi_colex_eq_lpmf(3.5, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_year = 1 (invalid input)
  expect_error(make_occupancy_multi_colex_eq_lpmf(2, 1), "max_year must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_multi_colex_eq_lpmf(2, -2), "max_year must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_multi_colex_eq_lpmf(2, 3.5), "max_year must be an integer greater than 1")
})

test_that("make_occupancy_multi_autologistic_lpmf works correctly", {
  # Test with max_rep = 2, max_year = 3
  stan_code_23 <- make_occupancy_multi_autologistic_lpmf(2, 3)
  expect_match(stan_code_23, "array[] int vint6", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint7", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint8", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint9", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint10", fixed = TRUE)
  expect_no_match(stan_code_23, "array[] int vint11", fixed = TRUE)
  
  expect_match(stan_code_23, "unit_index_array[,1] = vint6", fixed = TRUE)
  expect_match(stan_code_23, "unit_index_array[,2] = vint7", fixed = TRUE)
  expect_match(stan_code_23, "unit_index_array[,3] = vint8", fixed = TRUE)
  expect_no_match(stan_code_23, "unit_index_array[,4]", fixed = TRUE)
  
  expect_match(stan_code_23, "array[vint2[1], 2] int visit_index_array;", fixed = TRUE)
  expect_no_match(stan_code_23, "array[vint2[1], 3] int visit_index_array;", fixed = TRUE)
  
  expect_match(stan_code_23, "visit_index_array[,1] = vint9", fixed = TRUE)
  expect_match(stan_code_23, "visit_index_array[,2] = vint10", fixed = TRUE)
  expect_no_match(stan_code_23, "visit_index_array[,3] = vint11", fixed = TRUE)
  expect_no_match(stan_code_23, "visit_index_array[,1] = vint10", fixed = TRUE)
  
  # Test with max_rep = 1 (invalid input)
  expect_error(make_occupancy_multi_autologistic_lpmf(1, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = "foo" (invalid input)
  expect_error(make_occupancy_multi_autologistic_lpmf("foo", 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = c(2, 3) (invalid input)
  expect_error(make_occupancy_multi_autologistic_lpmf(c(2, 3), 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_multi_autologistic_lpmf(-2, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_multi_autologistic_lpmf(3.5, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_year = 1 (invalid input)
  expect_error(make_occupancy_multi_autologistic_lpmf(2, 1), "max_year must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_multi_autologistic_lpmf(2, -2), "max_year must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_multi_autologistic_lpmf(2, 3.5), "max_year must be an integer greater than 1")
})

test_that("make_occupancy_multi_autologistic_eq_lpmf works correctly", {
  # Test with max_rep = 2, max_year = 3
  stan_code_23 <- make_occupancy_multi_autologistic_eq_lpmf(2, 3)
  expect_match(stan_code_23, "array[] int vint6", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint7", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint8", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint9", fixed = TRUE)
  expect_match(stan_code_23, "array[] int vint10", fixed = TRUE)
  expect_no_match(stan_code_23, "array[] int vint11", fixed = TRUE)
  
  expect_match(stan_code_23, "unit_index_array[,1] = vint6", fixed = TRUE)
  expect_match(stan_code_23, "unit_index_array[,2] = vint7", fixed = TRUE)
  expect_match(stan_code_23, "unit_index_array[,3] = vint8", fixed = TRUE)
  expect_no_match(stan_code_23, "unit_index_array[,4]", fixed = TRUE)
  
  expect_match(stan_code_23, "array[vint2[1], 2] int visit_index_array;", fixed = TRUE)
  expect_no_match(stan_code_23, "array[vint2[1], 3] int visit_index_array;", fixed = TRUE)
  
  expect_match(stan_code_23, "visit_index_array[,1] = vint9", fixed = TRUE)
  expect_match(stan_code_23, "visit_index_array[,2] = vint10", fixed = TRUE)
  expect_no_match(stan_code_23, "visit_index_array[,3] = vint11", fixed = TRUE)
  expect_no_match(stan_code_23, "visit_index_array[,1] = vint10", fixed = TRUE)
  
  # Test with max_rep = 1 (invalid input)
  expect_error(make_occupancy_multi_autologistic_eq_lpmf(1, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = "foo" (invalid input)
  expect_error(make_occupancy_multi_autologistic_eq_lpmf("foo", 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = c(2, 3) (invalid input)
  expect_error(make_occupancy_multi_autologistic_eq_lpmf(c(2, 3), 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_multi_autologistic_eq_lpmf(-2, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_multi_autologistic_eq_lpmf(3.5, 2), "max_rep must be an integer greater than 1")
  
  # Test with max_year = 1 (invalid input)
  expect_error(make_occupancy_multi_autologistic_eq_lpmf(2, 1), "max_year must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_multi_autologistic_eq_lpmf(2, -2), "max_year must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_multi_autologistic_eq_lpmf(2, 3.5), "max_year must be an integer greater than 1")
})

test_that("make_occupancy_single_threaded_lpmf works correctly", {
  # Test with max_rep = 2
  stan_code_21 <- make_occupancy_single_threaded_lpmf(2, 1)
  expect_match(stan_code_21, "array[] int vint4", fixed = TRUE)
  expect_match(stan_code_21, "array[] int vint5", fixed = TRUE)
  expect_no_match(stan_code_21, "array[] int vint6", fixed = TRUE)
  
  # Test with max_rep = 4
  stan_code_42 <- make_occupancy_single_threaded_lpmf(4, 2)
  expect_match(stan_code_42, "array[] int vint4", fixed = TRUE)
  expect_match(stan_code_42, "array[] int vint5", fixed = TRUE)
  expect_match(stan_code_42, "array[] int vint6", fixed = TRUE)
  expect_match(stan_code_42, "array[] int vint7", fixed = TRUE)
  expect_no_match(stan_code_42, "array[] int vint8", fixed = TRUE)
  
  # Test with max_rep = 1 (invalid input)
  expect_error(make_occupancy_single_threaded_lpmf(1, 1), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_single_threaded_lpmf(-2, 1), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_single_threaded_lpmf(3.5, 1), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = "foo" (invalid input)
  expect_error(make_occupancy_single_threaded_lpmf("foo", 1), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = c(2, 3) (invalid input)
  expect_error(make_occupancy_single_threaded_lpmf(c(2, 3), 1), "max_rep must be an integer greater than 1")
})

test_that("make_occupancy_single_partial_sum works correctly", {
  # Test with max_rep = 2
  stan_code_2 <- make_occupancy_single_partial_sum(2)
  expect_match(stan_code_2, "array[] int vint4", fixed = TRUE)
  expect_match(stan_code_2, "array[] int vint5", fixed = TRUE)
  expect_no_match(stan_code_2, "array[] int vint6", fixed = TRUE)
  
  expect_match(stan_code_2, "index_array[,1] = vint4", fixed = TRUE)
  expect_match(stan_code_2, "index_array[,2] = vint5", fixed = TRUE)
  expect_no_match(stan_code_2, "index_array[,3", fixed = TRUE)
  
  # Test with max_rep = 4
  stan_code_4 <- make_occupancy_single_partial_sum(4)
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
  expect_error(make_occupancy_single_partial_sum(1), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = -2 (invalid input)
  expect_error(make_occupancy_single_partial_sum(-2), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = 3.5 (invalid input)
  expect_error(make_occupancy_single_partial_sum(3.5), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = "foo" (invalid input)
  expect_error(make_occupancy_single_partial_sum("foo"), "max_rep must be an integer greater than 1")
  
  # Test with max_rep = c(2, 3) (invalid input)
  expect_error(make_occupancy_single_partial_sum(c(2, 3)), "max_rep must be an integer greater than 1")
})

