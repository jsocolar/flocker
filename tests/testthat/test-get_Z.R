test_that("get_Z gives valid returns", {
  Z_single <- get_Z(example_flocker_model_single)
  expect_true(all(Z_single >= 0))
  expect_true(all(Z_single <= 1))
  expect_true(any(Z_single == 1))
  
  Z_single_C <- get_Z(example_flocker_model_single_C)
  expect_true(all(Z_single_C >= 0))
  expect_true(all(Z_single_C <= 1))
  expect_true(any(Z_single_C == 1))
  
  Z_augmented <- get_Z(example_flocker_model_aug)
  expect_true(all(Z_augmented >= 0))
  expect_true(all(Z_augmented <= 1))
  expect_true(any(Z_augmented == 1))
  
  Z_multi_colex_ex <- get_Z(example_flocker_model_multi_colex_ex)
  expect_true(all(Z_multi_colex_ex >= 0, na.rm = T))
  expect_true(all(Z_multi_colex_ex <= 1, na.rm = T))
  expect_true(any(Z_multi_colex_ex == 1, na.rm = T))
  
  Z_multi_colex_eq <- get_Z(example_flocker_model_multi_colex_eq)
  expect_true(all(Z_multi_colex_eq >= 0, na.rm = T))
  expect_true(all(Z_multi_colex_eq <= 1, na.rm = T))
  expect_true(any(Z_multi_colex_eq == 1, na.rm = T))
  
  Z_multi_auto_ex <- get_Z(example_flocker_model_multi_auto_ex)
  expect_true(all(Z_multi_auto_ex >= 0, na.rm = T))
  expect_true(all(Z_multi_auto_ex <= 1, na.rm = T))
  expect_true(any(Z_multi_auto_ex == 1, na.rm = T))
  
  Z_multi_auto_eq <- get_Z(example_flocker_model_multi_auto_eq)
  expect_true(all(Z_multi_auto_eq >= 0, na.rm = T))
  expect_true(all(Z_multi_auto_eq <= 1, na.rm = T))
  expect_true(any(Z_multi_auto_eq == 1, na.rm = T))
  
  
  
  Z_single_nohist <- get_Z(example_flocker_model_single, history_condition = FALSE)
  expect_true(all(Z_single_nohist >= 0))
  expect_true(all(Z_single_nohist <= 1))
  expect_true(sum(Z_single == 1) > sum(Z_single_nohist == 1))
  
  Z_single_C_nohist <- get_Z(example_flocker_model_single_C, history_condition = FALSE)
  expect_true(all(Z_single_C_nohist >= 0))
  expect_true(all(Z_single_C_nohist <= 1))
  expect_true(sum(Z_single_C == 1) > sum(Z_single_C_nohist == 1))
  
  Z_augmented_nohist <- get_Z(example_flocker_model_aug, history_condition = FALSE)
  expect_true(all(Z_augmented_nohist >= 0))
  expect_true(all(Z_augmented_nohist <= 1))
  expect_true(sum(Z_augmented == 1) > sum(Z_augmented_nohist == 1))
  
  Z_multi_colex_ex_nohist <- get_Z(example_flocker_model_multi_colex_ex, history_condition = FALSE)
  expect_true(all(Z_multi_colex_ex_nohist >= 0, na.rm = T))
  expect_true(all(Z_multi_colex_ex_nohist <= 1, na.rm = T))
  expect_true(sum(Z_multi_colex_ex == 1, na.rm = TRUE) > sum(Z_multi_colex_ex_nohist == 1, na.rm = TRUE))
  
  Z_multi_colex_eq_nohist <- get_Z(example_flocker_model_multi_colex_eq, history_condition = FALSE)
  expect_true(all(Z_multi_colex_eq_nohist >= 0, na.rm = T))
  expect_true(all(Z_multi_colex_eq_nohist <= 1, na.rm = T))
  expect_true(sum(Z_multi_colex_eq == 1, na.rm = TRUE) > sum(Z_multi_colex_eq_nohist == 1, na.rm = TRUE))
  
  Z_multi_auto_ex_nohist <- get_Z(example_flocker_model_multi_auto_ex, history_condition = FALSE)
  expect_true(all(Z_multi_auto_ex_nohist >= 0, na.rm = T))
  expect_true(all(Z_multi_auto_ex_nohist <= 1, na.rm = T))
  expect_true(sum(Z_multi_auto_ex == 1, na.rm = TRUE) > sum(Z_multi_auto_ex_nohist == 1, na.rm = TRUE))
  
  Z_multi_auto_eq_nohist <- get_Z(example_flocker_model_multi_auto_eq, history_condition = FALSE)
  expect_true(all(Z_multi_auto_eq_nohist >= 0, na.rm = T))
  expect_true(all(Z_multi_auto_eq_nohist <= 1, na.rm = T))
  expect_true(sum(Z_multi_auto_eq == 1, na.rm = TRUE) > sum(Z_multi_auto_eq_nohist == 1, na.rm = TRUE))
  
})


test_that("get_Z with sampling gives valid returns", {
  Z_single <- get_Z(example_flocker_model_single, sample = TRUE)
  expect_true(all(Z_single %in% c(0, 1, NA)))
  
  Z_single_C <- get_Z(example_flocker_model_single_C, sample = TRUE)
  expect_true(all(Z_single_C %in% c(0, 1, NA)))
  
  Z_augmented <- get_Z(example_flocker_model_aug, sample = TRUE)
  expect_true(all(Z_augmented %in% c(0, 1, NA)))
  
  Z_multi_colex_ex <- get_Z(example_flocker_model_multi_colex_ex, sample = TRUE)
  expect_true(all(Z_multi_colex_ex %in% c(0, 1, NA)))
  
  Z_multi_colex_eq <- get_Z(example_flocker_model_multi_colex_eq, sample = TRUE)
  expect_true(all(Z_multi_colex_eq %in% c(0, 1, NA)))
  
  Z_multi_auto_ex <- get_Z(example_flocker_model_multi_auto_ex, sample = TRUE)
  expect_true(all(Z_multi_auto_ex %in% c(0, 1, NA)))
  
  Z_multi_auto_eq <- get_Z(example_flocker_model_multi_auto_eq, sample = TRUE)
  expect_true(all(Z_multi_auto_eq %in% c(0, 1, NA)))

})

test_that("new_data works as expected", {
  fd1 <- simulate_flocker_data(n_sp = 5)
  mfd1 <- make_flocker_data(fd1$obs, fd1$unit_covs, fd1$event_covs)
  expect_silent(get_Z(example_flocker_model_single, new_data = mfd1))
  expect_silent(get_Z(example_flocker_model_single, history_condition = FALSE, new_data = fd1$unit_covs))
  expect_error(get_Z(example_flocker_model_single, history_condition = TRUE, new_data = fd1$unit_covs))
})


test_that("incorrect types cause an error", {
  expect_error(forward_sim("init", c(0.5, 0.5), c(0.5, 0.5)))
  expect_error(forward_sim(0.5, "colo", c(0.5, 0.5)))
  expect_error(forward_sim(0.5, c(0.5, 0.5), "ex"))
})

test_that("negative probabilities cause an error", {
  expect_silent(forward_sim(.4, c(.5, .5), c(.5, .5)))
  expect_error(forward_sim(.4, c(.5, NA), c(.5, .5)))
  expect_error(forward_sim(-0.1, c(0.5, .5), c(0.5, .5)))
  expect_error(forward_sim(0.5, c(.5, -0.5), c(.5, 0.5)))
  expect_error(forward_sim(0.5, c(NA, 0.5), c(NA, -0.5)))
})

test_that("correct inputs produce correct outputs", {
  # Simple case with no NA values and no sampling
  expect_equal(forward_sim(0.5, c(0.5, 0.5), c(0.5, 0.5)), c(0.5, 0.5))
  # Simple case with sampling, the result is stochastic, so we can't test exact values, just types and lengths
  result <- forward_sim(0.5, c(0.5, 0.5), c(0.5, 0.5), sample = TRUE)
  expect_type(result, "integer")
  expect_length(result, 2)
})

test_that("edge cases are handled", {
  # Case with all NA values should return NA
  expect_error(forward_sim(NA, c(NA, NA), c(NA, NA)))
  # Case with empty vectors
  expect_error(forward_sim(0.5, numeric(0), numeric(0)))
})

# Test for forward_backward_algorithm function
test_that("forward_backward_algorithm returns correct results", {
  el0 <- c(0.1, 0.2, 0.3)
  el1 <- c(0.9, 0.8, 0.7)
  init <- 0.5
  colo <- c(0.6, 0.7, 0.8)
  ex <- c(0.4, 0.3, 0.2)
  result <- forward_backward_algorithm(el0, el1, init, colo, ex)
  expect_type(result, "double")
  expect_length(result, length(el0))
})

# Test for forward_backward_sampling function
test_that("forward_backward_sampling returns correct results", {
  el0 <- c(0.1, 0.2, 0.3)
  el1 <- c(0.9, 0.8, 0.7)
  init <- 0.5
  colo <- c(0.6, 0.7, 0.8)
  ex <- c(0.4, 0.3, 0.2)
  result <- forward_backward_sampling(el0, el1, init, colo, ex)
  expect_type(result, "double")
  expect_length(result, length(el0))
})

# Test for forward_algorithm function
test_that("forward_algorithm returns correct results", {
  el0 <- c(0.1, 0.2, 0.3)
  el1 <- c(0.9, 0.8, 0.7)
  init <- 0.5
  colo <- c(0.6, 0.7, 0.8)
  ex <- c(0.4, 0.3, 0.2)
  result <- forward_algorithm(el0, el1, init, colo, ex)
  expect_type(result, "double")
  expect_equal(dim(result), c(length(el0), 2))
})

# Test for backward_algorithm function
test_that("backward_algorithm returns correct results", {
  el0 <- c(0.1, 0.2, 0.3)
  el1 <- c(0.9, 0.8, 0.7)
  colo <- c(0.6, 0.7, 0.8)
  ex <- c(0.4, 0.3, 0.2)
  result <- backward_algorithm(el0, el1, colo, ex)
  expect_type(result, "double")
  expect_equal(dim(result), c(length(el0), 2))
})

