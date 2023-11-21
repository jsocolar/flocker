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
  fd1 <- simulate_flocker_data()
  mfd1 <- make_flocker_data(fd1$obs, fd1$unit_covs, fd1$event_covs)
  expect_silent(get_Z(example_flocker_model_single, new_data = mfd1))
#  expect_silent(get_Z(example_flocker_model_single, history_condition = FALSE, new_data = fd1$unit_covs))
  expect_error(get_Z(example_flocker_model_single, history_condition = TRUE, new_data = fd1$unit_covs))
})