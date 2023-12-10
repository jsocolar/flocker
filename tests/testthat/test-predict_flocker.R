test_that("predict flocker works as expected", {
  p_single <- predict_flocker(example_flocker_model_single)
  expect_true(all(p_single %in% c(0,1)))
  
  p_single_C <- predict_flocker(example_flocker_model_single_C)
  expect_true(all(p_single_C %in% c(0,1)))
  
  p_augmented <- predict_flocker(example_flocker_model_aug)
  expect_true(all(p_augmented %in% c(0,1)))
  
  p_multi_colex_ex <- predict_flocker(example_flocker_model_multi_colex_ex)
  expect_true(all(p_multi_colex_ex %in% c(0,1,NA)))
  
  p_multi_colex_eq <- predict_flocker(example_flocker_model_multi_colex_eq)
  expect_true(all(p_multi_colex_eq %in% c(0,1,NA)))
  
  p_multi_auto_ex <- predict_flocker(example_flocker_model_multi_auto_ex)
  expect_true(all(p_multi_auto_ex %in% c(0,1,NA)))

  p_multi_auto_eq <- predict_flocker(example_flocker_model_multi_auto_eq)
  expect_true(all(p_multi_auto_eq %in% c(0,1,NA)))
  
  p_single_hist <- predict_flocker(example_flocker_model_single, history_condition = TRUE)
  expect_true(all(p_single_hist %in% c(0,1)))
  
  p_single_C_hist <- predict_flocker(example_flocker_model_single_C, history_condition = TRUE)
  expect_true(all(p_single_C_hist %in% c(0,1)))
  
  p_augmented_hist <- predict_flocker(example_flocker_model_aug, history_condition = TRUE)
  expect_true(all(p_augmented_hist %in% c(0,1)))
  
  p_multi_colex_ex_hist <- predict_flocker(example_flocker_model_multi_colex_ex, history_condition = TRUE)
  expect_true(all(p_multi_colex_ex_hist %in% c(0,1,NA)))
  
  p_multi_colex_eq_hist <- predict_flocker(example_flocker_model_multi_colex_eq, history_condition = TRUE)
  expect_true(all(p_multi_colex_eq_hist %in% c(0,1,NA)))
  
  p_multi_auto_ex_hist <- predict_flocker(example_flocker_model_multi_auto_ex, history_condition = TRUE)
  expect_true(all(p_multi_auto_ex_hist %in% c(0,1,NA)))
  
  p_multi_auto_eq_hist <- predict_flocker(example_flocker_model_multi_auto_eq, history_condition = TRUE)
  expect_true(all(p_multi_auto_eq_hist %in% c(0,1,NA)))
  
})

test_that("new_data works as expected", {
  fd1 <- simulate_flocker_data(n_sp = 5, n_pt = 5)
  mfd1 <- make_flocker_data(fd1$obs, fd1$unit_covs, fd1$event_covs)
  expect_silent(predict_flocker(example_flocker_model_single, new_data = mfd1))
  expect_error(predict_flocker(example_flocker_model_single, history_condition = FALSE, new_data = fd1$unit_covs))
  
  fd2 <- simulate_flocker_data(n_sp = 5, n_pt = 5, n_season = 3, multiseason = "colex", multi_init = "explicit")
  mfd2 <- make_flocker_data(fd2$obs, fd2$unit_covs, fd2$event_covs, type = "multi")
  expect_silent(predict_flocker(example_flocker_model_multi_colex_ex, new_data = mfd2))
  expect_silent(predict_flocker(example_flocker_model_multi_auto_eq, new_data = mfd2))
})

test_that("mixed checks work", {
  p_single_mixed <- predict_flocker(
    example_flocker_model_single, mixed = TRUE, 
    allow_new_levels = TRUE, sample_new_levels = "gaussian"
    )
  expect_true(all(p_single_mixed %in% c(0,1)))
  
  sfd <- simulate_flocker_data(n_pt = 5, n_sp = 5, n_rep = 2)
  fd <- make_flocker_data(sfd$obs, sfd$unit_covs, sfd$event_covs)
  p_single_mixed <- predict_flocker(
    example_flocker_model_single, mixed = TRUE, 
    allow_new_levels = TRUE, sample_new_levels = "gaussian",
    new_data = fd
  )
  expect_true(all(p_single_mixed %in% c(0,1)))
})
