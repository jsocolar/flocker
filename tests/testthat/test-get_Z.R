test_that("get_Z gives valid returns", {
  Z_single <- get_Z(example_flocker_model_single2)
  expect_true(all(Z_single >= 0))
  expect_true(all(Z_single <= 1))
  expect_true(any(Z_single == 1))
  
  Z_single_C <- get_Z(example_flocker_model_single_C)
  expect_true(all(Z_single_C >= 0))
  expect_true(all(Z_single_C <= 1))
  expect_true(any(Z_single_C == 1))
  
  testthat::skip_on_cran()
  
  Z_augmented <- get_Z(example_flocker_model_aug)
  expect_named(Z_augmented, c("unit", "level2"))
  expect_true(all(Z_augmented$unit >= 0))
  expect_true(all(Z_augmented$unit <= 1))
  expect_true(any(Z_augmented$unit == 1))
  expect_true(all(Z_augmented$level2 >= 0))
  expect_true(all(Z_augmented$level2 <= 1))
  
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
  
  
  
  Z_single_nohist <- get_Z(example_flocker_model_single2, history_condition = FALSE)
  expect_true(all(Z_single_nohist >= 0))
  expect_true(all(Z_single_nohist <= 1))
  expect_true(sum(Z_single == 1) > sum(Z_single_nohist == 1))
  
  Z_single_C_nohist <- get_Z(example_flocker_model_single_C, history_condition = FALSE)
  expect_true(all(Z_single_C_nohist >= 0))
  expect_true(all(Z_single_C_nohist <= 1))
  expect_true(sum(Z_single_C == 1) > sum(Z_single_C_nohist == 1))
  
  Z_augmented_nohist <- get_Z(example_flocker_model_aug, history_condition = FALSE)
  expect_true(all(Z_augmented_nohist$unit >= 0))
  expect_true(all(Z_augmented_nohist$unit <= 1))
  expect_true(
    sum(Z_augmented$unit == 1) > sum(Z_augmented_nohist$unit == 1)
  )
  
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
  Z_single <- get_Z(example_flocker_model_single2, sample = TRUE)
  expect_true(all(Z_single %in% c(0, 1, NA)))
  
  Z_single_C <- get_Z(example_flocker_model_single_C, sample = TRUE)
  expect_true(all(Z_single_C %in% c(0, 1, NA)))
  
  testthat::skip_on_cran()
  
  Z_augmented <- get_Z(example_flocker_model_aug, sample = TRUE)
  expect_true(all(Z_augmented$unit %in% c(0, 1, NA)))
  expect_true(all(Z_augmented$level2 %in% c(0, 1, NA)))
  for(sp in seq_len(nrow(Z_augmented$level2))) {
    group_Z <- matrix(
      Z_augmented$level2[sp, ],
      nrow = dim(Z_augmented$unit)[1],
      ncol = dim(Z_augmented$unit)[3],
      byrow = TRUE
    )
    expect_true(all(Z_augmented$unit[, sp, ] <= group_Z))
  }
  
  Z_multi_colex_ex <- get_Z(example_flocker_model_multi_colex_ex, sample = TRUE)
  expect_true(all(Z_multi_colex_ex %in% c(0, 1, NA)))
  
  Z_multi_colex_eq <- get_Z(example_flocker_model_multi_colex_eq, sample = TRUE)
  expect_true(all(Z_multi_colex_eq %in% c(0, 1, NA)))
  
  Z_multi_auto_ex <- get_Z(example_flocker_model_multi_auto_ex, sample = TRUE)
  expect_true(all(Z_multi_auto_ex %in% c(0, 1, NA)))
  
  Z_multi_auto_eq <- get_Z(example_flocker_model_multi_auto_eq, sample = TRUE)
  expect_true(all(Z_multi_auto_eq %in% c(0, 1, NA)))

})

test_that("two-level state probabilities are calculated exactly", {
  psi <- matrix(
    c(0.2, 0.4,
      0.5, 0.7,
      0.3, 0.6,
      0.8, 0.9),
    nrow = 4, byrow = TRUE
  )
  theta <- array(0.5, dim = c(4, 2, 2))
  Omega <- matrix(c(0.35, 0.55, 0.65, 0.75), nrow = 2)
  group_id <- c(1L, 1L, 2L, 2L)
  group_known_present <- c(1L, 0L)
  obs <- matrix(
    c(1, 0,
      0, 0,
      0, 0,
      0, 0),
    nrow = 4, byrow = TRUE
  )

  unconditioned <- get_twolevel_states_from_components(
    psi, NULL, Omega, group_id, group_known_present,
    sample = FALSE, history_condition = FALSE
  )
  expect_equal(unconditioned$level2, Omega)
  expect_equal(unconditioned$unit, psi * Omega[group_id, , drop = FALSE])

  conditioned <- get_twolevel_states_from_components(
    psi, theta, Omega, group_id, group_known_present,
    sample = FALSE, history_condition = TRUE, obs = obs
  )
  unit_lik_available <- 1 - 0.75 * psi
  p_y_available <- apply(unit_lik_available[3:4, , drop = FALSE], 2, prod)
  expected_level2 <- rbind(
    c(1, 1),
    Omega[2, ] * p_y_available /
      ((1 - Omega[2, ]) + Omega[2, ] * p_y_available)
  )
  z_given_available <- psi * 0.25 / unit_lik_available
  z_given_available[1, ] <- 1
  expected_unit <- z_given_available *
    expected_level2[group_id, , drop = FALSE]
  expect_equal(conditioned$level2, expected_level2)
  expect_equal(conditioned$unit, expected_unit)

  set.seed(1)
  sampled <- get_twolevel_states_from_components(
    psi, theta, Omega, group_id, group_known_present,
    sample = TRUE, history_condition = TRUE, obs = obs
  )
  expect_true(all(sampled$level2 %in% 0:1))
  expect_true(all(sampled$unit %in% 0:1))
  expect_true(all(sampled$unit <= sampled$level2[group_id, , drop = FALSE]))
})

test_that("generic two-level states remain in original unit order", {
  unit_order <- c(3L, 1L, 4L, 2L)
  original_group <- c(2L, 2L, 1L, 1L)
  psi <- matrix(
    c(0.2, 0.3,
      0.4, 0.5,
      0.6, 0.7,
      0.8, 0.9),
    nrow = 4, byrow = TRUE
  )
  Omega <- matrix(c(0.25, 0.35, 0.75, 0.85), nrow = 2)
  Omega_by_unit <- Omega[original_group, , drop = FALSE]
  unit_lps <- structure(
    list(
      linpred_occ = psi,
      linpred_Omega = Omega_by_unit
    ),
    unit_level = TRUE
  )
  packed_data <- data.frame(
    ff_n_unit = c(4L, -99L, -99L, -99L),
    ff_n_group = c(2L, -99L, -99L, -99L),
    ff_n_unit_group = c(2L, 2L, -99L, -99L),
    ff_group_index = c(1L, 3L, 2L, 4L),
    ff_group_known_present = c(1L, 0L, -99L, -99L)
  )
  components <- format_twolevel_postprocessing(
    unit_lps, det_lps = NULL, obs = NULL, lik_type = "twolevel_single",
    flocker_data_data = packed_data,
    flocker_metadata = list(
      unit_order = unit_order,
      level2_group_names = c("group1", "group2")
    )
  )
  states <- get_twolevel_states_from_components(
    components$psi, components$theta, components$Omega,
    components$group_id, components$group_known_present,
    sample = FALSE, history_condition = FALSE
  )

  expect_equal(states$level2, Omega)
  expect_equal(states$unit, psi * Omega_by_unit)
})

test_that("get_Z returns both levels for augmented models", {
  testthat::skip_on_cran()

  n_group <- example_flocker_model_aug$data$ff_n_group[1]
  known_present <- example_flocker_model_aug$data$ff_group_known_present[
    seq_len(n_group)
  ]

  conditioned <- get_Z(
    example_flocker_model_aug, draw_ids = 1:2
  )
  expected_names <- c(
    paste0("observed_species", seq_len(10)),
    paste0("pseudospecies", seq_len(10))
  )
  expect_named(conditioned, c("unit", "level2"))
  expect_equal(dim(conditioned$level2), c(n_group, 2L))
  expect_equal(rownames(conditioned$level2), expected_names)
  expect_equal(dimnames(conditioned$unit)[[2]], expected_names)
  expect_true(all(conditioned$level2 >= 0 & conditioned$level2 <= 1))
  expect_true(all(conditioned$level2[known_present == 1, ] == 1))

  unconditioned <- get_Z(
    example_flocker_model_aug, draw_ids = 1:2,
    history_condition = FALSE
  )
  Omega <- fitted_flocker(
    example_flocker_model_aug, components = "Omega", draw_ids = 1:2,
    response = TRUE, unit_level = TRUE
  )
  expected <- group_level_Omega(
    Omega, "augmented", example_flocker_model_aug$data,
    get_flocker_metadata(example_flocker_model_aug)
  )
  expect_equal(as.vector(unconditioned$level2), as.vector(expected))

  one_draw <- get_Z(
    example_flocker_model_aug, draw_ids = 1
  )
  expect_equal(dim(one_draw$level2), c(n_group, 1L))

  sampled <- get_Z(
    example_flocker_model_aug, draw_ids = 1:2, sample = TRUE
  )
  expect_true(all(sampled$unit %in% 0:1))
  expect_true(all(sampled$level2 %in% 0:1))
})

test_that("new_data works as expected", {
  fd1 <- simulate_flocker_data(n_sp = 5)
  mfd1 <- make_flocker_data(fd1$obs, fd1$unit_covs, fd1$event_covs, quiet = TRUE)
  expect_silent(get_Z(example_flocker_model_single2, new_data = mfd1))
  expect_silent(get_Z(example_flocker_model_single2, history_condition = FALSE, new_data = fd1$unit_covs))
  expect_error(get_Z(example_flocker_model_single2, history_condition = TRUE, new_data = fd1$unit_covs))
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
