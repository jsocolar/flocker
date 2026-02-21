## code to prepare `example_flocker_model_xx` datasets

setup_cache_path <- file.path(tempdir(), "flocker_test_cache.rds")

if (file.exists(setup_cache_path)) {
  cache <- readRDS(setup_cache_path)
  list2env(cache, envir = .GlobalEnv)
} else {
  set.seed(1)
  
  #### single-season rep-constant model ####
  suppressWarnings({
    fd <- simulate_flocker_data(n_pt = 20, n_sp = 5)
    mfd_single_C <- make_flocker_data(
      fd$obs, 
      fd$unit_covs
    )
    example_flocker_model_single_C <- flock(
      f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
      f_det = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
      flocker_data = mfd_single_C,
      cores = 1,
      iter = 8,
      save_warmup = FALSE
    )
    
    #### data augmented model ####
    fd <- simulate_flocker_data(augmented = TRUE,
                                n_sp = 10, n_pt = 20)
    mfd_aug <- make_flocker_data(
      fd$obs, 
      fd$unit_covs, 
      fd$event_covs,
      type = "augmented",
      n_aug = 10
    )
    example_flocker_model_aug <- flock(
      f_occ = ~ 0 + Intercept + (1 | ff_species),
      f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | ff_species),
      flocker_data = mfd_aug,
      augmented = TRUE,
      cores = 1,
      iter = 8,
      save_warmup = FALSE
    )
    
    #### multiseason colex with explicit inits ####
    fd <- simulate_flocker_data(
      n_pt = 20, n_sp = 1, n_season = 6, 
      multiseason = "colex", multi_init = "explicit",
      ragged_rep = TRUE, missing_seasons = TRUE
    )
    mfd_multi_colex_ex <- make_flocker_data(
      fd$obs, fd$unit_covs, fd$event_covs, type = "multi"
    )
    example_flocker_model_multi_colex_ex <- flock(
      f_occ = ~ 0 + Intercept + uc1,
      f_det = ~ 0 + Intercept + uc1 + ec1,
      f_col = ~ 0 + Intercept + uc1,
      f_ex = ~ 0 + Intercept + uc1,
      flocker_data = mfd_multi_colex_ex,
      multiseason = "colex",
      multi_init = "explicit",
      cores = 1,
      iter = 8,
      save_warmup = FALSE
    )
    
    #### multiseason colex with equilibrium inits ####
    fd <- simulate_flocker_data(
      n_pt = 20, n_sp = 1, n_season = 6, 
      multiseason = "colex", multi_init = "equilibrium",
      ragged_rep = TRUE, missing_seasons = TRUE
    )
    mfd_multi_colex_eq <- make_flocker_data(
      fd$obs, fd$unit_covs, fd$event_covs, type = "multi"
    )
    example_flocker_model_multi_colex_eq <- flock(
      f_det = ~ 0 + Intercept + uc1 + ec1,
      f_col = ~ 0 + Intercept + uc1,
      f_ex = ~ 0 + Intercept + uc1,
      flocker_data = mfd_multi_colex_eq,
      multiseason = "colex",
      multi_init = "equilibrium",
      cores = 1,
      iter = 8,
      save_warmup = FALSE
    )
    
    #### multiseason autologistic with explicit inits ####
    fd <- simulate_flocker_data(
      n_pt = 20, n_sp = 1, n_season = 6, 
      multiseason = "autologistic", multi_init = "explicit",
      ragged_rep = TRUE, missing_seasons = TRUE
    )
    mfd_multi_auto_ex <- make_flocker_data(
      fd$obs, fd$unit_covs, fd$event_covs, type = "multi"
    )
    example_flocker_model_multi_auto_ex <- flock(
      f_occ = ~ 0 + Intercept + uc1,
      f_det = ~ 0 + Intercept + uc1 + ec1,
      f_col = ~ 0 + Intercept + uc1,
      f_auto = ~ 0 + Intercept + uc1,
      flocker_data = mfd_multi_auto_ex,
      multiseason = "autologistic",
      multi_init = "explicit",
      cores = 1,
      iter = 8,
      save_warmup = FALSE
    )
    
    #### multiseason autologistic with equilibrium inits ####
    fd <- simulate_flocker_data(
      n_pt = 20, n_sp = 1, n_season = 6, 
      multiseason = "autologistic", multi_init = "equilibrium",
      ragged_rep = TRUE, missing_seasons = TRUE
    )
    mfd_multi_auto_eq <- make_flocker_data(
      fd$obs, fd$unit_covs, fd$event_covs, type = "multi"
    )
    example_flocker_model_multi_auto_eq <- flock(
      f_det = ~ 0 + Intercept + uc1 + ec1,
      f_col = ~ 0 + Intercept + uc1,
      f_auto = ~ 0 + Intercept + uc1,
      flocker_data = mfd_multi_auto_eq,
      multiseason = "autologistic",
      multi_init = "equilibrium",
      cores = 1,
      iter = 8,
      save_warmup = FALSE
    )
  })
  
  cache <- list(
    mfd_single_C = mfd_single_C,
    example_flocker_model_single_C = example_flocker_model_single_C,
    mfd_aug = mfd_aug,
    example_flocker_model_aug = example_flocker_model_aug,
    mfd_multi_colex_ex = mfd_multi_colex_ex,
    example_flocker_model_multi_colex_ex = example_flocker_model_multi_colex_ex,
    mfd_multi_colex_eq = mfd_multi_colex_eq,
    example_flocker_model_multi_colex_eq = example_flocker_model_multi_colex_eq,
    mfd_multi_auto_ex = mfd_multi_auto_ex,
    example_flocker_model_multi_auto_ex = example_flocker_model_multi_auto_ex,
    mfd_multi_auto_eq = mfd_multi_auto_eq,
    example_flocker_model_multi_auto_eq = example_flocker_model_multi_auto_eq
  )
  
  saveRDS(cache, setup_cache_path)
}

set.seed(1)