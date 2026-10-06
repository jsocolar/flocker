## code to prepare `example_flocker_model_xx` datasets

# Set FLOCKER_TEST_CACHE_DIR to reuse fitted fixtures across R sessions. Bump
# this version whenever fixture definitions or fit-producing code changes in a
# way not captured by the package-version checks below.
fixture_cache_version <- 1L
fixture_cache_dir <- Sys.getenv("FLOCKER_TEST_CACHE_DIR", unset = "")
if (!nzchar(fixture_cache_dir)) {
  fixture_cache_dir <- tempdir()
}
if (
  !dir.exists(fixture_cache_dir) &&
  !dir.create(fixture_cache_dir, recursive = TRUE, showWarnings = FALSE)
) {
  stop("Unable to create the flocker test-fixture cache directory")
}

full_fixture_cache <- identical(Sys.getenv("NOT_CRAN"), "true")
fixture_cache_variant <- if(full_fixture_cache) "full" else "reduced"
setup_cache_path <- file.path(
  fixture_cache_dir,
  paste0(
    "flocker_test_fixtures_v", fixture_cache_version, "_",
    fixture_cache_variant, ".rds"
  )
)
fixture_cache_signature <- list(
  fixture_cache_version = fixture_cache_version,
  full_fixture_cache = full_fixture_cache,
  R = paste(R.version$major, R.version$minor, sep = "."),
  platform = R.version$platform,
  brms = as.character(utils::packageVersion("brms")),
  rstan = as.character(utils::packageVersion("rstan")),
  StanHeaders = as.character(utils::packageVersion("StanHeaders"))
)

fixture_cache_record <- if(file.exists(setup_cache_path)) {
  tryCatch(readRDS(setup_cache_path), error = function(e) NULL)
} else {
  NULL
}
cache <- if(
  is.list(fixture_cache_record) &&
  identical(fixture_cache_record$signature, fixture_cache_signature) &&
  is.list(fixture_cache_record$fixtures)
) {
  fixture_cache_record$fixtures
} else {
  NULL
}

if (!is.null(cache)) {
  list2env(cache, envir = .GlobalEnv)
} else {
  set.seed(1)
  
  suppressWarnings({
    example_flocker_model_single2 <- flock(
      f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
      f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | species),
      flocker_data = mfd_single,
      cores = 2,
      chains = 2,
      iter = 8,
      warmup = 6,
      save_warmup = FALSE,
      refresh = 0,
      silent = 2
    )
  })

  #### generic two-level single-season model ####
  obs_twolevel <- matrix(
    c(1, 0,
      0, 0,
      0, 0,
      0, 0,
      0, 1,
      0, 0),
    nrow = 6, byrow = TRUE
  )
  unit_covs_twolevel <- data.frame(
    group = factor(rep(letters[1:3], each = 2))
  )
  mfd_twolevel <- make_flocker_data(
    obs_twolevel, unit_covs_twolevel, type = "twolevel_single",
    level2_group = "group", quiet = TRUE
  )
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    example_flocker_model_twolevel <- NULL
  } else {
    suppressWarnings({
      example_flocker_model_twolevel <- flock(
        f_occ = ~ 1,
        f_det = ~ 1,
        f_meta = ~ 1,
        flocker_data = mfd_twolevel,
        cores = 2,
        chains = 2,
        iter = 8,
        warmup = 6,
        save_warmup = FALSE,
        refresh = 0,
        silent = 2
      )
    })
  }
  
  #### single-season rep-constant model ####
  suppressWarnings({
    fd <- simulate_flocker_data(n_pt = 20, n_sp = 5)
    mfd_single_C <- make_flocker_data(
      fd$obs, 
      fd$unit_covs,
      quiet = TRUE
    )
    example_flocker_model_single_C <- flock(
      f_occ = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
      f_det = ~ 0 + Intercept + uc1 + (1 + uc1 | species),
      flocker_data = mfd_single_C,
      cores = 2,
      iter = 8,
      save_warmup = FALSE,
      refresh = 0,
      silent = 2
    )
    
    #### data augmented model ####
    fd <- simulate_flocker_data(augmented = TRUE,
                                n_sp = 10, n_pt = 20)
    dimnames(fd$obs)[[3]] <- paste0("observed_species", seq_len(10))
    mfd_aug <- make_flocker_data(
      fd$obs, 
      fd$unit_covs, 
      fd$event_covs,
      type = "augmented",
      n_aug = 10,
      quiet = TRUE
    )
    
    if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
      example_flocker_model_aug <- NULL
    } else {
      example_flocker_model_aug <- flock(
        f_occ = ~ 0 + Intercept + (1 | ff_species),
        f_det = ~ 0 + Intercept + uc1 + ec1 + (1 + uc1 + ec1 | ff_species),
        flocker_data = mfd_aug,
        augmented = TRUE,
        cores = 2,
        iter = 8,
        save_warmup = FALSE,
        refresh = 0,
        silent = 2
      )
    }
    
    #### multiseason colex with explicit inits ####
    fd <- simulate_flocker_data(
      n_pt = 20, n_sp = 1, n_season = 6, 
      multiseason = "colex", multi_init = "explicit",
      ragged_rep = TRUE, missing_seasons = TRUE
    )
    mfd_multi_colex_ex <- make_flocker_data(
      fd$obs, fd$unit_covs, fd$event_covs, type = "multi",
      quiet = TRUE
    )
    
    if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
      example_flocker_model_multi_colex_ex <- NULL
    } else {
      example_flocker_model_multi_colex_ex <- flock(
        f_occ = ~ 0 + Intercept + uc1,
        f_det = ~ 0 + Intercept + uc1 + ec1,
        f_col = ~ 0 + Intercept + uc1,
        f_ex = ~ 0 + Intercept + uc1,
        flocker_data = mfd_multi_colex_ex,
        multiseason = "colex",
        multi_init = "explicit",
        cores = 2,
        iter = 8,
        save_warmup = FALSE,
        refresh = 0,
        silent = 2
      )
    }
    
    
    #### multiseason colex with equilibrium inits ####
    fd <- simulate_flocker_data(
      n_pt = 20, n_sp = 1, n_season = 6, 
      multiseason = "colex", multi_init = "equilibrium",
      ragged_rep = TRUE, missing_seasons = TRUE
    )
    mfd_multi_colex_eq <- make_flocker_data(
      fd$obs, fd$unit_covs, fd$event_covs, type = "multi",
      quiet = TRUE
    )
    if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
      example_flocker_model_multi_colex_eq <- NULL
    } else {
      example_flocker_model_multi_colex_eq <- flock(
        f_det = ~ 0 + Intercept + uc1 + ec1,
        f_col = ~ 0 + Intercept + uc1,
        f_ex = ~ 0 + Intercept + uc1,
        flocker_data = mfd_multi_colex_eq,
        multiseason = "colex",
        multi_init = "equilibrium",
        cores = 2,
        iter = 8,
        save_warmup = FALSE,
        refresh = 0,
        silent = 2
      )
    }
    
    
    #### multiseason autologistic with explicit inits ####
    fd <- simulate_flocker_data(
      n_pt = 20, n_sp = 1, n_season = 6, 
      multiseason = "autologistic", multi_init = "explicit",
      ragged_rep = TRUE, missing_seasons = TRUE
    )
    mfd_multi_auto_ex <- make_flocker_data(
      fd$obs, fd$unit_covs, fd$event_covs, type = "multi",
      quiet = TRUE
    )
    if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
      example_flocker_model_multi_auto_ex <- NULL
    } else {
      example_flocker_model_multi_auto_ex <- flock(
        f_occ = ~ 0 + Intercept + uc1,
        f_det = ~ 0 + Intercept + uc1 + ec1,
        f_col = ~ 0 + Intercept + uc1,
        f_auto = ~ 0 + Intercept + uc1,
        flocker_data = mfd_multi_auto_ex,
        multiseason = "autologistic",
        multi_init = "explicit",
        cores = 2,
        iter = 8,
        save_warmup = FALSE,
        refresh = 0,
        silent = 2
      )
    }
    
    
    #### multiseason autologistic with equilibrium inits ####
    fd <- simulate_flocker_data(
      n_pt = 20, n_sp = 1, n_season = 6, 
      multiseason = "autologistic", multi_init = "equilibrium",
      ragged_rep = TRUE, missing_seasons = TRUE
    )
    mfd_multi_auto_eq <- make_flocker_data(
      fd$obs, fd$unit_covs, fd$event_covs, type = "multi",
      quiet = TRUE
    )
    if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
      example_flocker_model_multi_auto_eq <- NULL
    } else {
      example_flocker_model_multi_auto_eq <- flock(
        f_det = ~ 0 + Intercept + uc1 + ec1,
        f_col = ~ 0 + Intercept + uc1,
        f_auto = ~ 0 + Intercept + uc1,
        flocker_data = mfd_multi_auto_eq,
        multiseason = "autologistic",
        multi_init = "equilibrium",
        cores = 2,
        iter = 8,
        save_warmup = FALSE,
        refresh = 0,
        silent = 2
      )
    }
  })
  
  cache <- list(
    example_flocker_model_single2 = example_flocker_model_single2,
    mfd_twolevel = mfd_twolevel,
    example_flocker_model_twolevel = example_flocker_model_twolevel,
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
  
  if (
    !dir.exists(fixture_cache_dir) &&
    !dir.create(fixture_cache_dir, recursive = TRUE, showWarnings = FALSE)
  ) {
    stop("Unable to create the flocker test-fixture cache directory")
  }
  fixture_cache_tmp <- tempfile(
    pattern = "flocker_test_fixtures_",
    tmpdir = fixture_cache_dir,
    fileext = ".rds"
  )
  saveRDS(
    list(signature = fixture_cache_signature, fixtures = cache),
    fixture_cache_tmp
  )
  if (file.exists(setup_cache_path)) {
    unlink(setup_cache_path)
  }
  if (!file.rename(fixture_cache_tmp, setup_cache_path)) {
    unlink(fixture_cache_tmp)
    stop("Unable to write the flocker test-fixture cache")
  }
}

set.seed(1)
