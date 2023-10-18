## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----simulate, eval = FALSE---------------------------------------------------
#  library(flocker)
#  multi_data <- simulate_flocker_data(
#    n_season = 3,
#    n_pt = 300,
#    n_sp = 1,
#    multiseason = "colex",
#    multi_init = "explicit",
#    seed = 1
#    )
#  
#  fd <- make_flocker_data(
#    multi_data$obs,
#    multi_data$unit_covs,
#    multi_data$event_covs,
#    type = "multi",
#    quiet = TRUE
#    )
#  

## ----colex-explicit, eval = FALSE---------------------------------------------
#  multi_colex <- flock(
#    f_occ = ~ uc1,
#    f_det = ~ uc1 + ec1,
#    f_col = ~ uc1,
#    f_ex = ~ uc1,
#    flocker_data = fd,
#    multiseason = "colex",
#    multi_init = "explicit",
#    backend = "cmdstanr",
#    cores = 4
#  )
#  

## ----colex-equilibrium, eval = FALSE------------------------------------------
#  multi_colex_eq <- flock(
#    f_det = ~ uc1 + ec1,
#    f_col = ~ uc1,
#    f_ex = ~ uc1,
#    flocker_data = fd,
#    multiseason = "colex",
#    multi_init = "equilibrium",
#    backend = "cmdstanr",
#    cores = 4
#  )
#  

## ----auto-explicit, eval = FALSE----------------------------------------------
#  multi_auto <- flock(
#    f_occ = ~ uc1,
#    f_det = ~ uc1 + ec1,
#    f_col = ~ uc1,
#    f_auto = ~ 1,
#    flocker_data = fd,
#    multiseason = "autologistic",
#    multi_init = "explicit",
#    backend = "cmdstanr",
#    cores = 4
#  )
#  

## ----auto-equilirium, eval = FALSE--------------------------------------------
#  multi_auto_eq <- flock(
#    f_det = ~ uc1 + ec1,
#    f_col = ~ uc1,
#    f_auto = ~ 1,
#    flocker_data = fd,
#    multiseason = "autologistic",
#    multi_init = "equilibrium",
#    backend = "cmdstanr",
#    cores = 4
#  )
#  

## ----Z, eval = FALSE----------------------------------------------------------
#  Z <- get_Z(multi_colex_eq)
#  

## ----predict, eval = FALSE----------------------------------------------------
#  pp <- predict_flocker(multi_colex_eq)
#  

## ----loo, eval = FALSE--------------------------------------------------------
#  loo_out <- loo_compare_flocker(
#    list(multi_colex, multi_colex_eq, multi_auto, multi_auto_eq)
#  )
#  loo_out
#  

