#' Compute unit-wise or series-wise log-likelihood matrix for a flocker_fit object
#' @param flocker_fit A flocker_fit object
#' @param new_data optional new data at which to compute log likelihood
#' @param allow_new_levels allow new levels for random effect terms in 
#'    'new_data'? Will error if set to 'FALSE' and new levels are provided in 
#'    'new_data'.
#' @param sample_new_levels If new_data is provided and contains random effect 
#'    levels not present in the original data, how should predictions be 
#'    handled? See '?brms::prepare_predictions' for options.
#' @param draw_ids the draw ids to compute log-likelihoods for. Defaults to 
#'    using the full posterior. 
#' @return A posterior log-likelihood matrix, where iterations are rows and 
#'    units, series, or species are columns.
#' @details In single-season models, rows are units (e.g. points or 
#'   species-points; suitable for leave-one-unit-out CV). In multiseason models, 
#'   rows are series (i.e. points or species-points, suitable for 
#'   leave-one-series-out CV). In augmented models, rows are species (suitable
#'   for leave-one-species-out CV).
#' @export
#' @examples 
#' \dontrun{
#' log_lik_flocker(example_flocker_model_single)
#' }

log_lik_flocker <- function(
    flocker_fit, new_data = NULL, allow_new_levels = FALSE,
    sample_new_levels = "uncertainty", draw_ids = NULL
    ) {
  assertthat::assert_that(is_flocker_fit(flocker_fit))
  ndraws <- brms::ndraws(flocker_fit)
  if(!is.null(draw_ids)){
    assertthat::assert_that(
      max(draw_ids) <= ndraws,
      msg = "some requested draw ids greater than total number of available draws"
    )
    ndraws <- length(draw_ids)
  }
  
  if(!is.null(new_data)){
    assertthat::assert_that(
      is_flocker_data(new_data),
      msg = paste0(
        "new_data, if provided, must be a flocker_data object formatted ",
        "by `make_flocker_data`"
      )
    )
    a <- attributes(flocker_fit)
    fit_datatype <- a$data_type
    assertthat::assert_that(
      new_data$type == fit_datatype,
      msg = "the data type in new_data does not match the data type in flocker_fit"
    )
  }
  
  lik_type <- type_flocker_fit(flocker_fit)
  if (is.null(new_data) && "flocker_data" %in% names(attributes(flocker_fit))) {
    original_data <- attributes(flocker_fit)$flocker_data$data
  } else {
    original_data <- flocker_fit$data
  }
  
  if (lik_type %in% c("single")) {
    lps <- fitted_flocker(
      flocker_fit, draw_ids = draw_ids, new_data = new_data, 
      allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, 
      response = TRUE, unit_level = FALSE
    )
    
    psi_all <- lps$linpred_occ[ , 1, ] # first index is unit, second is visit, third is draw
    theta_all <- lps$linpred_det
    if (is.null(new_data)) {
      gp <- get_positions(flocker_fit)
      the_data <- original_data
    } else {
      gp <- get_positions(new_data)
      the_data <- new_data$data
    }
    
    obs <- new_matrix(gp, the_data$ff_y[gp])

    # get emission likelihoods
    el_0 <- el_1 <- matrix(nrow = nrow(psi_all), ncol = ncol(psi_all))
    
    for(i in seq_len(ncol(psi_all))){
      el_0[ , i] <- emission_likelihood(0, obs, theta_all[,,i])
      el_1[ , i] <- emission_likelihood(1, obs, theta_all[,,i])
    }
      
    assertthat::assert_that(identical(dim(psi_all), dim(el_0)))
      
    ll <- log((psi_all * el_1) + (1 - psi_all) * el_0) |>
      t()
  } else if (lik_type == "single_C") {
    ll <- brms::log_lik(
      flocker_fit, newdata = new_data$data, draw_ids = draw_ids, #note that if new_data is NULL, new_data$data is also NULL
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels)
  } else if (lik_type %in% c("twolevel_single", "augmented")) {
    if(is.null(new_data)){
      gp <- get_positions(flocker_fit)
      obs <- new_array(gp, original_data$ff_y[gp])
      the_data <- original_data
    } else {
      gp <- get_positions(new_data)
      obs <- new_array(gp, new_data$data$ff_y[gp])
      the_data <- new_data$data
    }

    lps <- fitted_flocker(
      flocker_fit, draw_ids = draw_ids, new_data = new_data, 
      allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, 
      response = TRUE, unit_level = FALSE
    )
    if(lik_type == "augmented") {
      n_point <- dim(obs)[1]
      n_visit <- dim(obs)[2]
      n_species <- dim(obs)[3]
      n_unit <- the_data$ff_n_unit[1]
      unit_rows <- seq_len(n_unit)
      site_id <- the_data$ff_site[unit_rows]
      species_id <- the_data$ff_group[unit_rows]
      psi_all_array <- lps$linpred_occ[ , 1, , ]
      psi_all <- matrix(NA_real_, nrow = n_unit, ncol = ndraws)
      theta_all <- array(NA_real_, dim = c(n_unit, n_visit, ndraws))
      obs_use <- matrix(NA_real_, nrow = n_unit, ncol = n_visit)
      for(i in unit_rows) {
        psi_all[i, ] <- psi_all_array[site_id[i], species_id[i], ]
        obs_use[i, ] <- obs[site_id[i], , species_id[i]]
        theta_all[i, , ] <- lps$linpred_det[site_id[i], , species_id[i], ]
      }
    } else {
      n_unit <- the_data$ff_n_unit[1]
      unit_rows <- seq_len(n_unit)
      orig_unit <- the_data$ff_orig_unit[unit_rows]
      psi_all <- lps$linpred_occ[ , 1, ][orig_unit, , drop = FALSE]
      theta_all <- lps$linpred_det[orig_unit, , , drop = FALSE]
      obs_use <- obs[orig_unit, , drop = FALSE]
    }
    Omega <- lps$linpred_Omega
    group_id <- the_data$ff_group[seq_len(the_data$ff_n_unit[1])]
    group_known_present <- the_data$ff_group_known_present[
      seq_len(the_data$ff_n_group[1])
    ]
    ll <- log_lik_twolevel_single_from_components(
      psi_all, theta_all, Omega, group_id, group_known_present, obs_use
    ) |>
      t()
  } else if (lik_type %in% c("multi_colex")) {
    if(is.null(new_data)){
      gp <- get_positions(flocker_fit)
      obs <- new_array(gp, flocker_fit$data$ff_y[gp])
    } else {
      gp <- get_positions(new_data)
      obs <- new_array(gp, new_data$data$ff_y[gp])
    }
    
    lps1 <- fitted_flocker(
      flocker_fit, components = c("det"),
      new_data = new_data,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels,
      draw_ids = draw_ids
    )
    lps2 <- fitted_flocker(
      flocker_fit, components = c("occ", "colo", "ex"),
      new_data = new_data,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels,
      draw_ids = draw_ids, unit_level = TRUE
    )
    init <- lps2$linpred_occ[,1,]
    colo <- lps2$linpred_col
    ex <- lps2$linpred_ex
    det <- lps1$linpred_det
    ll <- log_lik_dynamic(init, colo, ex, obs, det) |>
      t()
  } else if (lik_type %in% c("multi_colex_eq")) {
    if(is.null(new_data)){
      gp <- get_positions(flocker_fit)
      obs <- new_array(gp, flocker_fit$data$ff_y[gp])
    } else {
      gp <- get_positions(new_data)
      obs <- new_array(gp, new_data$data$ff_y[gp])
    }
    
    lps1 <- fitted_flocker(
      flocker_fit, components = c("det"),
      new_data = new_data,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels,
      draw_ids = draw_ids
    )
    lps2 <- fitted_flocker(
      flocker_fit, components = c("col", "ex"),
      new_data = new_data,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels,
      draw_ids = draw_ids, unit_level = TRUE
    )
    colo <- lps2$linpred_col
    ex <- lps2$linpred_ex
    init <- colo[,1,] / (colo[,1,] + ex[,1,])
    det <- lps1$linpred_det
    ll <- log_lik_dynamic(init, colo, ex, obs, det) |>
      t()
  } else if (lik_type == "multi_autologistic") {
    if(is.null(new_data)){
      gp <- get_positions(flocker_fit)
      obs <- new_array(gp, flocker_fit$data$ff_y[gp])
    } else {
      gp <- get_positions(new_data)
      obs <- new_array(gp, new_data$data$ff_y[gp])
    }
    
    lps1 <- fitted_flocker(
      flocker_fit, components = c("det"), 
      new_data = new_data,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels,
      draw_ids = draw_ids
    )
    lps2 <- fitted_flocker(
      flocker_fit, components = c("occ", "colo", "auto"),
      new_data = new_data,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels,
      draw_ids = draw_ids, unit_level = TRUE
    )
    init <- lps2$linpred_occ[,1,]
    colo <- lps2$linpred_col
    ex <- 1 - boot::inv.logit(boot::logit(colo) + lps2$linpred_auto)
    det <- lps1$linpred_det
    ll <- log_lik_dynamic(init, colo, ex, obs, det) |>
      t()
  } else if (lik_type == "multi_autologistic_eq") {
    if(is.null(new_data)){
      gp <- get_positions(flocker_fit)
      obs <- new_array(gp, flocker_fit$data$ff_y[gp])
    } else {
      gp <- get_positions(new_data)
      obs <- new_array(gp, new_data$data$ff_y[gp])
    }
    
    lps1 <- fitted_flocker(
      flocker_fit, components = c("det"), 
      new_data = new_data,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels,
      draw_ids = draw_ids
    )
    lps2 <- fitted_flocker(
      flocker_fit, components = c("col", "auto"),
      new_data = new_data,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels,
      draw_ids = draw_ids, unit_level = TRUE
    )
    colo <- lps2$linpred_col
    ex <- 1 - boot::inv.logit(boot::logit(colo) + lps2$linpred_auto)
    init <- colo[,1,] / (colo[,1,] + ex[,1,])
    det <- lps1$linpred_det
    ll <- log_lik_dynamic(init, colo, ex, obs, det) |>
      t()
  }
  ll
}

#' Compute grouped log-likelihoods for two-level single-season models
#' @noRd
log_lik_twolevel_single_from_components <- function(
    psi_all, theta_all, Omega, group_id, group_known_present, obs
    ) {
  n_unit <- nrow(psi_all)
  n_draw <- ncol(psi_all)
  n_group <- nrow(Omega)
  el_0 <- el_1 <- matrix(NA, nrow = n_unit, ncol = n_draw)
  
  for(i in seq_len(n_draw)){
    el_0[ , i] <- emission_likelihood(0, obs, theta_all[,,i])
    el_1[ , i] <- emission_likelihood(1, obs, theta_all[,,i])
  }
  
  unit_lik_available <- (1 - psi_all) * el_0 + psi_all * el_1
  out <- matrix(NA, nrow = n_group, ncol = n_draw)
  
  for(g in seq_len(n_group)) {
    rows <- group_id == g
    p_y_available <- apply(unit_lik_available[rows, , drop = FALSE], 2, prod)
    p_y_unavailable <- as.numeric(group_known_present[g] == 0)
    out[g, ] <- log(
      Omega[g, ] * p_y_available +
        (1 - Omega[g, ]) * p_y_unavailable
    )
  }
  
  out
}

#' A log-likelihood function for the rep-constant occupancy model, sufficient for
#' \code{brms::loo()} to work. 
#' @param i Observation id
#' @param prep Output of \code{brms::prepare_predictions}. See brms custom families
#' vignette at 
#' https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html
#' @return The log-likelihood for observation i
#' @noRd
log_lik_occupancy_single_C <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  occ <- brms::get_dpar(prep, "occ", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  return(occupancy_single_C_lpmf(y, mu, occ, trials))
}

#' An R implementation of the rep constant lpmf without the binomial coefficient.
#' @param y number of detections
#' @param mu logit-scale detection probability
#' @param occ logit-scale occupancy probability
#' @param trials number of reps
#' @return The log-likelihood
#' @details By omitting the binomial coefficient, the log-likelihood is directly
#'   comparable to log-likelihoods from rep-varying models.
#' @noRd
occupancy_single_C_lpmf <- Vectorize(
  function (y, mu, occ, trials) {
    if (y == 0) {
      out <- 
        matrixStats::logSumExp(
          c(log1m_inv_logit(occ), log_inv_logit(occ) + trials * log1m_inv_logit(mu))
        )
    } else {
      out <- log_inv_logit(occ) + 
        y * log_inv_logit(mu) + 
        (trials - y) * log1m_inv_logit(mu)
      
    }
    out
  }
)

#' An R implementation of occupancy_C_lpmf including the binomial coefficient.
#' Not currently in use.
#' @param y number of detections
#' @param mu logit-scale detection probability
#' @param occ logit-scale occupancy probability
#' @param trials number of reps
#' @return The log-likelihood
#' @noRd
occupancy_single_C_lpmf_with_coef <- Vectorize(
  function (y, mu, occ, trials) {
    if (y == 0) {
      out <- 
        matrixStats::logSumExp(
          c(log1m_inv_logit(occ), log_inv_logit(occ) + trials * log1m_inv_logit(mu))
        )
    } else {
      out <- log_inv_logit(occ) + 
        log(choose(trials, y)) +
        y * log_inv_logit(mu) + 
        (trials - y) * log1m_inv_logit(mu)
      
    }
    out
  }
)

#' log_lik_dynamic
#' @param init matrix of initial occupancy probabilities (rows are series and 
#'   columns are draws)
#' @param colo array of colonization probabilities (rows are series, columns are
#'  timesteps, slices are draws)
#' @param ex array of extinction probabilities (rows are series, columns are
#'  timesteps, slices are draws)
#' @param obs Array of observations. Rows are series, columns are visits, and 
#'  slices are timesteps.
#' @param det Array of detection probabilities. Rows are series, columns are 
#'  visits, slices are timesteps, and slice_2s are draws.
#' @return A log-likelihood matrix. Rows are sites, 
#'   columns are draws.
#' @noRd
log_lik_dynamic <- function(init, colo, ex, obs, det){
  assertthat::assert_that(is.matrix(init))
  assertthat::assert_that(is.array(colo))
  assertthat::assert_that(is.array(ex))
  assertthat::assert_that(is.array(obs))
  assertthat::assert_that(is.array(det))
  
  nsite <- nrow(init)
  ndraw <- ncol(init)
  
  assertthat::assert_that(identical(dim(obs), dim(det)[1:3]))
  assertthat::assert_that(isTRUE(dim(det)[4] == ndraw))
  assertthat::assert_that(identical(dim(colo), dim(ex)))
  assertthat::assert_that(isTRUE(nrow(colo) == nsite))
  assertthat::assert_that(isTRUE(nslice(colo) == ndraw))
  assertthat::assert_that(isTRUE(nrow(obs) == nrow(colo)))
  assertthat::assert_that(isTRUE(nslice(obs) == ncol(colo)))
  
  out <- matrix(nrow = nsite, ncol = ndraw)
  
  for(i in 1:nsite){
    for(j in 1:ndraw){
      el0 <- emission_likelihood(0, t(obs[i,,]), t(det[i,,,j]))
      el1 <- emission_likelihood(1, t(obs[i,,]), t(det[i,,,j]))
      forward_output <- forward_algorithm(el0, el1, init[i,j], colo[i,,j], ex[i,,j]) |>
        rowSums()
      # forward_output should be strictly decreasing except when there is a season
      # with no visits. In testing, this occasionally caused numerical issues where
      # forward_output seemed to increase minimally, thus the tolerance of 1e-15
      # below.
      assertthat::assert_that(
        all(diff(forward_output) <= 1e-15, na.rm = TRUE),
        msg = "forward_output tolerance error; please report a bug at www.github.com/jsocolar/flocker"
        )
      out[i,j] <- log(min(forward_output, na.rm = TRUE))
    }
  }
  assertthat::assert_that(!(NA %in% out))
  out
}
