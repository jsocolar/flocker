#' Compute pointwise log-likelihood matrix for a flocker_fit object
#' @param flocker_fit A flocker_fit object
#' @param new_data Optional new data at which to compute log likelihood. If
#'   supplied, it must be a flocker_data object produced by
#'   `make_flocker_data()`.
#' @param allow_new_levels allow new levels for random effect terms in 
#'    'new_data'? Will error if set to 'FALSE' and new levels are provided in 
#'    'new_data'.
#' @param sample_new_levels If new_data is provided and contains random effect 
#'    levels not present in the original data, how should predictions be 
#'    handled? See '?brms::prepare_predictions' for options.
#' @param draw_ids the draw ids to compute log-likelihoods for. Defaults to 
#'    using the full posterior. 
#' @return A posterior log-likelihood matrix, where iterations are rows and
#'   likelihood units are columns. For two-level models, columns are named with
#'   the level-two group names.
#' @details In one-level single-season models, columns are closure units (e.g.
#'   points or species-points; suitable for leave-one-unit-out CV). In
#'   multiseason models, columns are series (i.e. points or species-points;
#'   suitable for leave-one-series-out CV). In `twolevel_single` models,
#'   columns are level-two groups. In augmented models, columns are species
#'   (suitable for leave-one-species-out CV).
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
  
  if (lik_type %in% c("single")) {
    lps <- fitted_flocker(
      flocker_fit, draw_ids = draw_ids, new_data = new_data, 
      allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, 
      response = TRUE, unit_level = FALSE
    )
    
    psi_all <- first_column_draw_matrix(lps$linpred_occ)
    theta_all <- lps$linpred_det
    if (is.null(new_data)) {
      gp <- get_positions(flocker_fit)
      the_data <- flocker_fit$data
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
    components <- prepare_twolevel_postprocessing(
      flocker_fit, lik_type, draw_ids, new_data, allow_new_levels,
      sample_new_levels, include_detection = TRUE
    )
    ll <- log_lik_twolevel_single_from_components(
      components$occ_lp, components$det_lp, components$Omega_lp,
      components$group_id, components$group_known_present, components$obs
    ) |>
      t()
    colnames(ll) <- components$group_names
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
    init <- first_column_draw_matrix(lps2$linpred_occ)
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
    init_colo <- first_column_draw_matrix(colo)
    init_ex <- first_column_draw_matrix(ex)
    init <- init_colo / (init_colo + init_ex)
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
    init <- first_column_draw_matrix(lps2$linpred_occ)
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
    init_colo <- first_column_draw_matrix(colo)
    init_ex <- first_column_draw_matrix(ex)
    init <- init_colo / (init_colo + init_ex)
    det <- lps1$linpred_det
    ll <- log_lik_dynamic(init, colo, ex, obs, det) |>
      t()
  }
  ll
}

#' Compute grouped log-likelihoods for two-level single-season models
#' @noRd
log_lik_twolevel_single_from_components <- function(
    occ_lp, det_lp, Omega_lp, group_id, group_known_present, obs
    ) {
  assertthat::assert_that(inherits(obs, "matrix"))
  assertthat::assert_that(
    all(is.na(obs) | obs %in% 0:1)
  )
  assertthat::assert_that(is.matrix(occ_lp))
  assertthat::assert_that(length(dim(det_lp)) == 3)
  assertthat::assert_that(identical(dim(det_lp)[1:2], dim(obs)))
  assertthat::assert_that(is.matrix(Omega_lp))

  n_unit <- nrow(occ_lp)
  n_visit <- ncol(obs)
  n_draw <- ncol(occ_lp)
  n_group <- nrow(Omega_lp)
  assertthat::assert_that(dim(det_lp)[3] == n_draw)
  assertthat::assert_that(identical(dim(Omega_lp), c(n_group, n_draw)))

  log_el_0 <- matrixStats::rowSums2(log1p(-obs), na.rm = TRUE)
  log_unit_lik_available <- matrix(
    NA_real_, nrow = n_unit, ncol = n_draw
  )
  for(i in seq_len(n_draw)) {
    det_lp_draw <- matrix(
      det_lp[, , i, drop = FALSE], nrow = n_unit, ncol = n_visit
    )
    assertthat::assert_that(
      all(which(is.na(det_lp_draw)) %in% which(is.na(obs)))
    )
    log_event_lik_1 <- ifelse(
      obs == 1,
      log_inv_logit(det_lp_draw),
      log1m_inv_logit(det_lp_draw)
    )
    log_el_1 <- matrixStats::rowSums2(log_event_lik_1, na.rm = TRUE)
    log_unit_lik_available[, i] <- matrixStats::rowLogSumExps(
      cbind(
        log1m_inv_logit(occ_lp[, i]) + log_el_0,
        log_inv_logit(occ_lp[, i]) + log_el_1
      )
    )
  }

  out <- matrix(NA_real_, nrow = n_group, ncol = n_draw)
  for(g in seq_len(n_group)) {
    rows <- group_id == g
    log_p_y_available <- matrixStats::colSums2(
      log_unit_lik_available[rows, , drop = FALSE]
    )
    log_p_group_available <-
      log_inv_logit(Omega_lp[g, ]) + log_p_y_available
    out[g, ] <- if(group_known_present[g] == 1) {
      log_p_group_available
    } else {
      matrixStats::rowLogSumExps(
        cbind(log1m_inv_logit(Omega_lp[g, ]), log_p_group_available)
      )
    }
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
