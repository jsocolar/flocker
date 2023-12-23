#' Get posterior distribution of the Z matrix
#' @param flocker_fit A flocker_fit object
#' @param draw_ids Vector of indices of the posterior draws to be 
#'     used. If `NULL` (the default) all draws are used in their native order. 
#' @param history_condition Should the posterior distribution for Z directly 
#'     condition on the observed detection history (`TRUE`) or not (`FALSE`)?
#'     For example, at sites with at least one detection, the true occupancy 
#'     state conditioned on the history is one with absolute certainty. Without 
#'     directly conditioning on the history, the occupancy state is controlled 
#'     by the posterior distribution for the occupancy probability psi.
#' @param sample Should the return be posterior probabilities of occupancy (FALSE),
#'     or bernoulli samples from those probabilities (TRUE)
#' @param new_data Optional new data at which to predict the Z matrix. Can be 
#'     the output of `make_flocker_data` or the `unit_covs` input to 
#'     `make_flocker_data` provided that `history_condition` is `FALSE` and the 
#'     occupancy model is a single-season, non-augmented model.
#' @param allow_new_levels allow new levels for random effect terms in `new_data`?
#'     Will error if set to `FALSE` and new levels are provided in `new_data`.
#' @param sample_new_levels If `new_data` is provided and contains random effect
#'     levels not present in the original data, how should predictions be
#'     handled? Passed directly to `brms::prepare_predictions`, which see. 
#' @return The posterior Z matrix in the shape of the first visit in `obs` as
#'     passed to make_flocker_data, with posterior iterations stacked along the
#'     final dimension
#' @export
#' @examples
#' \donttest{
#' get_Z(example_flocker_model_single)
#' }
get_Z <- function (flocker_fit, draw_ids = NULL, history_condition = TRUE, 
                   sample = FALSE,
                   new_data = NULL, allow_new_levels = FALSE, 
                   sample_new_levels = "uncertainty") {
  # Input checking and processing
  assertthat::assert_that(
    is_one_logical(history_condition),
    msg = "history_condition must be a single logical value"
  )
  assertthat::assert_that(
    !history_condition | is_flocker_data(new_data) | is.null(new_data),
    msg = "conditioning on observed detection histories while using the `new_data` argument requires passing a `flocker_data` object, not a dataframe"
  )
  lik_type <- type_flocker_fit(flocker_fit)
  is_multi <- fdtl()$data_output_type[fdtl()$model_type == lik_type] == "multi"
  is_aug <- fdtl()$data_output_type[fdtl()$model_type == lik_type] == "augmented"
  assertthat::assert_that(
    !is_multi | is_flocker_data(new_data) | is.null(new_data),
    msg = "using the `new_data` argument for a multiseason model requires passing a `flocker_data` object, not a dataframe"
  )
  assertthat::assert_that(
    !is_aug | is_flocker_data(new_data) | is.null(new_data),
    msg = "using the `new_data` argument for a data-augmented model requires passing a `flocker_data` object, not a dataframe"
  )
  
  if (is.null(draw_ids)) {
    n_iter <- brms::ndraws(flocker_fit)
  } else {
    n_iter <- length(draw_ids)
  }
  
  if(history_condition) {
    use_components <- c("occ", "det", "col", "ex", "auto", "Omega")
    
    if(lik_type != "single_C"){
      if(is.null(new_data)){
        gp <- get_positions(flocker_fit)
        obs <- new_array(gp, flocker_fit$data$ff_y[gp])
      } else {
        gp <- get_positions(new_data)
        obs <- new_array(gp, flocker_fit$data$ff_y[gp])
      }
    } else { # lik_type is single_C
      if(is.null(new_data)){
        obs <- flocker_fit$data[c("ff_n_suc", "ff_n_trial")]
      } else {
        obs <- new_data[c("ff_n_suc", "ff_n_trial")]
      }
    }
  } else { # history_condition is FALSE
    use_components <- c("occ", "col", "ex", "auto", "Omega")
    obs <- NULL
  }
  
  if (lik_type %in% c("single")) {
    lps <- fitted_flocker(
      flocker_fit, components = use_components, draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, response = FALSE, unit_level = FALSE
    )
    Z <- get_Z_single(lps, sample, history_condition, obs)
  } else if (lik_type == "single_C") {
    lps <- fitted_flocker(
      flocker_fit, components = use_components, draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, response = FALSE, unit_level = FALSE
    )
    Z <- get_Z_single_C(lps, sample, history_condition, obs)
  } else if (lik_type == "augmented") {
    lps <- fitted_flocker(
      flocker_fit,
      components = use_components, 
      draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, response = FALSE, unit_level = FALSE
    )
    Z <- get_Z_augmented(lps, sample, history_condition, obs)
  } else if (lik_type %in% c("multi_colex")) {
    lps2 <- fitted_flocker(
      flocker_fit, components = c("occ", "colo", "ex"),
      draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, response = TRUE, unit_level = TRUE
    )
    init <- lps2$linpred_occ[,1,]
    colo <- lps2$linpred_col
    ex <- lps2$linpred_ex
    if(history_condition){
      lps1 <- fitted_flocker(
        flocker_fit, components = c("det"),
        draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels, response = TRUE, unit_level = FALSE
      )
      det <- lps1$linpred_det
    } else {
      det <- NULL
    }
    Z <- get_Z_dynamic(init, colo, ex, history_condition, sample, obs, det)
  } else if (lik_type %in% c("multi_colex_eq")) {
    if(history_condition){
      lps1 <- fitted_flocker(
        flocker_fit, components = c("det"),
        draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels, response = TRUE, unit_level = FALSE
      )
      det <- lps1$linpred_det
    } else {
      det <- NULL
    }
    lps2 <- fitted_flocker(
      flocker_fit, components = c("col", "ex"),
      draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, response = TRUE, unit_level = TRUE
    )
    colo <- lps2$linpred_col
    ex <- lps2$linpred_ex
    init <- colo[,1,] / (colo[,1,] + ex[,1,])
    Z <- get_Z_dynamic(init, colo, ex, history_condition, sample, obs, det)
  } else if (lik_type == "multi_autologistic") {
    if(history_condition){
      lps1 <- fitted_flocker(
        flocker_fit, components = c("det"),
        draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels, response = TRUE, unit_level = FALSE
      )
      det <- lps1$linpred_det
    } else {
      det <- NULL
    }
    lps2 <- fitted_flocker(
      flocker_fit, components = c("occ", "colo", "auto"),
      draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, response = TRUE, unit_level = TRUE
    )
    init <- lps2$linpred_occ[,1,]
    colo <- lps2$linpred_col
    ex <- 1 - boot::inv.logit(boot::logit(colo) + lps2$linpred_auto)
    Z <- get_Z_dynamic(init, colo, ex, history_condition, sample, obs, lps1$linpred_det)
  } else if (lik_type == "multi_autologistic_eq") {
    if(history_condition){
      lps1 <- fitted_flocker(
        flocker_fit, components = c("det"),
        draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels, response = TRUE, unit_level = FALSE
      )
      det <- lps1$linpred_det
    } else {
      det <- NULL
    }
    lps2 <- fitted_flocker(
      flocker_fit, components = c("col", "auto"),
      draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
      sample_new_levels = sample_new_levels, response = TRUE, unit_level = TRUE
    )
    colo <- lps2$linpred_col
    ex <- 1 - boot::inv.logit(boot::logit(colo) + lps2$linpred_auto)
    init <- colo[,1,] / (colo[,1,] + ex[,1,])
    Z <- get_Z_dynamic(init, colo, ex, history_condition, sample, obs, lps1$linpred_det)
  }
  class(Z) <- c("postZ", class(Z))
  Z
}

#' get Z matrix for single-season model
#' @param lps the linear predictors from the model
#' @param sample logical: return fitted probabilities or bernoulli samples
#' @param history_condition logical: condition on the observed history?
#' @param obs if history_condition is true, the observed histories
#' @return a matrix of fitted Z probabilities or sampled Z values. Rows are
#'   units and columns are posterior iterations.
#' @noRd
get_Z_single <- function(lps, sample, history_condition, obs = NULL){
  if(length(dim(lps$linpred_occ)) == 3) { # from flockerdata
    lpo <- lps$linpred_occ[ , 1, ] # first index is unit, second is visit, third is draw
  } else { # from data.frame
    assertthat::assert_that(length(dim(lps$linpred_occ)) == 2)
    lpo <- lps$linpred_occ
  }
  
  n_unit <- nrow(lpo)
  psi_all <- boot::inv.logit(lpo)
  
  if (!history_condition){
    if(sample) {
      Z <- new_matrix(psi_all, stats::rbinom(length(psi_all), 1, psi_all))
    } else {
      Z <- psi_all
    }
  } else {
    theta_all <- boot::inv.logit(lps$linpred_det) # lpo must be from flockerdata, since history_condition is TRUE
    
    # get emission likelihoods
    el_0 <- el_1 <- new_matrix(psi_all)
    
    for(i in seq_len(ncol(psi_all))){
      el_0[ , i] <- emission_likelihood(0, obs, theta_all[,,i])
      el_1[ , i] <- emission_likelihood(1, obs, theta_all[,,i])
    }
    
    hc <- Z_from_emission(el_0, el_1, psi_all)
    
    if(sample) {
      Z <- new_matrix(psi_all, stats::rbinom(length(hc), 1, hc))
    } else {
      Z <- hc
    }
  }
  Z
}

#' get Z matrix for single-season rep-constant model
#' @param lps the linear predictors from the model
#' @param sample logical: return fitted probabilities or bernoulli samples
#' @param history_condition logical: condition on the observed history?
#' @param obs if history_condition is true, the observed histories
#' @return a matrix of fitted Z probabilities or sampled Z values. Rows are
#'   units and columns are posterior iterations.
#' @noRd
get_Z_single_C <- function(lps, sample, history_condition, obs = NULL){
  if(length(dim(lps$linpred_occ)) == 3) { # from flockerdata
    lpo <- lps$linpred_occ[ , 1, ] # first index is unit, second is visit, third is draw
  } else { # from data.frame
    assertthat::assert_that(length(dim(lps$linpred_occ)) == 2)
    lpo <- lps$linpred_occ
  }
  
  n_unit <- nrow(lpo)
  psi_all <- boot::inv.logit(lpo)
  
  if (!history_condition){
    if(sample) {
      Z <- new_matrix(psi_all, stats::rbinom(length(psi_all), 1, psi_all))
    } else {
      Z <- psi_all
    }
  } else {
    theta_all <- boot::inv.logit(lps$linpred_det[ , 1, ])
    
    # get emission likelihoods
    el_0 <- el_1 <- new_matrix(psi_all)
    
    theta_c_all <- 1 - theta_all
    n_fail <- obs$ff_n_trial - obs$ff_n_suc
    for(i in seq_len(ncol(psi_all))){
      el_0[ , i] <- as.numeric(obs$ff_n_suc == 0)
      # we don't need to worry about the binomial coefficient here, because
      # the only time we care about the actual emission likelihood for a 1 is
      # when the history is all zeros (otherwise we care only that the
      # likelihood is nonzero)
      el_1[ , i] <- theta_all[,i]^obs$ff_n_suc * theta_c_all[,i]^n_fail
    }
    
    hc <- Z_from_emission(el_0, el_1, psi_all)
    
    if(sample) {
      Z <- new_matrix(psi_all, stats::rbinom(length(hc), 1, hc))
    } else {
      Z <- hc
    }
  }
  Z
}

#' get Z matrix for data-augmented model
#' @param lps the linear predictors from the model
#' @param sample logical: return fitted probabilities or bernoulli samples
#' @param history_condition logical: condition on the observed history?
#' @param obs if history_condition is true, the observed histories
#' @return an array of fitted Z probabilities or sampled Z values. Rows are
#'   units and columns are posterior iterations.
#' @noRd
get_Z_augmented <- function(lps, sample, history_condition, obs = NULL){
  lpo <- lps$linpred_occ[ , 1, , ] # first index is point, second is visit, third is species, fourth is draw
  n_point <- nrow(lpo)
  n_species <- ncol(lpo)
  n_unit <- n_point * n_species
  psi_all <- boot::inv.logit(lpo)
  Omega <- boot::inv.logit(lps$linpred_Omega[1,1,,])
  
  if (!history_condition){
    if(sample) {
      Z1 <- new_array(psi_all, stats::rbinom(length(psi_all), 1, psi_all))
      Z2 <- new_matrix(Omega, stats::rbinom(length(Omega), 1, Omega))
    } else {
      Z1 <- psi_all
      Z2 <- Omega
    }
  } else {
    theta_all <- boot::inv.logit(lps$linpred_det)
    
    # get emission likelihoods
    el_0 <- el_1 <- new_array(psi_all)
    
    message("computing emission probabilities")
    pb <- utils::txtProgressBar(max = ncol(psi_all))
    for(j in seq_len(ncol(psi_all))){
      utils::setTxtProgressBar(pb, j)
      for(i in seq_len(nslice(psi_all))){
        el_0[ , j, i] <- emission_likelihood(0, obs[,,j], theta_all[,,j,i])
        el_1[ , j, i] <- emission_likelihood(1, obs[,,j], theta_all[,,j,i])
      }
    }
    # history-conditioned probabilities, given species is in metacommunity
    hc <- Z_from_emission(el_0, el_1, psi_all)
    if(sample) {
      Z1 <- new_array(psi_all, stats::rbinom(length(hc), 1, hc))
    } else {
      Z1 <- hc
    }
    # history-conditioned Omegas
    log_lik_absent <- apply(log(el_0), c(2,3), sum)
    log_lik_present <- apply(log(el_1), c(2,3), sum)
    hc2 <- exp(log_lik_present - apply(abind::abind(log_lik_absent, log_lik_present, along = 3), c(1,2), matrixStats::logSumExp))
    if(sample){
      Z2 <- new_matrix(hc2, stats::rbinom(length(hc2), 1, hc2))
    } else {
      Z2 <- hc2
    }
  }
  Z <- new_array(Z1)
  for(i in seq_len(nrow(Z))){
    Z[i,,] <- Z1[i,,] * Z2
  }
  Z
}

#' get Z matrix for dynamic model
#' @param init matrix of initial occupancy probabilities (rows are sites and 
#'   columns are draws)
#' @param colo array of colonization probabilities (rows are sites, columns are
#'  timesteps, slices are draws)
#' @param ex array of extinction probabilities (rows are sites, columns are
#'  timesteps, slices are draws)
#' @param history_condition should we condition on the observed history
#' @param sample Should the return be posterior probabilities of occupancy (FALSE),
#'  valid posterior predictive samples (TRUE). 
#' @param obs Array of observations. Must be NULL if history_condition
#'   is FALSE. Rows are sites, columns are visits, and slices are timesteps.
#' @param det Array of detection probabilities. Ignored (and defaults to NULL) if
#'   history_condition is FALSE. Rows are sites, columns are visits, slices
#'   are timesteps, and slice_2s are draws.
#' @return an array of occupancy probabilities or Z samples. Rows are sites, 
#'   columns are timesteps, and slices are draws.
#' @details History-conditioning is via the forward-backward algorithm (forward-
#'  filtering-backward-sampling when sample is TRUE)
#' @noRd
get_Z_dynamic <- function(
    init, colo, ex, history_condition, sample, obs = NULL, det = NULL
    ){
  assertthat::assert_that(is_one_logical(history_condition))
  assertthat::assert_that(is.matrix(init))
  assertthat::assert_that(is.array(colo))
  assertthat::assert_that(is.array(ex))
  
  nsite <- nrow(init)
  ntimestep <- ncol(colo)
  ndraw <- ncol(init)
  
  if(history_condition){
    assertthat::assert_that(identical(dim(obs), dim(det)[1:3]))
    assertthat::assert_that(isTRUE(dim(det)[4] == ndraw))
  } else {
    assertthat::assert_that(is.null(obs))
    det <- NULL
  }
  assertthat::assert_that(identical(dim(colo), dim(ex)))
  assertthat::assert_that(isTRUE(nrow(colo) == nsite))
  assertthat::assert_that(isTRUE(nslice(colo) == ndraw))
  
  out <- new_array(colo)
  
  if(history_condition){
    assertthat::assert_that(isTRUE(nrow(obs) == nrow(colo)))
    assertthat::assert_that(isTRUE(nslice(obs) == ncol(colo)))
    for(i in 1:nsite){
      for(j in 1:ndraw){
        el0 <- emission_likelihood(0, t(obs[i,,]), t(det[i,,,j]))
        el1 <- emission_likelihood(1, t(obs[i,,]), t(det[i,,,j]))
        if(sample){
          out[i,,j] <- forward_backward_sampling(el0, el1, init[i,j], colo[i,,j], ex[i,,j])
        } else {
          out[i,,j] <- forward_backward_algorithm(el0, el1, init[i,j], colo[i,,j], ex[i,,j])
        }
      }
    }
  } else {
    assertthat::assert_that(is.null(obs), msg = "obs should be NULL if history_condition is false")
    assertthat::assert_that(is.null(det), msg = "det should be NULL if history_condition is false")
    for(i in 1:nsite){
      for(j in 1:ndraw){
        out[i,,j] <- forward_sim(init[i,j], colo[i,,j], ex[i,,j], sample)
      }
    }
  }
  out
}

#' Forward simulation for a dynamic model
#' @param init initial occupancy probability
#' @param colo vector of colonization probabilities
#' @param ex vector of extinction probabilities
#' @param sample Should the return be posterior probabilities of occupancy (FALSE),
#'     valid posterior predictive samples (TRUE)
#' @return vector of occupancy probabilities or posterior predictive samples
#' @noRd
forward_sim <- function(init, colo, ex, sample = FALSE){
  assertthat::assert_that(identical(is.na(colo), is.na(ex)))
  assertthat::assert_that(identical(length(colo), length(ex)))
  length_out <- length(colo)
  if(!identical(colo, NA)){
    assertthat::assert_that(!all(is.na(colo)))
    colo <- colo[1:max(which(!is.na(colo)))]
    ex <- ex[1:max(which(!is.na(colo)))]
    assertthat::assert_that(!any(is.na(colo)))
  }
  assertthat::assert_that(is.numeric(init))
  assertthat::assert_that(length(init) == 1)
  assertthat::assert_that(init >= 0 & init <= 1)
  assertthat::assert_that(all(colo >= 0, na.rm = TRUE))
  assertthat::assert_that(all(colo <= 1, na.rm = TRUE))
  assertthat::assert_that(all(ex >= 0, na.rm = TRUE))
  assertthat::assert_that(all(ex <= 1, na.rm = TRUE))
  
  out <- rep(NA, length_out)
  if(sample){
    out[1] <- stats::rbinom(1, 1, init)
    if(length(colo) > 1){
      for(i in 2:length(colo)){
        if(out[i - 1] == 0){
          out[i] <- stats::rbinom(1, 1, colo[i])
        } else {
          out[i] <- stats::rbinom(1, 1, 1 - ex[i])
        }
      }
    }
  } else { # not sampling
    out[1] <- init
    if(length(colo) > 1){
      for(i in 2:length(colo)){
        out[i] <- (1 - out[i - 1]) * colo[i] +
          (out[i - 1]) * (1 - ex[i])
      }
    }
  }
  out
}

#' forward backward algorithm
#' @param el0 the emission probabilities given non-occupancy
#' @param el1 the emission probabilities given occupancy
#' @param init the initial occupancy probability
#' @param colo a vector giving the colonization probs
#' @param ex a vector giving the extinction probs
#' @return a vector of posterior Z probabilities
#' @noRd
forward_backward_algorithm <- function(el0, el1, init, colo, ex) {
  # Run the forward algorithm to get forward probabilities
  forward_probs <- forward_algorithm(el0, el1, init, colo, ex)
  
  # Run the backward algorithm to get backward probabilities
  backward_probs <- backward_algorithm(el0, el1, colo, ex)
  
  # Compute the smoothed state probabilities
  smoothed_probs <- forward_probs * backward_probs
  
  # Normalize the smoothed probabilities of occupancy
  normalized_probs <- apply(smoothed_probs, 1, function(x) x[2] / sum(x))
  
  normalized_probs
}

#' Forward filtering backward sampling algorithm
#' @inheritParams forward_backward_algorithm
#' @noRd
forward_backward_sampling <- function(el0, el1, init, colo, ex) {
  # Run the forward algorithm to get forward probabilities
  forward_probs <- forward_algorithm(el0, el1, init, colo, ex)
  
  # Initialize the sampled state sequence
  n_tot <- n <- length(el0)
  sampled_states <- rep(NA, n_tot)
  
  while(TRUE){
    probs <- forward_probs[n, ]
    if(any(is.na(probs))){
      assertthat::assert_that(all(is.na(probs)))
      n <- n - 1
    } else {
      break
    }
  }

  
  # Sample the last state based on the forward probabilities
  sampled_states[n] <- sample(c(0, 1), size = 1, prob = probs)
    # we use sample rather than stats::rbinom to avoid the need to normalize the probs
  
  if(n > 1){
    # get transition probabilities matrix
    trans <- matrix(c(1 - colo, colo, ex, 1 - ex), nrow = n_tot)[2:n, ]
    
    # Iterate over the timesteps in reverse, starting from the second-to-last timestep
    for (t in (n - 1):1) {
      # Compute the conditional probabilities for the previous state, given the sampled state at t+1
      cond_probs <- trans[t, (2 * sampled_states[t + 1] + 1):(2 * sampled_states[t + 1] + 2)] * forward_probs[t,]
      
      # Sample the current state based on the conditional probabilities
      sampled_states[t] <- sample(c(0, 1), size = 1, prob = cond_probs)
    }
  }
  sampled_states
}

#' Forward algorithm
#' @inheritParams forward_backward_algorithm
#' @return an n x 2 matrix giving the forward probabilities associated with 
#'   non-occupancy (first column) and occupancy (second column)
#' @noRd
forward_algorithm <- function(el0, el1, init, colo, ex) {
  n <- length(el0)
  assertthat::assert_that(length(el1) == n)
  assertthat::assert_that(length(colo) == n)
  assertthat::assert_that(length(ex) == n)
  assertthat::assert_that(length(init) == 1)
  init <- c(1-init, init)
  
  # Initialize the forward probabilities matrix
  forward_probs <- matrix(NA, nrow = n, ncol = 2)
  
  # Compute the forward probabilities for the first timestep
  forward_probs[1,] <- init * c(el0[1], el1[1])
  
  if(n > 1){
    # get transition probabilties matrix
    trans <- matrix(c(1 - colo, colo, ex, 1 - ex), nrow = n)[2:n, ]
  
    # Iterate over the remaining timesteps
    for (t in 2:n) {
      # Extract the transition probabilities for the current timestep
      trans_t <- matrix(c(trans[t - 1, 1], trans[t - 1, 2], trans[t - 1, 3], trans[t - 1, 4]), nrow = 2, byrow = TRUE)
      
      # Compute the forward probabilities for the current timestep
      forward_probs[t,] <- (forward_probs[t - 1,] %*% trans_t) * c(el0[t], el1[t])
    }
  }
  forward_probs
}

#' Backward algorithm
#' @inheritParams forward_backward_algorithm
#' @return an n x 2 matrix giving the backward probabilities
#' @noRd
backward_algorithm <- function(el0, el1, colo, ex) {
  n <- length(el0)
  assertthat::assert_that(length(el1) == n)
  assertthat::assert_that(length(colo) == n)
  assertthat::assert_that(length(ex) == n)
  assertthat::assert_that(n > 1)
  
  # get transition probabilties matrix
  trans <- matrix(c(1 - colo, colo, ex, 1 - ex), nrow = n)[2:n, ]
  
  # Initialize the backward probabilities matrix
  backward_probs <- matrix(NA, nrow = n, ncol = 2)
  
  # Initialize the backward probabilities for the last timestep
  backward_probs[n,] <- c(1, 1)
  
  # Iterate over the timesteps in reverse
  for (t in (n - 1):1) {
    # Extract the transition probabilities for the current timestep
    trans_t <- matrix(c(trans[t, 1], trans[t, 2], trans[t, 3], trans[t, 4]), nrow = 2, byrow = TRUE)
    
    # Compute the backward probabilities for the current timestep
    if(any(is.na(trans_t))){
      assertthat::assert_that(all(is.na(trans_t)))
      backward_probs[t,] <- c(1,1)
    } else {
      backward_probs[t,] <- (trans_t %*% c(el0[t + 1], el1[t + 1])) * backward_probs[t + 1,]
    }
  }
  backward_probs
}
