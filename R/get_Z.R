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
#'     `make_flocker_data` provided that `history_condition` is `FALSE`
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
#' \dontrun{
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
  
  if (is.null(draw_ids)) {
    n_iter <- brms::ndraws(flocker_fit)
  } else {
    n_iter <- length(draw_ids)
  }
  
  lps <- fitted_flocker(
    flocker_fit, draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
    sample_new_levels = sample_new_levels, response = FALSE, unit_level = FALSE
  )
  
  if(history_condition) {
    if(is.null(new_data)){
      gp <- get_positions(flocker_fit)
      obs <- matrix(flocker_fit$data$ff_y[gp], nrow = nrow(gp), ncol = ncol(gp))
    } else {
      gp <- get_positions(new_data)
      obs <- matrix(new_data$ff_y[gp], nrow = nrow(gp), ncol = ncol(gp))
    }
  }
  
  lik_type <- type_flocker_fit(flocker_fit)
  
  if (lik_type == "single") {
    Z <- get_Z_single(lps, sample, history_condition, obs)
  } else if (lik_type == "single_C") {
    Z <- get_Z_single(lps, sample, history_condition, obs)
  } else if (lik_type == "augmented") {
    stop("get_Z not yet implmented for augmented models")
  } else if (lik_type == "multi_colex") {
    init <- lps$occ
    Z <- get_Z_dynamic(init, lps$colo, lps$ex, history_condition, obs, lps$det)
  } else if (lik_type == "multi_colex_eq") {
    init <- lps$colo[,1,] / (lps$colo[,1,] + lps$ex[,1,])
    Z <- get_Z_dynamic(init, lps$colo, lps$ex, history_condition, obs, lps$det)
  } else if (lik_type == "multi_autologistic") {
    init <- lps$occ
    colo <- lps$colo
    ex <- 1 - boot::inv.logit(boot::logit(lps$colo) + lps$auto)
    Z <- get_Z_dynamic(init, colo, ex, history_condition, obs, lps$det)
  } else if (lik_type == "multi_autologistic_eq") {
    colo <- lps$colo
    ex <- 1 - boot::inv.logit(boot::logit(lps$colo) + lps$auto)
    init <- colo[,1,] / (colo[,1,] + ex[,1,])
    Z <- get_Z_dynamic(init, colo, ex, history_condition, obs, lps$det)
  } else if (lik_type == "single_fp") {
    Z <- get_Z_single(lps, sample, history_condition, obs)
  } else if (lik_type == "multi_colex_fp") {
    init <- lps$occ
    Z <- get_Z_dynamic(init, lps$colo, lps$ex, history_condition, obs, lps$det)
  } else if (lik_type == "multi_colex_eq_fp") {
    init <- lps$colo[,1,] / (lps$colo[,1,] + lps$ex[,1,])
    Z <- get_Z_dynamic(init, lps$colo, lps$ex, history_condition, obs, lps$det)
  }
  class(Z) <- c("postZ", "matrix")
  Z
}

#' get Z matrix for single-season model
#' @param lps the linear predictors from the model
#' @param sample logical: return fitted probabilities or bernoulli samples
#' @param history_condition logical: condition on the observed history?
#' @param obs if history_condition is true, the observed histories
#' @return a matrix of fitted Z probabilities or sampled Z values. Rows are
#'   units and columns are posterior iterations.
get_Z_single <- function(lps, sample, history_condition, obs = NULL){
  lpo <- lps$linpred_occ[ , 1, ] # first index is unit, second is visit, third is draw
  n_unit <- nrow(lpo)
  psi_all <- boot::inv.logit(lpo)
  
  if (!history_condition){
    if(sample) {
      Z <- new_matrix(psi_all, stats::rbinom(length(psi_all), 1, psi_all))
    } else {
      Z <- psi_all
    }
  } else {
    theta_all <- boot::inv.logit(lps$linpred_det)
    
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

# PLACEHOLDER: Z MATRIX FOR AUGMENTED MODEL


#' get Z matrix for dynamic model
#' @param init matrix of initial occupancy probabilities (rows are sites and 
#'   columns are draws)
#' @param colo array of colonization probabilities (rows are sites, columns are
#'  timesteps, slices are draws)
#' @param ex array of extinction probabilities (rows are sites, columns are
#'  timesteps, slices are draws)
#' @param history_condition should we condition on the observed history
#' @param obs Array of observations. Ignored (and defaults to NULL) if history_condition
#'   is FALSE. Rows are sites, columns are visits, and slices are timesteps.
#' @param det Array of detection probabilities. Ignored (and defaults to NULL) if
#'   history_condition is FALSE. Rows are sites, columns are visits, slices
#'   are timesteps, and slice_2s are draws.
#' @param sample Should the return be posterior probabilities of occupancy (FALSE),
#'  valid posterior predictive samples (TRUE). 
#' @return an array of occupancy probabilities or Z samples. Rows are sites, 
#'   columns are timesteps, and slices are draws.
#' @details History-conditioning is via the forward-backward algorithm (forward
#'  filter backward sample when sample is TRUE)
get_Z_dynamic <- function(
    init, colo, ex, history_condition, obs = NULL, det = NULL
    ){
  assertthat::assert_that(is_one_logical(history_condition))
  assertthat::assert_that(is.matrix(init))
  assertthat::assert_that(is.array(colo))
  assertthat::assert_that(is.array(ex))
  
  nsite <- nrow(init)
  ntimestep <- ncol(colo)
  ndraw <- ncol(init)
  
  assertthat::assert_that(identical(dim(obs), dim(det)[1:3]))
  assertthat::assert_that(identical(dim(colo), dim(ex)))
  assertthat::assert_that(isTRUE(nrow(colo) == nsite))
  assertthat::assert_that(isTRUE(nslice(colo)) == ndraw)
  assertthat::assert_that(isTRUE(dim(det)[4] == ndraw))
  
  out <- array(dim = dim(colo))
  
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
forward_sim <- function(init, colo, ex, sample = FALSE){
  assertthat::assert_that(identical(is.na(colo), is.na(ex)))
  assertthat::assert_that(identical(length(colo), length(ex)))
  assertthat::assert_that(identical(colo, NA) | !(NA %in% colo))
  assertthat::assert_that(is.numeric(init))
  assertthat::assert_that(lenth(init) == 1)
  assertthat::assert_that(init >= 0 & init <= 1)
  
  out <- rep(NA, length(colo) + 1)
  if(sample){
    out[1] <- rbinom(1, 1, init)
    if(length(colo) > 1){
      for(i in 2:length(colo)){
        if(out[i - 1] == 0){
          out[i] <- rbinom(1, 1, colo[i])
        } else {
          out[i] <- rbinom(1, 1, 1 - ex[i])
        }
      }
    }
  } else { # not sampling
    out[1] <- init
    if(length(colo) > 1){
      for(i in 2:length(colo)){
        out[i] <- (out[i - 1] == 0) * colo[i] +
          (out[i - 1] == 1) * (1 - ex[i])
      }
    }
  }
  out
}

#' forward backward algorithm
#' @param e0 the emission probabilities given non-occupancy
#' @param e1 the emission probabilities given occupancy
#' @param init the initial state probabilities
#' @param colo a vector giving the colonization probs
#' @param ex a vector giving the extinction probs
#' @return a vector of posterior Z probabilities
forward_backward_algorithm <- function(e0, e1, init, colo, ex) {
  # Run the forward algorithm to get forward probabilities
  forward_probs <- forward_algorithm(obs, init, trans)
  
  # Run the backward algorithm to get backward probabilities
  backward_probs <- backward_algorithm(obs, trans)
  
  # Compute the smoothed state probabilities
  smoothed_probs <- forward_probs * backward_probs
  
  # Normalize the smoothed probabilities of occupancy
  normalized_probs <- apply(smoothed_probs, 1, function(x) x[2] / sum(x))
  
  normalized_probs
}

#' Forward filtering backward sampling algorithm
#' @inheritParams forward_backward_algorithm
forward_backward_sampling <- function(e0, e1, init, colo, ex) {
  # Run the forward algorithm to get forward probabilities
  forward_probs <- forward_algorithm(e0, e1, init, colo, ex)
  
  # Initialize the sampled state sequence
  n <- length(e0)
  sampled_states <- numeric(n)
  
  # Sample the last state based on the forward probabilities
  sampled_states[n] <- sample(c(0, 1), size = 1, prob = forward_probs[n,])
    # we use sample rather than stats::rbinom to avoid the need to normalize the probs
  
  if(n > 1){
    # get transition probabilities matrix
    trans <- matrix(c(1 - colo, colo, ex, 1 - ex), nrow = n)[2:n, ]
    
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
forward_algorithm <- function(e0, e1, init, colo, ex) {
  n <- length(e0)
  assertthat::assert_that(length(e1) == n)
  assertthat::assert_that(length(colo) == n)
  assertthat::assert_that(length(ex) == n)
  assertthat::assert_that(length(init) == 2)
  
  # Initialize the forward probabilities matrix
  forward_probs <- matrix(NA, nrow = n, ncol = 2)
  
  # Compute the forward probabilities for the first timestep
  forward_probs[1,] <- init * c(e0[1], e1[1])
  
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
backward_algorithm <- function(e0, e1, colo, ex) {
  n <- length(e0)
  assertthat::assert_that(length(e1) == n)
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
    backward_probs[t,] <- (trans_t %*% c(el0[t + 1], el1[t + 1])) * backward_probs[t + 1,]
  }
  backward_probs
}
