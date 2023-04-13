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
    lpo <- lps$linpred_occ[ , 1, ] # first index is unit, second is visit, third is draw
    n_unit <- nrow(lpo)
    psi_all <- boot::inv.logit(lpo)
    Z <- matrix(nrow = n_unit, ncol = n_iter)
    
    assertthat::assert_that(identical(dim(psi_all), dim(Z)))
    
    if (!history_condition){
      if(sample) {
        Z <- matrix(stats::rbinom(length(psi_all), 1, psi_all), n_unit, n_iter)
      } else {
        Z <- psi_all
      }
    } else {
      theta_all <- boot::inv.logit(lps$linpred_det)
      
      # get emission likelihoods
      el_0 <- el_1 <- matrix(nrow = nrow(psi_all), ncol = ncol(psi_all))
      
      for(i in seq_len(ncol(psi_all))){
        el_0[ , i] <- mapply(function(a, b){emission_likelihood(0, a, b)}, asplit(obs, 1), asplit(theta_all[ , , i], 1))
        el_1[ , i] <- mapply(function(a, b){emission_likelihood(1, a, b)}, asplit(obs, 1), asplit(theta_all[ , , i], 1))
      }
      
      assertthat::assert_that(identical(dim(psi_all), dim(el_0)))
      
      hc <- (psi_all * el_1)/(psi_all * el_1 + (1 - psi_all) * el_0)
      
      if(sample) {
        Z <- matrix(stats::rbinom(length(hc), 1, hc), n_unit, n_iter)
      } else {
        Z <- hc
      }
    }
  } else if (lik_type == "single_C") {
    stop("get_Z not yet implmented for rep-constant models")
  } else if (lik_type == "augmented") {
    stop("get_Z not yet implmented for augmented models")
  } else if (lik_type == "multi_colex") {
    stop("get_Z not yet implmented for dynamic models")
  } else if (lik_type == "multi_colex_eq") {
    stop("get_Z not yet implmented for dynamic models")
  } else if (lik_type == "multi_autologistic") {
    stop("get_Z not yet implmented for dynamic models")
  } else if (lik_type == "multi_autologistic_eq") {
    stop("get_Z not yet implmented for dynamic models")
  } else if (lik_type == "single_fp") {
    stop("get_Z not yet implmented for false positive")
  } else if (lik_type == "multi_colex_fp") {
    stop("get_Z not yet implmented for dynamic models")
  } else if (lik_type == "multi_colex_eq_fp") {
    stop("get_Z not yet implmented for dynamic models")
  }
  
  class(Z) <- c("postZ", "matrix")
  Z
}
