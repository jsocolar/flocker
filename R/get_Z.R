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
#' @param new_data Optional new data at which to predict the Z matrix. Can be 
#'     the output of `make_flocker_data` or the `unit_covs` input to 
#'     `make_flocker_data` provided that `history_condition` is `FALSE`
#' @param allow_new_levels allow new levels for random effect terms in `new_data`?
#'     Will error if set to `FALSE` and new levels are provided in `new_data`.
#' @param sample_new_levels If `new_data` is provided and contains random effect
#'     levels not present in the original data, how should predictions be
#'     handled? Passed directly to `brms::prepare_predictions`, which see. 
#' @return The posterior Z matrix in the shape of the first slice of `obs` as
#'     passed to make_flocker_data, with posterior iterations stacked along the
#'     final dimension
#' @export

get_Z <- function (flocker_fit, draw_ids = NULL, history_condition = TRUE, 
                   new_data = NULL, allow_new_levels = FALSE, 
                   sample_new_levels = "uncertainty") {
  # Input checking and processing
  
  assertthat::assert_that(
    is_one_logical(history_condition),
    msg = "history_condition must be a single logical value"
  )
  assertthat::assert_that(
    !history_condition | is_flocker_data(new_data) | is.null(new_data),
    msg = "conditioning on observed detection histories while using the `new_data` argument requires passing a `flocker_data` object"
  )
  
  if(is.null(ndraws)){
    ndraws <- brms::niterations(flocker_fit)
  }
  
  lps <- fitted_flocker(
    flocker_fit, draw_ids = draw_ids, new_data = new_data, allow_new_levels = allow_new_levels, 
    sample_new_levels = sample_new_levels, response = FALSE, unit_level = FALSE
  )
  
  lik_type <- type_flocker_fit(flocker_fit)
  
  if (lik_type == "single") {
    lpo <- lps$linpred_occ[1 , ]
    psi_all <- boot::inv.logit(lpo)
    Z <- matrix(data = NA, nrow = n_unit, ncol = n_iter)
    lpd <- lps$linpred_det
    
    if(!hist_condition){
      
    }
  } else if (lik_type == "single_C") {
    
  } else if (lik_type == "augmented") {
    
  } else if (lik_type == "multi_colex") {
    
  } else if (lik_type == "multi_colex_eq") {
    
  } else if (lik_type == "multi_autologistic") {
    
  } else if (lik_type == "multi_autologistic_eq") {
    
  } else if (lik_type == "single_fp") {
    
  } else if (lik_type == "multi_colex_fp") {
    
  } else if (lik_type == "multi_colex_eq_fp") {
    
  }
  
  class(Z) <- c("postZ", "matrix")
  Z
  
  
  
  if (lik_type == "V") {
    
  } else {
    n_unit <- nrow(new_data)
  }

  

  
  message("computing Z")
  if (!history_condition) {
    Z <- psi_all
  } else {
    pb <- utils::txtProgressBar(min = 0, max = n_unit, style = 3)
    antitheta_all <- boot::inv.logit(-lpd)
    if (lik_type == "V") {
      index_matrix <- new_data[grep("^rep_index", names(new_data))] 
      Q <- as.integer(new_data$Q[1:n_unit]) # get Q
      for (i in 1:n_unit) {
        if (Q[i] == 0L) { # if Q is 0
          psi <- psi_all[ , i]
          antitheta <- antitheta_all[ , as.integer(index_matrix[1:new_data$n_rep[i], i])]
          numerator <- psi * matrixStats::rowProds(antitheta)
          denominator <- numerator + (1 - psi)
          Z[, i] <- numerator/denominator
        }
        utils::setTxtProgressBar(pb, i)
      }
    } else { # lik_type == "C"
      Q <- as.integer(new_data$n_suc > 0)
      Z <- matrix(data = 1, nrow = n_iter, ncol = n_unit)
      for (i in 1:n_unit) {
        if (Q[i] == 0L) {
          psi <- psi_all[, i]
          antitheta <- antitheta_all[, i]
          numerator <- psi * antitheta^new_data$n_trial[i] 
          denominator <- numerator + (1 - psi) 
          Z[, i] <- numerator/denominator
        }
        utils::setTxtProgressBar(pb, i)
      }
    }
    close(pb)
  }

  class(Z) <- c("postZ", "matrix")
  return(Z)
}
