#' Get posterior distribution of Z matrix
#' @param flocker_fit A flocker_fit object
#' @param n_iter The number of posterior iterations desired. If `NULL`, use
#'     all available posterior iterations.
#' @param history_condition Should the posterior distribution for Z directly 
#'     condition on the observed detection history (`TRUE`) or not (`FALSE`)?
#'     For example, at sites with at least one detection, the true occupancy 
#'     state conditioned on the history is one with absolute certainty. Without 
#'     directly conditioning on the history, the occupancy state is controlled 
#'     by the posterior distribution for the occupancy probability psi.
#' @param new_data Optional new data at which to predict the Z matrix. Generally
#'     the output of `make_flocker_data`. Cannot be used if `history_condition` is
#'     set to `TRUE`.
#' @param sample_new_levels If new_data is provided and contains random effect
#'     levels not present in the original data, how should predictions be
#'     handled? Passed directly to brms::prepare_predictions, which see. 
#' @return The posterior Z matrix. Rows are iterations and columns
#'     are closure-units. Values are samples from the posterior distribution
#'     of occupancy probability.
#' @export

get_Z <- function (flocker_fit, n_iter = NULL, history_condition = TRUE, 
                   new_data = NULL, sample_new_levels = "uncertainty") {
  # Input checking and processing
  assertthat::assert_that(
    is_flocker_fit(flocker_fit),
    msg = "flocker_fit must be an object of class `flocker_fit`"
  )
  total_iter <- brms::niterations(flocker_fit)*brms::nchains(flocker_fit)
  if (is.null(n_iter)) {
    n_iter <- total_iter
  }
  n_iter <- as.integer(n_iter)
  assertthat::assert_that(
    is_one_pos_int(n_iter),
    msg = "n_iter, if supplied, must be a single positive integer"
  )
  assertthat::assert_that(
    n_iter <= total_iter,
    msg = "requested number of iterations exceeds iterations contained in flocker_fit"
  )
  assertthat::assert_that(
    is_one_logical(history_condition),
    msg = "history_condition must be a single logical value"
  )
  if (is.null(new_data)) {
    new_data <- flocker_fit$data
  } else {
    assertthat::assert_that(
      "flocker_data" %in% class(new_data),
      msg = "new_data is provided but is not a `flocker_data` object"
    )
    assertthat::assert_that(
      new_data$type == attributes(flocker_fit)$data_type & 
        new_data$fp == attributes(flocker_data)$fp,
      msg = "new_data is formatted for a different model type than flocker_fit"
    )
    new_data <- new_data$data
  }
  
  message("computing linear predictors")
  
  lps <- fitted_flocker(
    flocker_fit, ndraws = n_iter, new_data = new_data, allow_new_levels = TRUE, 
    sample_new_levels = sample_new_levels, response = F
    )
  lik_type <- type_flocker_fit(flocker_fit)
  
  if (lik_type == "single") {
    n_unit <- new_data$n_unit[1]
    lpo <- lps$linpred_occ[1:n_unit, ]
    psi_all <- boot::inv.logit(lpo)
    Z <- matrix(data = 1, nrow = n_iter, ncol = n_unit)
    lpd <- lps$linpred_det
    
    
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
