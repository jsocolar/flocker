#' Get posterior distribution of Z matrix
#' @param flocker_fit A flocker_fit object
#' @param n_iter The number of posterior iterations desired. If `NULL`, use
#'     all available posterior iterations.
#' @param hist_condition Should the posterior distribution for Z directly 
#'     condition on the observed detection history (`TRUE`) or not (`FALSE`)?
#'     For example, at sites with at least one detection, the true occupancy 
#'     state conditioned on the history is one with absolute certainty. Without 
#'     directly conditioning on the history, the occupancy state is controlled 
#'     by the posterior distribution for the occupancy probability psi. Of 
#'     course even without conditioning directly on the detection history, we 
#'     still condition indirectly on the observed history via the fitted value 
#'     of psi, which itself depends on all of the observed detection histories.
#' @param new_data Optional new data at which to predict the Z matrix. Generally
#'     the output of `make_flocker_data`. Cannot be used if `hist_condition` is
#'     set to `TRUE`.
#' @param sample_new_levels If new_data is provided and contains random effect
#'     levels not present in the original data, how should predictions be
#'     handled? Passed directly to brms::prepare_predictions, which see. 
#' @return The posterior Z matrix. Rows are iterations and columns
#'     are closure-units. Values are samples from the posterior distribution
#'     of occupancy probability.
#' @export

get_Z <- function (flocker_fit, n_iter = NULL, hist_condition = TRUE, 
                   new_data = NULL, sample_new_levels = "uncertainty") {
  if (!is_flocker_fit(flocker_fit)) {
    stop("flocker_fit must be an object of class `flocker_fit`")
  }
  
  total_iter <- brms::niterations(flocker_fit)*brms::nchains(flocker_fit)
  
  if (is.null(n_iter)) {
    n_iter <- total_iter
  }
  
  n_iter <- as.integer(n_iter)
  
  if ((length(n_iter) != 1) | is.na(n_iter)) {
    stop("iter, if supplied, must be a single integer")
  }
  
  if (n_iter > total_iter) {
    stop("requested more iterations than contained in flocker_fit")
  }
  
  if (!is.logical(hist_condition)) {
    stop("hist_condition must be logical")
  }
  
  if (is.null(new_data)) {
    new_data <- flocker_fit$data
  } else if (hist_condition) {
    stop("cannot condition on Q if new_data is provided")
  } else if (!("flocker_data" %in% class(new_data))) {
    stop("new_data, if provided, must be a `flocker_data` object as produced by 
         `make_flocker_data`")
  } else {
    new_data <- new_data$data
  }
  
  iter <- (1:n_iter)*floor(total_iter/n_iter)
  
  lpo <- brms::posterior_linpred(flocker_fit, dpar = "occ", draw_ids = iter, 
                                 newdata = new_data, allow_new_levels = TRUE,
                                 sample_new_levels = sample_new_levels)
  lpd <- brms::posterior_linpred(flocker_fit, dpar = "mu", draw_ids = iter,
                                 newdata = new_data, allow_new_levels = TRUE,
                                 sample_new_levels = sample_new_levels)

  lik_type <- attributes(flocker_fit)$lik_type
  
  if (lik_type == "V") {
    n_unit <- new_data$n_unit[1]
    lpo <- lpo[, 1:n_unit]
  } else {
    n_unit <- nrow(new_data)
  }
  
  psi_all <- boot::inv.logit(lpo)
  Z <- matrix(data = 1, nrow = length(iter), ncol = n_unit)
  
  if (!hist_condition) {
    Z <- psi_all
  } else {
    antitheta_all <- boot::inv.logit(-lpd)
    if (lik_type == "V") {
      index_matrix <- new_data[grep("^rep_index", names(new_data))]
      Q <- as.integer(new_data$Q[1:n_unit])
      for (i in 1:n_unit) {
        if (Q[i] == 0L) {
          psi <- psi_all[ , i]
          antitheta <- antitheta_all[ , as.integer(index_matrix[1:flocker_data$data$n_rep[i], i])]
          numerator <- psi * matrixStats::rowProds(antitheta)
          denominator <- numerator + (1 - psi)
          Z[, i] <- numerator/denominator
        }
      }
    } else {
      Q <- as.integer(new_data$n_suc > 0)
      Z <- matrix(data = 1, nrow = length(iter), ncol = n_unit)
      for (i in 1:n_unit) {
        if (Q[i] == 0L) {
          psi <- psi_all[, i]
          antitheta <- antitheta_all[, i]
          numerator <- psi * antitheta^new_data$n_trial[i]
          denominator <- numerator + (1 - psi)
          Z[, i] <- numerator/denominator
        }
      }
    }
  }

  class(Z) <- c("postZ", "matrix")
  return(Z)
}
