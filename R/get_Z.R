#' Get posterior distribution of Z matrix
#' @param flocker_fit A flocker_fit object
#' @param n_iter The number of posterior iterations desired
#' @return The posterior Z matrix. Rows are iterations and columns
#'     are closure-units. Values are samples from the posterior distribution
#'     of occupancy probability.
#' @export

get_Z <- function (flocker_fit, n_iter = NULL) {
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
  
  iter <- (1:n_iter)*floor(total_iter/n_iter)
  
  lpo <- brms::posterior_linpred(flocker_fit, dpar = "occ", draw_ids = iter)
  lpd <- brms::posterior_linpred(flocker_fit, dpar = "mu", draw_ids = iter)
  
  lik_type <- attributes(flocker_fit)$lik_type
  
  if (lik_type == "V") {
    n_unit <- flocker_fit$data$n_unit[1]
    Q <- as.integer(flocker_fit$data$Q[1:n_unit])
    Z <- matrix(data = 1, nrow = length(iter), ncol = n_unit)
    psi_all <- boot::inv.logit(lpo)
    antitheta_all <- boot::inv.logit(-lpd)
    
    index_matrix <- flocker_fit$data[grep("^rep_index", names(flocker_fit$data))]
    
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
    n_unit <- nrow(flocker_fit$data)
    Q <- as.integer(flocker_fit$data$n_suc > 0)
    Z <- matrix(data = 1, nrow = length(iter), ncol = n_unit)
    psi_all <- boot::inv.logit(lpo)
    antitheta_all <- boot::inv.logit(-lpd)
    for (i in 1:n_unit) {
      if (Q[i] == 0L) {
        psi <- psi_all[, i]
        antitheta <- antitheta_all[, i]
        numerator <- psi * antitheta^flocker_fit$data$n_trial[i]
          denominator <- numerator + (1 - psi)
          Z[, i] <- numerator/denominator
      }
    }
  }
  class(Z) <- c("postZ", "matrix")
  return(Z)
}
