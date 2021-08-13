#' Compute sitewise log-likelihood matrix for a flocker_fit object with 
#' visit-variable covariates
#' @param flocker_fit_V A flocker_fit with visit-variable covariates
#' @return A sitewise posterior log-likelihood matrix
#' @export

log_lik_V <- function (flocker_fit_V) {
  if (!("flocker_fit" %in% class(flocker_fit_V))) {
    stop("flocker_fit_V must be an object of class flocker_fit.")
  }
  if (attributes(flocker_fit_V)$lik_type != "V") {
    stop("flocker_fit_V works only for flocker_fits with visit-variable covariates")
  }
  n_site <- flocker_fit_V$data$nsite[1]
  visit_index_matrix <- 
    flocker_fit_V$data[1:n_site, grepl("visit_index", names(flocker_fit_V$data))]
  resp <- flocker_fit_V$data$y
  Q <- flocker_fit_V$data$Q[1:n_site]
  lpo <- brms::posterior_linpred(flocker_fit_V, dpar = "occ")[ , 1:n_site]
  lpd <- brms::posterior_linpred(flocker_fit_V, dpar = "mu")
  
  n_iter <- dim(flocker_fit_V$fit)[1]
  log_lik_matrix <- matrix(NA, nrow = n_iter, ncol = n_site)
  
  for (iter in 1:n_iter) {
    log_lik_matrix[iter, ] <-
      log_lik_V_iter(n_site, visit_index_matrix, resp, Q, lpo[iter, ], lpd[iter, ])
  }
  return(log_lik_matrix)
}


#' Compute sitewise log likelihood for a single posterior iteration for a 
#' flocker_fit object with visit-variable covariates
#' @param n_site number of sites or species-sites in the data
#' @param visit_index_matrix a matrix or dataframe giving the first n_site rows
#'     of the matrix of visit indices in the data (subsequent rows are redundant)
#' @param resp a vector giving the integer response data
#' @param Q a vector giving the first n_site elements of Q (whether there is at
#'     least one detection at the site)
#' @param lpo_iter a vector giving the first n_site elements of the logit-scale 
#'     linear predictor for occupancy for the iteration of interest
#' @param lpd_iter a vector giving the first n_site elements of the logit-scale 
#'     linear predictor for detection for the iteration of interest
#' @return a vector of sitewise log-likelihood values for the iteration of interest

log_lik_V_iter <- function (n_site, visit_index_matrix, resp, Q, 
                            lpo_iter, lpd_iter) {
  log_lik_iter <- rep(NA, n_site)
  for (i in 1:n_site) {
    visit_subset_i <- as.integer(visit_index_matrix[i, ])
    resp_i <- resp[visit_subset_i]
    Q_i <- Q[i]
    lpo_i <- lpo_iter[i]
    lpd_i <- lpd_iter[visit_subset_i]
    log_lik_iter[i] <-
      log_lik_V_site(resp_i, Q_i, lpo_i, lpd_i)
  }
  return(log_lik_iter)
}


#' Compute the likelihood for a single site in a single posterior iteration for a
#' flocker_fit object with visit-variable covariates
#' @param resp_i the response vector (detection/non-detection) at the site
#' @param Q_i whether there is at least one detection, equivalent to
#'     (\code{as.integer(sum(resp_i) > 1)})
#' @param lpo_i the logit-scale linear predictor for the site and iteration
#' @param lpd_i a vector (or scalar if the site has just one visit) giving the 
#'     logit-scale linear predictor for the site and iteration
log_lik_V_site <- function (resp_i, Q_i, lpo_i, lpd_i) {
  log_lik_occ <- log_inv_logit(lpo_i)
  if (Q_i == 1) {
    log_lik <- log_lik_occ
    for (v in 1:length(resp_i)) {
      log_lik <- log_lik + 
        as.numeric(resp_i[v] == 1) * log_inv_logit(lpd_i[v]) +
        as.numeric(resp_i[v] == 0) * log1m_inv_logit(lpd_i[v])
    }
  } else {
    log_lik_occpart <- log_lik_occ
    for (v in 1:length(resp_i)) {
      log_lik_occpart <- log_lik_occpart + log1m_inv_logit(lpd_i[v])
    }
    log_lik <- matrixStats::logSumExp(c(log1m_inv_logit(lpo_i), log_lik_occpart))
  }
  return(log_lik)
}

