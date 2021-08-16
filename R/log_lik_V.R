#' Compute sitewise log-likelihood matrix for a flocker_fit object with 
#' visit-variable covariates
#' @param flocker_fit_V A flocker_fit with visit-variable covariates
#' @return A sitewise posterior log-likelihood matrix
#' @export

log_lik_V <- function(flocker_fit_V) {
  if (!("flocker_fit" %in% class(flocker_fit_V))) {
    stop("flocker_fit_V must be an object of class flocker_fit.")
  }
  if (attributes(flocker_fit_V)$lik_type != "V") {
    stop("flocker_fit_V works only for flocker_fits with visit-variable covariates")
  }
  
  # dimensions
  n_site <- flocker_fit_V$data$nsite[1]
  n_visit <- max(flocker_fit_V$data$nvisit)
  n_iter <- dim(flocker_fit_V$fit)[1]
  
  visit_index_matrix <- 
    as.matrix(flocker_fit_V$data[1:n_site, grepl("visit_index", names(flocker_fit_V$data))])
  
  lpo_t <- t(brms::posterior_linpred(flocker_fit_V, dpar = "occ"))
  lpd_t <- t(brms::posterior_linpred(flocker_fit_V, dpar = "mu"))
  
  # create long-format dataframe (with iterations stacked down rows)
  # note: missed site visits are inserted as -99s
  all_iters <- data.frame(resp = -99,
                          site_index = rep(1:n_site, n_visit), 
                          visit_index = c(visit_index_matrix),
                          Q = rep(flocker_fit_V$data$Q[1:n_site], n_visit), 
                          # note: everything above this is getting duplicated n_iter times
                          iter = rep(1:n_iter, each = n_site*n_visit), 
                          lpo = NA, 
                          lpd = NA)
  all_iters$resp[all_iters$visit_index != -99] <- flocker_fit_V$data$y
  all_iters$lpo[all_iters$visit_index != -99] <- c(lpo_t)
  all_iters$lpd[all_iters$visit_index != -99] <- c(lpd_t)
  
  # calculate visit-level component of likelihood
  all_iters$ll <- calc_log_lik_partial(all_iters$resp, all_iters$Q, all_iters$lpd)
  
  # spread this to wide format (i.e. 1 col per visit)
  visit_index <- rep(rep(1:n_visit, each=n_site), n_iter)
  
  ll_partial_V <- do.call("cbind", 
                          lapply(1:n_visit, function(x) matrix(all_iters$ll[visit_index == x])))
  
  ll_partial_S <- data.frame(Q = rep(all_iters$Q[1:n_site], n_iter),
                             lpo = all_iters$lpo[visit_index == 1], # note: duplicated across visits 
                             iter = all_iters$iter[visit_index == 1]) 
  
  # finish likelihood calculation
  Q_index <- as.logical(ll_partial_S$Q)
  ll_partial_S$log_lik <- NA
  ll_partial_S$log_lik[Q_index] <- log_inv_logit(ll_partial_S$lpo[Q_index]) + 
    rowSums(ll_partial_V[Q_index,])
  ll_partial_S$log_lik[!Q_index] <- matrixStats::rowLogSumExps(
    cbind(log1m_inv_logit(ll_partial_S$lpo[!Q_index]),
          log_inv_logit(ll_partial_S$lpo[!Q_index]) + rowSums(ll_partial_V[!Q_index,])))
  
  # unstack to matrix [n_iter, n_site]
  log_lik_mat <- t(unstack(ll_partial_S[c("log_lik", "iter")], log_lik ~ iter))
  
  return(log_lik_mat)
}

#' Compute the part of the log-likelihood relating to visits. To be used 
#' internally in log_lik_V(). Missed visits are returned as 0s
#' @param resp the response vector (detection/non-detection) at the site. Missing 
#' visits are represented as -99
#' @param Q whether there is at least one detection at a species:point combination
#' @param lpd the logit-scale linear predictor

calc_log_lik_partial <- function(resp, Q, lpd) {
  Q_index <- as.logical(Q)
  ll <- rep(NA, length(Q))
  
  ll[Q_index] <- ifelse(as.logical(resp[Q_index]),
                        log_inv_logit(lpd[Q_index]),
                        log1m_inv_logit(lpd[Q_index]))
  
  ll[!Q_index] <- log1m_inv_logit(lpd[!Q_index])
  
  ll[resp == -99] <- 0
  return(ll)
}

