#' Compute unitwise log-likelihood matrix for a flocker_fit object with 
#' visit-variable covariates
#' @param flocker_fit_V A flocker_fit with visit-variable covariates
#' @return A unitwise posterior log-likelihood matrix
#' @export

log_lik_V <- function(flocker_fit_V) {
  if (!("flocker_fit" %in% class(flocker_fit_V))) {
    stop("flocker_fit_V must be an object of class flocker_fit.")
  }
  if (attributes(flocker_fit_V)$lik_type != "V") {
    stop("flocker_fit_V works only for flocker_fits with visit-variable covariates")
  }
  
  # dimensions
  n_unit <- flocker_fit_V$data$n_unit[1]
  max_visit <- max(flocker_fit_V$data$n_visit)
  n_iter <- prod(dim(flocker_fit_V$fit)[1:2])
  
  visit_index_matrix <- 
    as.matrix(flocker_fit_V$data[1:n_unit, grepl("visit_index", names(flocker_fit_V$data))])
  
  lpo_t <- t(brms::posterior_linpred(flocker_fit_V, dpar = "occ"))
  lpd_t <- t(brms::posterior_linpred(flocker_fit_V, dpar = "mu"))
  
  # create long-format dataframe (with iterations stacked down rows)
  # note: missed visits are inserted as -99s
  all_iters <- data.frame(resp = -99,
                          unit_index = rep(1:n_unit, max_visit), 
                          visit_index = c(visit_index_matrix),
                          Q = rep(flocker_fit_V$data$Q[1:n_unit], max_visit), 
                          # note: everything above this is getting duplicated n_iter times
                          iter = rep(1:n_iter, each = n_unit*max_visit), 
                          lpo = NA, 
                          lpd = NA)
  all_iters$resp[all_iters$visit_index != -99] <- flocker_fit_V$data$y
  all_iters$lpo[all_iters$visit_index != -99] <- c(lpo_t)
  all_iters$lpd[all_iters$visit_index != -99] <- c(lpd_t)
  
  # calculate visit-level component of likelihood
  all_iters$ll <- calc_log_lik_partial(all_iters$resp, all_iters$Q, all_iters$lpd)
  
  # spread this to wide format (i.e. 1 col per visit)
  visit_index <- rep(rep(1:max_visit, each=n_unit), n_iter)
  
  ll_partial_V <- do.call("cbind", 
                          lapply(1:max_visit, function(x) matrix(all_iters$ll[visit_index == x])))
  
  ll_partial_S <- data.frame(Q = rep(all_iters$Q[1:n_unit], n_iter),
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
  
  # unstack to matrix [n_iter, n_unit]
  log_lik_mat <- t(unstack(ll_partial_S[c("log_lik", "iter")], log_lik ~ iter))
  
  return(log_lik_mat)
}

#' Compute the part of the log-likelihood relating to visits. To be used 
#' internally in log_lik_V(). Missed visits are returned as 0s
#' @param resp the response vector (detection/non-detection) at the unit. Missing 
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




#' A log-likelihood function for the visit-constant occupancy model, sufficient for
#' \code{brms::loo(vc_fit)} to work. 
#' @param i Posterior iteration
#' @param prep Output of \code{brms::prepare_predictions}. See brms custom families
#' vignette at 
#' https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html
#' @return The log-likelihood for observation i

log_lik_occupancy_C <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  occ <- brms::get_dpar(prep, "occ", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  return(occupancy_C_lpmf(y, mu, occ, trials))
}

#' An R implementation of the visit constant lpmf without the binomial coefficient
#' @param y number of detections
#' @param mu logit-scale detection probability
#' @param occ logit-scale occupancy probability
#' @param trials number of visits
#' @return The log-likelihood

occupancy_C_lpmf <- Vectorize(
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
    return(out)
  }
)



#' An R implementation of occupancy_C_lpmf including the binomial coefficient.
#' Not currently in use.
#' @param y number of detections
#' @param mu logit-scale detection probability
#' @param occ logit-scale occupancy probability
#' @param trials number of visits
#' @return The log-likelihood

occupancy_C_lpmf_with_coef <- Vectorize(
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
    return(out)
  }
)


