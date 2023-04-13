#' Compute unit-wise or series-wise log-likelihood matrix for a flocker_fit object
#' @param flocker_fit A flocker_fit object
#' @param draw_ids the draw ids to compute log-likelihoods for. Defaults to using the full posterior. 
#' @return A unit-wise or series-wise posterior log-likelihood matrix, where iterations are rows and 
#'    units/series are columns
#' @export

log_lik_flocker <- function(flocker_fit, draw_ids = NULL) {
  assertthat::assert_that(is_flocker_fit(flocker_fit))
  ndraws <- brms::ndraws(flocker_fit)
  if(!is.null(draw_ids)){
    assertthat::assert_that(
      max(draw_ids) <= ndraws,
      msg = "some requested draw ids greater than total number of available draws"
    )
  }
  
  lps <- fitted_flocker(
    flocker_fit, draw_ids = draw_ids, new_data = NULL, allow_new_levels = FALSE, 
    response = TRUE, unit_level = FALSE
  )
  
  lik_type <- type_flocker_fit(flocker_fit)
  
  if (lik_type == "single") {
    psi_all <- lps$linpred_occ[ , 1, ] # first index is unit, second is visit, third is draw
    n_unit <- nrow(psi_all)
    theta_all <- lps$linpred_det
    

    gp <- get_positions(flocker_fit)
    obs <- matrix(flocker_fit$data$ff_y[gp], nrow = nrow(gp), ncol = ncol(gp))
      
    # get emission likelihoods
    el_0 <- el_1 <- matrix(nrow = nrow(psi_all), ncol = ncol(psi_all))
      
    for(i in seq_len(ncol(psi_all))){
      el_0[ , i] <- mapply(function(a, b){emission_likelihood(0, a, b)}, asplit(obs, 1), asplit(theta_all[ , , i], 1))
      el_1[ , i] <- mapply(function(a, b){emission_likelihood(1, a, b)}, asplit(obs, 1), asplit(theta_all[ , , i], 1))
    }
      
    assertthat::assert_that(identical(dim(psi_all), dim(el_0)))
      
    ll <- log((psi_all * el_1) + (1 - psi_all) * el_0) |>
      t()
  } else {
    stop("log_lik_flocker is currently implemented only for standard single-season models")
  }
  ll
}



#' A log-likelihood function for the rep-constant occupancy model, sufficient for
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


#' An R implementation of the rep constant lpmf without the binomial coefficient
#' @param y number of detections
#' @param mu logit-scale detection probability
#' @param occ logit-scale occupancy probability
#' @param trials number of reps
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
#' @param trials number of reps
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
