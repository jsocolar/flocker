#' Compute loo for flocker_fit objects
#' @param x a flocker_fit object or a list of flocker_fit objects
#' @param binom_coef Currently deprecated.
#' For visit-constant models, should the binomial coefficient
#' be included in the likelihood (\code{T}) or dropped (\code{F})? To compare
#' visit-constant models with alternative visit-varying specifications, the
#' binomial coefficient *must* be dropped.
#' @return a loo object or a list of loo objects
#' @export

loo_flock <- function(x, binom_coef = F) {
  if (!(is_flocker_fit(x) | is.list(x))) {
    stop("x must be a flocker_fit object or a list of flocker_fit objects")
  }
  if ("list" %in% class(x)) {
    for(i in 1:length(x)) {
      if (!is_flocker_fit(x[[i]])) {
        stop(paste0("x is a list, but x[[", i, "]] is not a flocker_fit object."))
      }
    }
  }
  
  if (is_flocker_fit(x)) {
    out <- loo_flock_onefit(x, binom_coef)
  } else {
    out <- list()
    for (i in 1:length(x)) {
      out[[i]] <- loo_flock_onefit(x[[i]], binom_coef)
    }
  }
  return(out)
}

#' Compute loo for a single flocker_fit object
#' @param x a flocker_fit object
#' @param binom_coef Currently deprecated.
#' For visit-constant models, should the binomial coefficient
#' be included in the likelihood (\code{T}) or dropped (\code{F})? To compare
#' visit-constant models with alternative visit-varying specifications, the
#' binomial coefficient *must* be dropped.
#' @return a loo object

loo_flock_onefit <- function(x, binom_coef) {
  type <- type_flocker_fit(x)
  if (type == "C") {
    # if (binom_coef) {
    #   .GlobalEnv$occupancy_C_lpmf <- occupancy_C_lpmf_with_coef
    # } else {
    #   .GlobalEnv$occupancy_C_lpmf <- occupancy_C_lpmf_without_coef
    # }
    out <- brms::loo(x)
  } else if (type == "V") {
    chain_id <- rep(c(1:dim(x$fit)[2]), each = dim(x$fit)[1])
    ll <- log_lik_V(x)
    out <- brms::loo(ll, r_eff = loo::relative_eff(ll, chain_id = chain_id))
  }
  return(out)
}
