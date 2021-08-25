#' Compute loo for flocker_fit objects
#' @param x a flocker_fit object or a list of flocker_fit objects
#' @return a loo object or a list of loo objects
#' @export

loo_flock <- function(x) {
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
    out <- loo_flock_onefit(x)
  } else {
    out <- list()
    for (i in 1:length(x)) {
      out[[i]] <- loo_flock_onefit(x[[i]])
    }
  }
  names(out) <- names(x)
  return(out)
}

#' Compute loo for a single flocker_fit object
#' @param x a flocker_fit object
#' @return a loo object

loo_flock_onefit <- function(x) {
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



#' LOO comparisons for flocker models. 
#' @param x a list of flocker_fit objects.
#' @param model_names An optional vector of names for the models.
#' @export

loo_compare_flock <- function(x, model_names = NULL) {
  if (!("list" %in% class(x))) {
    stop("x must be a list of flocker_fit objects.")
  }
  if (length(x) < 2L) {
    stop("x must contain at least two flocker_fit objects.")
  }
  occupancy_loo <- loo_flock(x)
  if (!is.null(model_names)) {
    names(occupancy_loo) <- model_names
  }
  loo::loo_compare(occupancy_loo)
}
