#' Compute loo for flocker_fit objects
#' @param x a flocker_fit object or a list of flocker_fit objects
#' @param thin specify the amount of thinning required. 1 or NULL implies no thinning, 2 implies every other value, 3 every third, etc.
#' @return a loo object or a list of loo objects
#' @export
#' @examples 
#' \dontrun{
#' loo_flocker(example_flocker_model_single)
#' }
loo_flocker <- function(x, thin = NULL) {
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
    out <- loo_flocker_onefit(x, thin = thin)
  } else {
    out <- list()
    for (i in 1:length(x)) {
      out[[i]] <- loo_flocker_onefit(x[[i]], thin = thin)
    }
    names(out) <- names(x)
  }
  out
}

#' Compute loo for a single flocker_fit object
#' @param x a flocker_fit object
#' @param thin specify the amount of thinning required. 1 or NULL implies no thinning, 2 implies every other value, 3 every third, etc.
#' @return a loo object
#' @noRd
loo_flocker_onefit <- function(x, thin = NULL) {
  type <- type_flocker_fit(x)
  # do thinning 
  # thin == 1 | thin == NULL: retain all
  # thin == 2: retain every other 
  # thin == 3 retain every third, etc. 
  niter <- brms::niterations(x)
  nchains <- brms::nchains(x)
  
  if(is.null(thin)) {
    draw_ids <- 1:(niter*nchains) 
    chain_ids <- rep(c(1:nchains), each = niter)
  } else {
    iter_keep <- seq(1, niter, thin)
    draw_ids <- rep(seq(0, (niter-1)*nchains, niter), each=length(iter_keep)) + iter_keep
    chain_ids <- rep(1:nchains, each = length(iter_keep))
  }
  
  ll <- log_lik_flocker(x, draw_ids = draw_ids)
  loo::loo(ll, r_eff = loo::relative_eff(ll, chain_id = chain_ids))
}

#' LOO comparisons for flocker models. 
#' @param model_list a list of flocker_fit objects.
#' @param thin specify the amount of thinning required. 1 or NULL results in no 
#'    thinning, 2 retains every other value, 3 every third, etc.
#' @param model_names An optional vector of names for the models.
#' @return a `compare.loo` matrix
#' @export
#' @examples
#' \donttest{
#' ml <- rep(list(example_flocker_model_single), 3)
#' loo_compare_flocker(ml)
#' }
loo_compare_flocker <- function(model_list, model_names = NULL, thin = NULL) {
  if (!("list" %in% class(model_list))) {
    stop("model_list must be a list of flocker_fit objects.")
  }
  if (length(model_list) < 2L) {
    stop("model_list must contain at least two flocker_fit objects.")
  }
  occupancy_loo <- loo_flocker(model_list, thin = thin)
  if (!is.null(model_names)) {
    names(occupancy_loo) <- model_names
  }
  loo::loo_compare(occupancy_loo)
}
