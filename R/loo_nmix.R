#' Leave-one-unit-out cross-validation for N-mixture models
#' @param x An nmix fit or a nonempty list of nmix fits.
#' @param thin Positive integer thinning interval within each chain; NULL means 1.
#' @return A loo object, or a correspondingly named list of loo objects.
#' @details Each held-out observation is a unit's complete visit history.
#'   Abundance is marginalized, while model coefficients and group effects
#'   remain posterior draws. Inspect Pareto-k diagnostics before interpreting
#'   results. This does not implement leave-one-visit-out cross-validation,
#'   moment matching, or automatic refitting of problematic units.
#' @export
loo_nmix <- function(x, thin = NULL) {
  if (is.null(thin)) thin <- 1L
  if (!is.numeric(thin) || length(thin) != 1L || is.na(thin) ||
      !is.finite(thin) || thin < 1 || thin != floor(thin)) {
    stop("thin must be a positive integer or NULL.")
  }
  if (is_nmix_fit(x)) return(loo_nmix_onefit(x, thin))
  if (!is.list(x) || !length(x) || !all(vapply(x, is_nmix_fit, logical(1)))) {
    stop("x must be an nmix fit or a nonempty list of nmix fits.")
  }
  lapply(x, loo_nmix_onefit, thin = thin)
}

loo_nmix_onefit <- function(x, thin) {
  niter <- brms::niterations(x)
  nchains <- brms::nchains(x)
  keep <- seq(1L, niter, by = thin)
  ids <- unlist(lapply(seq_len(nchains), function(chain) (chain - 1L) * niter + keep))
  chains <- rep(seq_len(nchains), each = length(keep))
  ll <- log_lik_nmix(x, draw_ids = ids)
  # relative_eff needs likelihoods, not log likelihoods. Columnwise scaling
  # avoids underflow without changing their autocorrelation or relative ESS.
  likelihood <- exp(sweep(ll, 2L, apply(ll, 2L, max), "-"))
  out <- loo::loo(ll, r_eff = loo::relative_eff(likelihood, chain_id = chains))
  rownames(out$pointwise) <- colnames(ll)
  out
}

#' Compare N-mixture models using leave-one-unit-out cross-validation
#' @param model_list A list of at least two nmix fits to the same ordered counts.
#' @param model_names Optional distinct names for the models.
#' @inheritParams loo_nmix
#' @return The comparison table returned by \code{loo::loo_compare()}.
#' @export
loo_compare_nmix <- function(model_list, model_names = NULL, thin = NULL) {
  if (!is.list(model_list) || is_nmix_fit(model_list) || length(model_list) < 2L ||
      !all(vapply(model_list, is_nmix_fit, logical(1)))) {
    stop("model_list must contain at least two nmix fits.")
  }
  histories <- lapply(model_list, function(fit) {
    index <- get_flocker_metadata(fit)$original_to_packed
    array(fit$data$ff_y[index], dim(index), dimnames(index))
  })
  if (!all(vapply(histories, identical, logical(1), histories[[1L]]))) {
    stop("Models must use the same counts, missing visits, and unit/visit ordering and names.")
  }
  if (!is.null(model_names) && (!is.character(model_names) ||
      length(model_names) != length(model_list) || anyNA(model_names) ||
      any(!nzchar(model_names)) || anyDuplicated(model_names))) {
    stop("model_names must supply one distinct nonempty name per model.")
  }
  results <- loo_nmix(model_list, thin = thin)
  if (!is.null(model_names)) names(results) <- model_names
  loo::loo_compare(results)
}
