#' Fitted abundance and detection parameters
#'
#' Extract expected abundance and detection probabilities from a fitted
#' Poisson-binomial N-mixture model. These are model parameters, not draws of
#' latent abundance or predicted counts, and do not condition latent abundance
#' on the observed visit history.
#' @param nmix_fit A fitted model returned by \code{nmix()}.
#' @param components Components to return: \code{"abund"} and/or \code{"det"}.
#' @param new_data Optional formatted data from \code{make_nmix_data()}.
#'   Both model data and data made with \code{newdata_checks = TRUE} are accepted.
#'   If NULL, use the fitted data.
#' @param response Return expected abundance and detection probabilities
#'   (TRUE), or log abundance and logit detection (FALSE).
#' @param summarise Return posterior means and credible interval bounds instead
#'   of individual draws.
#' @param CI Two probabilities defining the credible interval.
#' @param draw_ids Optional indices of posterior draws to use.
#' @param re_formula Group-level effects to include; NULL includes all effects,
#'   and NA excludes them.
#' @param allow_new_levels Allow new group levels in \code{new_data}.
#' @param sample_new_levels How to sample effects for new group levels; passed
#'   to \code{brms::prepare_predictions()}.
#' @return A named list. \code{abund} is a unit-by-draw matrix; \code{det} is a
#'   unit-by-visit-by-draw array, with NA at unsampled visits. Original unit and
#'   visit names are preserved. With \code{summarise = TRUE}, the final dimension
#'   instead contains \code{mean} and the requested quantiles. Dimensions are
#'   retained for single units, visits, and draws.
#' @export
fitted_nmix <- function(nmix_fit, components = c("abund", "det"),
                         new_data = NULL, response = TRUE, summarise = FALSE,
                         CI = c(0.05, 0.95), draw_ids = NULL,
                         re_formula = NULL, allow_new_levels = FALSE,
                         sample_new_levels = "uncertainty") {
  if (!is_nmix_fit(nmix_fit)) {
    stop("nmix_fit must be a fitted model returned by nmix().")
  }
  components <- match.arg(components, c("abund", "det"), several.ok = TRUE)
  draw_ids <- validate_fitted_params(
    nmix_fit, new_data, response, summarise, CI, draw_ids
  )
  object <- if (is.null(new_data)) nmix_fit else new_data
  metadata <- get_flocker_metadata(object)
  prediction_data <- object$data
  if ("Intercept" %in% names(nmix_fit$data) &&
      !"Intercept" %in% names(prediction_data)) prediction_data$Intercept <- 1

  # Prepare once so correlated effects for new levels are sampled jointly
  # across both distributional parameters.
  prep <- brms::prepare_predictions(
    nmix_fit, newdata = prediction_data, draw_ids = draw_ids,
    re_formula = re_formula, allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels, check_response = FALSE
  )
  out <- lapply(components, function(component) {
    dpar <- if (component == "abund") "abund" else "mu"
    lp <- t(brms::posterior_epred(prep, dpar = dpar, nlpar = NULL,
                                  sort = FALSE, scale = "linear", summary = FALSE))
    format_nmix_fitted(lp, component, metadata, response, summarise, CI, draw_ids)
  })
  stats::setNames(out, components)
}

# Reconstruct each component at its intrinsic ecological level.
format_nmix_fitted <- function(lp, component, metadata, response, summarise, CI,
                                draw_ids) {
  if (response) lp <- if (component == "abund") exp(lp) else stats::plogis(lp)
  values <- if (summarise) as.matrix(summarise_fun(lp, CI)) else lp
  unit_names <- metadata$original_dimnames[[1L]]
  if (is.null(unit_names)) unit_names <- paste0("site_", seq_len(metadata$original_dim[1L]))
  dn <- list(unit_names)
  if (component == "abund") {
    positions <- metadata$unit_start
  } else {
    positions <- metadata$original_to_packed
    visit_names <- metadata$original_dimnames[[2L]]
    if (is.null(visit_names)) visit_names <- paste0("visit_", seq_len(metadata$original_dim[2L]))
    dn <- c(dn, list(visit_names))
  }
  labels <- if (summarise) c("mean", paste0("Q", 100 * CI)) else paste0("draw_", draw_ids)
  out <- reshape_fun(values, positions)
  dimnames(out) <- c(dn, list(labels))
  out
}
