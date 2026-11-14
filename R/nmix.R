#' Fit a single-season N-mixture model
#'
#' Fit Poisson abundance with binomial detection, marginalizing abundance
#' exactly without a truncation bound.
#' @param f_abund A one-sided formula for log expected abundance. Covariates
#'   must be unit-level covariates supplied to \code{make_nmix_data()}.
#' @param f_det A one-sided formula for logit detection probability. Both unit
#'   and event covariates may be used. Combined \code{brmsformula} objects are
#'   not currently supported.
#' @param data A model-data object returned by \code{make_nmix_data()}.
#' @param data2 Additional data passed to brms, such as covariance matrices.
#' @param prior A \code{brmsprior} object. Detection uses the main brms
#'   parameter (omit \code{dpar}); abundance uses \code{dpar = "abund"}.
#' @param threads Must be NULL; threading is not currently supported.
#' @param ... Additional arguments passed to \code{brms::brm()}, or the
#'   corresponding brms function for code, data, or prior inspection.
#' @details Abundance covariates are repeated across visits for brms
#'   preprocessing, but the likelihood reads abundance only at each unit's
#'   first visit. Autoregressive terms are not supported.
#'   As in brms, custom transformations used by either formula must be visible
#'   from the detection formula's environment.
#'
#'   Detection intercepts and group/spline standard deviations default to
#'   \code{student_t(3, 0, 2.5)}, independently of observed counts, matching
#'   brms's ordinary binomial-scale defaults. Standard deviations are positive.
#'   Other defaults come from brms; in particular, the abundance intercept has
#'   a flat prior. Supply \code{prior} to override these defaults and inspect
#'   them with \code{get_nmix_prior()}.
#'
#'   The custom family uses identity links internally: its predictors already
#'   represent logit detection and log abundance. Native brms linear-predictor
#'   output therefore remains on those scales. Use \code{fitted_nmix()} for
#'   fitted parameters and \code{get_N()} for latent abundance. N-mixture
#'   prediction and pointwise log likelihoods are not yet implemented.
#' @return \code{nmix()} returns a \code{brmsfit} with additional class
#'   \code{nmix_fit} and flocker metadata.
#'   The inspection helpers return Stan code, Stan data, or a prior table.
#' @export
#' @examples
#' counts <- rbind(c(2, 3, 0), c(0, 1, NA))
#' nd <- make_nmix_data(counts, quiet = TRUE)
#' get_nmix_prior(~ 1, ~ 1, nd)
#' \dontrun{
#' nmix(~ 1, ~ 1, nd, chains = 4)
#' }
nmix <- function(f_abund, f_det, data, data2 = NULL, prior = NULL,
                   threads = NULL, ...) {
  nmix_("model", f_abund, f_det, data, data2, prior, threads, ...)
}

#' @rdname nmix
#' @export
nmix_stancode <- function(f_abund, f_det, data, data2 = NULL, prior = NULL,
                           threads = NULL, ...) {
  nmix_("code", f_abund, f_det, data, data2, prior, threads, ...)
}

#' @rdname nmix
#' @export
nmix_standata <- function(f_abund, f_det, data, data2 = NULL, prior = NULL,
                           threads = NULL, ...) {
  nmix_("data", f_abund, f_det, data, data2, prior, threads, ...)
}

#' @rdname nmix
#' @export
get_nmix_prior <- function(f_abund, f_det, data, data2 = NULL, prior = NULL,
                            threads = NULL, ...) {
  nmix_("prior", f_abund, f_det, data, data2, prior, threads, ...)
}

nmix_ <- function(output, f_abund, f_det, data, data2, prior, threads, ...) {
  dots <- list(...)
  validate_no_time_series_terms(list(f_abund, f_det, dots$autocor))
  if (!is.null(threads)) {
    stop("Threading is not currently supported for N-mixture models.")
  }
  for (name in c("f_abund", "f_det")) {
    f <- get(name)
    if (!inherits(f, "formula") || length(f) != 2L) {
      stop(name, " must be a one-sided formula beginning with ~.")
    }
  }
  if (!is_nmix_data(data)) {
    stop("nmix() requires N-mixture data from make_nmix_data().")
  }
  metadata <- get_flocker_metadata(data)
  validate_packed_schema(data)
  if (!identical(metadata$data_role, "model_data")) {
    stop("Fitting requires model_data; new_data cannot be used for fitting.")
  }
  if (any(!all.vars(f_abund) %in% c(metadata$unit_covs, names(data2)))) {
    stop("f_abund may use only unit covariates supplied to make_nmix_data().")
  }
  if (any(grepl("^ff_", all.vars(f_det)))) {
    stop("f_det cannot use reserved ff_ columns as predictors.")
  }
  if ("." %in% all.vars(f_det)) {
    stop("List detection covariates explicitly; ~ . would include auxiliary columns.")
  }
  if (any(c("formula", "family", "stanvars") %in% names(dots))) {
    stop("nmix() supplies formula, family, and stanvars internally.")
  }
  if (!is.null(prior) && !inherits(prior, "brmsprior")) {
    stop("prior must be a brmsprior object.")
  }
  auxiliary <- data$auxiliary
  required <- c("ff_nunit", "ff_ncoef", "ff_unit_offsets", "ff_coef_offsets",
                "ff_max_count", "ff_log_factorial", "ff_log_coef", "ff_coef_power")
  if (!is.list(auxiliary) || !all(required %in% names(auxiliary))) {
    stop("Missing N-mixture auxiliary data; recreate data with make_nmix_data().")
  }
  expected_offsets <- c(metadata$unit_start, nrow(data$data) + 1L)
  if (!identical(auxiliary$ff_unit_offsets, expected_offsets)) {
    stop("N-mixture unit offsets do not match the packed rows.")
  }

  # Construct calls rather than deparsing, preserving each formula's environment.
  det <- set_formula_response(f_det, quote(ff_y))
  abund <- set_formula_response(f_abund, quote(abund))
  formula <- brms::bf(det, abund)
  family <- nmix_poisson_binomial()
  stanvars <- nmix_poisson_binomial_stanvars(auxiliary)
  available <- brms::get_prior(formula, data = data$data, data2 = data2,
                               family = family, ...)
  # Provisional defaults pending a package-wide prior review. The identity
  # link is internal: count-derived locations/scales are wrong for detection.
  detection <- available[available$dpar == "" & available$nlpar == "" &
                           available$coef == "" & available$group == "" &
                           available$class %in% c("Intercept", "sd", "sds"), ]
  defaults <- NULL
  for (class in detection$class) {
    p <- brms::set_prior("student_t(3, 0, 2.5)", class = class)
    defaults <- if (is.null(defaults)) p else c(defaults, p)
  }
  # brms treats mu as an alias for the main parameter when specifying priors.
  if (!is.null(prior)) prior$dpar[prior$dpar == "mu"] <- ""
  prior <- if (is.null(defaults)) prior else if (is.null(prior)) defaults else
    c(defaults, prior, replace = TRUE)
  args <- c(list(f_use = formula, data = data$data, data2 = data2,
                 family = family, stanvars = stanvars, prior = prior), dots)

  # Formula transformations can create missing values even in validated inputs.
  # brms must not silently remove visits and invalidate the auxiliary offsets.
  standata <- do.call(flocker_fit_code_util, c(list(output = "data"), args))
  if (standata$N != nrow(data$data) ||
      !identical(as.integer(standata$Y), data$data$ff_y)) {
    stop("Formula processing changed the packed visits. Check for missing or invalid predictor values.")
  }
  if (output == "data") return(standata)
  if (output == "prior") {
    return(brms::validate_prior(prior, formula, data = data$data, data2 = data2,
                                family = family, ...))
  }
  out <- do.call(flocker_fit_code_util, c(list(output = output), args))
  if (output == "model") {
    class(out) <- c(class(out), "nmix_fit")
    attr(out, "flocker_metadata") <- metadata
    attr(out, "data_type") <- metadata$type
    attr(out, "flocker_version") <- flocker_version()
  }
  out
}
