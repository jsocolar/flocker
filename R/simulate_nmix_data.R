#' Simulate single-season Poisson-binomial N-mixture data
#'
#' Simulate independent closure units with Poisson abundance and binomial visit
#' counts. One unit covariate affects abundance and detection; one visit
#' covariate affects detection. Group effects are not included.
#' @param n_pt Positive integer number of closure units.
#' @param n_rep Positive integer maximum number of visits per unit.
#' @param params Named numeric vector with coefficients \code{abund_intercept},
#'   \code{abund_slope_unit}, \code{det_intercept}, \code{det_slope_unit}, and
#'   \code{det_slope_visit}. Abundance uses a log link and detection a logit link.
#'   Defaults are fixed coefficients, not draws from the fitting priors.
#' @param covariates Optional list containing \code{uc1}, a numeric vector of
#'   length \code{n_pt}, and/or \code{ec1}, a numeric matrix with \code{n_pt} rows
#'   and \code{n_rep} columns. Omitted covariates are independent standard normals.
#' @param ragged_rep If TRUE, sample each unit's number of visits uniformly from
#'   1 through \code{n_rep}; trailing unvisited counts and event covariates are NA.
#' @param seed Random seed; NULL uses and updates the current RNG state. A supplied
#'   seed leaves the caller's RNG state unchanged.
#' @return A list with \code{obs}, \code{unit_covs}, and \code{event_covs} ready
#'   for \code{make_nmix_data()}, plus realized \code{abundance}, expected
#'   abundance \code{lambda}, visit detection probabilities \code{detection},
#'   coefficients \code{params}, and complete generated \code{covariates}.
#' @export
#' @examples
#' sim <- simulate_nmix_data(n_pt = 20, ragged_rep = TRUE)
#' data <- make_nmix_data(sim$obs, sim$unit_covs, sim$event_covs, quiet = TRUE)
simulate_nmix_data <- function(n_pt = 200, n_rep = 3,
    params = c(abund_intercept = log(5), abund_slope_unit = .3,
               det_intercept = 0, det_slope_unit = .3, det_slope_visit = .3),
    covariates = NULL, ragged_rep = FALSE, seed = 1) {
  for (name in c("n_pt", "n_rep")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 1 || value != floor(value) || value > .Machine$integer.max) {
      stop(name, " must be a positive integer.")
    }
  }
  required <- c("abund_intercept", "abund_slope_unit", "det_intercept",
                 "det_slope_unit", "det_slope_visit")
  if (!is.numeric(params) || length(params) != length(required) ||
      anyDuplicated(names(params)) || !setequal(names(params), required) ||
      any(!is.finite(params))) stop("params must supply all five named finite coefficients.")
  if (!is_one_logical(ragged_rep)) stop("ragged_rep must be TRUE or FALSE.")
  if (!is.null(covariates) && (!is.list(covariates) ||
      (length(covariates) && (is.null(names(covariates)) ||
       anyNA(names(covariates)) || anyDuplicated(names(covariates)) ||
       any(!names(covariates) %in% c("uc1", "ec1")))))) {
    stop("covariates must be a named list containing uc1 and/or ec1.")
  }
  if (!is.null(covariates$uc1) && (!is.numeric(covariates$uc1) ||
      !is.null(dim(covariates$uc1)) || length(covariates$uc1) != n_pt ||
      any(!is.finite(covariates$uc1)))) stop("uc1 must be a finite numeric vector of length n_pt.")
  if (!is.null(covariates$ec1) && (!is.numeric(covariates$ec1) ||
      !is.matrix(covariates$ec1) || !identical(dim(covariates$ec1), as.integer(c(n_pt, n_rep))) ||
      any(!is.finite(covariates$ec1)))) stop("ec1 must be a finite numeric n_pt by n_rep matrix.")
  simulate <- function() {
    uc1 <- covariates$uc1
    if (is.null(uc1)) uc1 <- stats::rnorm(n_pt)
    ec1 <- covariates$ec1
    if (is.null(ec1)) ec1 <- matrix(stats::rnorm(n_pt * n_rep), n_pt, n_rep)
    units <- paste0("unit_", seq_len(n_pt))
    visits <- paste0("visit_", seq_len(n_rep))
    dimnames(ec1) <- list(units, visits)
    lambda <- exp(params[["abund_intercept"]] + params[["abund_slope_unit"]] * uc1)
    if (any(!is.finite(lambda))) stop("Abundance means must be finite; reduce coefficients or covariate magnitudes.")
    abundance <- stats::rpois(n_pt, lambda)
    detection <- stats::plogis(params[["det_intercept"]] +
      params[["det_slope_unit"]] * uc1 + params[["det_slope_visit"]] * ec1)
    obs <- matrix(stats::rbinom(n_pt * n_rep, rep(abundance, n_rep), detection),
                   n_pt, n_rep, dimnames = list(units, visits))
    event <- ec1
    if (ragged_rep) {
      n_visit <- sample.int(n_rep, n_pt, replace = TRUE)
      missing <- col(obs) > n_visit
      obs[missing] <- event[missing] <- detection[missing] <- NA_real_
    }
    names(abundance) <- names(lambda) <- units
    list(obs = obs, unit_covs = data.frame(uc1 = uc1, row.names = units),
         event_covs = list(ec1 = event), abundance = abundance, lambda = lambda,
         detection = detection, params = params, covariates = list(uc1 = uc1, ec1 = ec1))
  }
  if (is.null(seed)) simulate() else withr::with_seed(seed, simulate())
}
