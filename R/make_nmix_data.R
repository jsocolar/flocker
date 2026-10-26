#' Format single-season N-mixture count data
#'
#' @param obs A numeric matrix or data frame of nonnegative integer counts.
#'   Rows are closure units and columns are repeated visits. NA denotes an
#'   unsampled visit; missing visits must trail within each row. Every unit
#'   must have at least one observed visit. At least two columns are required
#'   unless \code{newdata_checks = TRUE}; the last column cannot be entirely NA.
#' @param unit_covs An optional data frame with one row per unit, containing
#'   covariates constant across visits. Missing covariates are not allowed.
#' @param event_covs An optional named list of matrices or data frames with
#'   the same dimensions as \code{obs}, containing visit-varying covariates.
#'   Missing covariates are allowed only at unsampled visits.
#' @param quiet Hide formatting messages?
#' @param newdata_checks Allow a single visit column and mark the result as
#'   new data rather than fitting data. All other input checks still apply.
#' @return A \code{flocker_data} list containing \code{data}, \code{metadata},
#'   and \code{auxiliary}. Metadata identify the model as \code{n_mix}.
#' @details Observed visits remain separate even without event covariates.
#'   Unit covariates repeat at their actual values on every retained visit;
#'   centering and spline construction therefore use visit-expanded data.
#'   The auxiliary list holds count-only coefficients and indexing for exact
#'   marginalization with binomial detection, separately from predictor rows.
#'   Only single-season, one-level data are supported. This function prepares
#'   N-mixture data; \code{flock()} continues to accept occupancy data only.
#' @export
#' @examples
#' counts <- rbind(c(2, 3, 0), c(0, 1, NA))
#' make_nmix_data(counts, quiet = TRUE)
make_nmix_data <- function(obs, unit_covs = NULL, event_covs = NULL,
                           quiet = FALSE, newdata_checks = FALSE) {
  make_flocker_data_(obs = obs, unit_covs = unit_covs, event_covs = event_covs,
                    quiet = quiet, newdata_checks = newdata_checks,
                    type = "single", model = "n_mix")
}
