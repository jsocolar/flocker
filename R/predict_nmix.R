#' Get posterior count predictions from an N-mixture model
#'
#' Generate replicated visit counts from a fitted Poisson-binomial N-mixture
#' model, sharing one sampled abundance across all visits within each unit.
#' @inheritParams fitted_nmix
#' @inheritParams get_N
#' @return A unit-by-visit-by-draw array of counts, retaining original names,
#'   singleton dimensions, and NA at unsampled visits.
#' @details By default, abundance is drawn from its fitted Poisson distribution.
#'   With \code{history_condition = TRUE}, abundance is instead drawn from its
#'   exact conditional distribution given the observed visit counts. In either
#'   case, new binomial counts are drawn: conditioning does not reproduce the
#'   observed counts, and replicated counts can be smaller than observed counts.
#'
#'   For \code{new_data}, supplied counts are treated as observed histories only
#'   when conditioning. To predict at a visit without an observed count, supply
#'   a nonmissing placeholder count and use \code{history_condition = FALSE};
#'   NA denotes an unsampled visit and remains NA in the output. Model parameters
#'   always come from the fitted posterior. Abundance and detection use the same
#'   jointly prepared effects, including for new group levels.
#' @export
predict_nmix <- function(nmix_fit, draw_ids = NULL, history_condition = FALSE,
                          new_data = NULL, re_formula = NULL,
                          allow_new_levels = FALSE,
                          sample_new_levels = "uncertainty") {
  if (!is_one_logical(history_condition)) {
    stop("history_condition must be TRUE or FALSE.")
  }
  lps <- fitted_nmix(
    nmix_fit, response = FALSE, draw_ids = draw_ids, new_data = new_data,
    re_formula = re_formula, allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels
  )
  abundance <- nmix_abundance_from_predictors(
    nmix_fit, lps, history_condition, TRUE, new_data
  )
  out <- lps$det
  # Matrix rows are units and columns are draws; repeat each draw's unit
  # abundances for every visit in the corresponding array slice.
  sizes <- abundance[, rep(seq_len(ncol(abundance)), each = dim(out)[2]), drop = FALSE]
  visited <- !is.na(out)
  out[visited] <- stats::rbinom(sum(visited), sizes[visited], stats::plogis(out[visited]))
  out
}
