#' Get posterior latent abundance
#'
#' Reconstruct latent abundance for a fitted Poisson-binomial N-mixture model,
#' optionally conditioning on each unit's observed visit counts.
#' @inheritParams fitted_nmix
#' @param history_condition Condition on the observed visit history (TRUE), or
#'   use the fitted Poisson abundance distribution without directly conditioning
#'   on that history (FALSE). Model parameters always come from the fitted
#'   posterior, including when this argument is FALSE.
#' @param sample Return integer abundance draws (TRUE), or the conditional mean
#'   abundance for each parameter draw (FALSE).
#' @return A unit-by-draw matrix, retaining unit names and singleton dimensions.
#'   With \code{sample = FALSE}, variation across columns reflects parameter
#'   uncertainty only; intervals over those means are not intervals for latent
#'   abundance. Use \code{sample = TRUE} to include latent abundance uncertainty.
#' @details With history conditioning, abundance is sampled exactly as a finite
#'   mixture of shifted Poisson distributions, using the same falling-factorial
#'   expansion as the likelihood. No abundance truncation is required. Missing
#'   visits are omitted. For \code{new_data}, counts supplied to
#'   \code{make_nmix_data()} are treated as observed histories when
#'   \code{history_condition = TRUE}; with FALSE their values are ignored.
#' @export
get_N <- function(nmix_fit, draw_ids = NULL, history_condition = TRUE,
                  sample = FALSE, new_data = NULL, re_formula = NULL,
                  allow_new_levels = FALSE,
                  sample_new_levels = "uncertainty") {
  if (!is_one_logical(history_condition)) {
    stop("history_condition must be TRUE or FALSE.")
  }
  if (!is_one_logical(sample)) stop("sample must be TRUE or FALSE.")
  lps <- fitted_nmix(
    nmix_fit, components = if (history_condition) c("abund", "det") else "abund",
    new_data = new_data, response = FALSE, draw_ids = draw_ids,
    re_formula = re_formula, allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels
  )
  out <- lps$abund
  if (!history_condition) {
    out[] <- if (sample) stats::rpois(length(out), exp(out)) else exp(out)
    return(out)
  }
  object <- if (is.null(new_data)) nmix_fit else new_data
  metadata <- get_flocker_metadata(object)
  # brms retains the data supplied by stanvars; new data carry their own
  # coefficients. Reuse these rather than repeating the quadratic precomputation.
  auxiliary <- if (is.null(new_data)) {
    lapply(nmix_fit$stanvars, function(x) x$sdata)
  } else {
    new_data$auxiliary
  }
  for (u in seq_len(nrow(out))) {
    positions <- metadata$original_to_packed[u, ]
    visited <- !is.na(positions)
    y <- object$data$ff_y[positions[visited]]
    det <- matrix(lps$det[u, visited, , drop = FALSE], nrow = length(y))
    coef_rows <- seq.int(auxiliary$ff_coef_offsets[u],
                         auxiliary$ff_coef_offsets[u + 1L] - 1L)
    out[u, ] <- nmix_poisson_binomial_abundance(
      y, lps$abund[u, ], det, sample, auxiliary$ff_log_coef[coef_rows]
    )
  }
  out
}

# For M = max(y), q = prod(1-p), and x = lambda*q, the expansion gives
# N | y = M + I + Poisson(x), with Pr(I=i) proportional to c_i*x^i.
# The default coefficients also allow this helper to be tested independently.
nmix_poisson_binomial_abundance <- function(
    y, log_abundance, det_lp, sample,
    log_coef = nmix_binomial_log_coefficients(y)) {
  power <- seq_along(log_coef) - 1L
  # Stable log(1 - inv_logit(eta)), even for very large positive eta.
  log_q <- colSums(-pmax(det_lp, 0) - log1p(exp(-abs(det_lp))))
  log_x <- log_abundance + log_q
  vapply(log_x, function(lx) {
    log_weight <- log_coef
    # Avoid 0 * -Inf if the Poisson mean underflows to zero.
    positive <- power > 0L
    log_weight[positive] <- log_weight[positive] + power[positive] * lx
    weight <- exp(log_weight - max(log_weight))
    weight <- weight / sum(weight)
    if (sample) {
      max(y) + power[sample.int(length(power), 1L, prob = weight)] +
        stats::rpois(1L, exp(lx))
    } else {
      max(y) + sum(power * weight) + exp(lx)
    }
  }, numeric(1))
}
