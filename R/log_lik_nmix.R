#' Marginal log likelihoods for N-mixture models
#' @inheritParams fitted_nmix
#' @return A draw-by-unit matrix, retaining draw and unit names and singleton
#'   dimensions. Each column is the joint likelihood of all observed visits
#'   within that unit, with abundance marginalized exactly.
#' @details These are marginal likelihoods, not probabilities conditional on
#'   reconstructed abundance. Missing visits are omitted. For new data, the
#'   supplied counts are evaluated under the fitted parameter posterior.
#' @export
log_lik_nmix <- function(nmix_fit, new_data = NULL, draw_ids = NULL,
                         re_formula = NULL, allow_new_levels = FALSE,
                         sample_new_levels = "uncertainty") {
  lps <- fitted_nmix(nmix_fit, response = FALSE, new_data = new_data,
    draw_ids = draw_ids, re_formula = re_formula,
    allow_new_levels = allow_new_levels, sample_new_levels = sample_new_levels)
  object <- if (is.null(new_data)) nmix_fit else new_data
  metadata <- get_flocker_metadata(object)
  auxiliary <- if (is.null(new_data)) {
    lapply(nmix_fit$stanvars, function(x) x$sdata)
  } else new_data$auxiliary
  out <- t(lps$abund)
  for (u in seq_len(ncol(out))) {
    positions <- metadata$original_to_packed[u, ]
    visited <- !is.na(positions)
    y <- object$data$ff_y[positions[visited]]
    det <- matrix(lps$det[u, visited, , drop = FALSE], nrow = length(y))
    rows <- seq.int(auxiliary$ff_coef_offsets[u], auxiliary$ff_coef_offsets[u + 1L] - 1L)
    out[, u] <- nmix_poisson_binomial_log_lik(
      y, lps$abund[u, ], det, auxiliary$ff_log_coef[rows])
  }
  out
}

# Same stable falling-factorial evaluation as the Stan unit likelihood,
# vectorized over posterior draws; no abundance cutoff is required.
nmix_poisson_binomial_log_lik <- function(y, log_abundance, det_lp,
    log_coef = nmix_binomial_log_coefficients(y)) {
  log_q <- log_seen <- rep(0, length(log_abundance))
  log_seen[] <- -Inf
  M <- max(y)
  lp <- M * log_abundance - sum(lgamma(y + 1))
  for (j in seq_along(y)) {
    log_p <- log_inv_logit(det_lp[j, ])
    log_miss <- log1m_inv_logit(det_lp[j, ])
    log_seen <- vapply(seq_along(log_q), function(i)
      matrixStats::logSumExp(c(log_seen[i], log_q[i] + log_p[i])), numeric(1))
    log_q <- log_q + log_miss
    if (y[j] > 0) lp <- lp + y[j] * log_p
    if (M > y[j]) lp <- lp + (M - y[j]) * log_miss
  }
  power <- seq_along(log_coef) - 1L
  polynomial <- vapply(log_abundance + log_q, function(x) {
    terms <- log_coef
    positive <- power > 0L
    terms[positive] <- terms[positive] + power[positive] * x
    matrixStats::logSumExp(terms)
  }, numeric(1))
  lp - exp(log_abundance + log_seen) + polynomial
}
