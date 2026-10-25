// Exact closed-population Poisson abundance / binomial detection likelihood.
// Haines' falling-factorial expansion, with visit-specific detection.
// Predictor inputs: mu is LOGIT detection probability; log_abundance is LOG
// expected abundance. Coefficients and factorial terms are computed once in R.
real n_mixture_poisson_unit_lpmf(array[] int y, vector mu,
                                real log_abundance, int max_count,
                                real log_factorial, vector log_coef,
                                vector coef_power) {
  real log_q = 0;
  real log_seen = negative_infinity();
  real lp = max_count * log_abundance - log_factorial;
  for (t in 1:num_elements(y)) {
    real log_p = log_inv_logit(mu[t]);
    real log_miss = log1m_inv_logit(mu[t]);
    // Disjoint first-detection probabilities compute log(1-q) without
    // subtracting nearly equal numbers, including very small detection.
    log_seen = log_sum_exp(log_seen, log_q + log_p);
    log_q += log_miss;
    // Algebraically sum(y * logit(p)) + max_count * log(q), but without
    // cancellation when detection is near one.
    lp += y[t] * log_p + (max_count - y[t]) * log_miss;
  }
  // Powers 0,...,m are data, precomputed alongside each unit's coefficients.
  return lp - exp(log_abundance + log_seen)
         + log_sum_exp(log_coef + coef_power * (log_abundance + log_q));
}

// Units occupy contiguous visit blocks. Abundance covariates repeat on each
// visit row; only abund at unit_start is read. Both offset arrays contain a
// terminal sentinel. Coefficient storage is separate from predictor rows.
// Predictor inputs: mu is logit detection; abund is log expected abundance.
real n_mixture_poisson_lpmf(array[] int y, vector mu, vector abund,
                           array[] int unit_offsets, array[] int coef_offsets,
                           array[] int max_count, vector log_factorial,
                           vector log_coef, vector coef_power) {
  real lp = 0;
  for (u in 1:num_elements(max_count)) {
    int unit_start = unit_offsets[u];
    int unit_end = unit_offsets[u + 1] - 1;
    lp += n_mixture_poisson_unit_lpmf(y[unit_start:unit_end] |
      mu[unit_start:unit_end], abund[unit_start], max_count[u], log_factorial[u],
      log_coef[coef_offsets[u]:(coef_offsets[u + 1] - 1)],
      coef_power[coef_offsets[u]:(coef_offsets[u + 1] - 1)]);
  }
  return lp;
}
