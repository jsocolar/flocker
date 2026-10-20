// Nested fixed-width packing: series starts index units; unit starts index rows.
// A unit with no visits has one row with y = -99. Its emission log probability
// is zero under either state, but its occupancy transition still contributes.
real occupancy_multi_loglik(array[] int y, vector mu, vector occ,
                            vector colo, vector ex,
                            array[] int n_series, array[] int series_starts,
                            array[] int n_unit, array[] int unit_starts,
                            array[] int Q) {
  int n = num_elements(y);
  int series;
  int units;
  real lp = 0;
  if (n < 1 || num_elements(mu) != n || num_elements(occ) != n
      || num_elements(colo) != n || num_elements(ex) != n
      || num_elements(n_series) != n || num_elements(series_starts) != n
      || num_elements(n_unit) != n || num_elements(unit_starts) != n
      || num_elements(Q) != n) {
    reject("Multiseason inputs must have the same nonzero packed-row length.");
  }
  series = n_series[1];
  units = n_unit[1];
  if (series < 1 || units < series || units > n
      || series_starts[1] != 1 || unit_starts[1] != 1) {
    reject("Invalid multiseason block counts or first starts.");
  }
  for (s in 1:series) {
    int series_start = series_starts[s];
    int series_end = s < series ? series_starts[s + 1] - 1 : units;
    real alpha_0 = 0;
    real alpha_1 = 0;
    if (series_start < 1 || series_end < series_start || series_end > units) {
      reject("Invalid series block.");
    }
    for (u in series_start:series_end) {
      int unit_start = unit_starts[u];
      int unit_end = u < units ? unit_starts[u + 1] - 1 : n;
      real emission = 0;
      if (unit_start < 1 || unit_end < unit_start || unit_end > n) {
        reject("Invalid unit block.");
      }
      if (y[unit_start] == -99) {
        if (unit_end != unit_start || Q[u] != 0) {
          reject("An unvisited unit must have one placeholder row and Q = 0.");
        }
      } else {
        emission = bernoulli_logit_lpmf(y[unit_start:unit_end] |
                                       mu[unit_start:unit_end]);
      }
      if (u == series_start) {
        alpha_0 = bernoulli_logit_lpmf(0 | occ[unit_start]);
        alpha_1 = bernoulli_logit_lpmf(1 | occ[unit_start]) + emission;
      } else {
        real previous_0 = alpha_0;
        real previous_1 = alpha_1;
        alpha_0 = log_sum_exp(previous_0 + bernoulli_logit_lpmf(0 | colo[unit_start]),
                              previous_1 + bernoulli_logit_lpmf(1 | ex[unit_start]));
        alpha_1 = log_sum_exp(previous_0 + bernoulli_logit_lpmf(1 | colo[unit_start]),
                              previous_1 + bernoulli_logit_lpmf(0 | ex[unit_start]))
                  + emission;
      }
      if (Q[u] == 1) alpha_0 = negative_infinity();
    }
    lp += log_sum_exp(alpha_0, alpha_1);
  }
  return lp;
}

real occupancy_multi_colex_lpmf(array[] int y, vector mu, vector occ,
                                vector colo, vector ex,
                                array[] int n_series, array[] int series_starts,
                                array[] int n_unit, array[] int unit_starts,
                                array[] int Q) {
  return occupancy_multi_loglik(y, mu, occ, colo, ex,
                                n_series, series_starts, n_unit, unit_starts, Q);
}

real occupancy_multi_colex_eq_lpmf(array[] int y, vector mu,
                                   vector colo, vector ex,
                                   array[] int n_series, array[] int series_starts,
                                   array[] int n_unit, array[] int unit_starts,
                                   array[] int Q) {
  // Equilibrium odds are colonization / extinction; compute on the log scale.
  vector[num_elements(colo)] occ = log_inv_logit(colo) - log_inv_logit(ex);
  return occupancy_multi_loglik(y, mu, occ, colo, ex,
                                n_series, series_starts, n_unit, unit_starts, Q);
}

real occupancy_multi_autologistic_lpmf(array[] int y, vector mu, vector occ,
                                       vector colo, vector autologistic,
                                       array[] int n_series, array[] int series_starts,
                                       array[] int n_unit, array[] int unit_starts,
                                       array[] int Q) {
  vector[num_elements(colo)] ex = -(colo + autologistic);
  return occupancy_multi_loglik(y, mu, occ, colo, ex,
                                n_series, series_starts, n_unit, unit_starts, Q);
}

real occupancy_multi_autologistic_eq_lpmf(array[] int y, vector mu,
                                          vector colo, vector autologistic,
                                          array[] int n_series, array[] int series_starts,
                                          array[] int n_unit, array[] int unit_starts,
                                          array[] int Q) {
  vector[num_elements(colo)] ex = -(colo + autologistic);
  vector[num_elements(colo)] occ = log_inv_logit(colo) - log_inv_logit(ex);
  return occupancy_multi_loglik(y, mu, occ, colo, ex,
                                n_series, series_starts, n_unit, unit_starts, Q);
}
