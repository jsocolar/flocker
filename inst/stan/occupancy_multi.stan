// Series and units store their row counts at their own first rows.
// Complete series can be sliced without rebasing any auxiliary columns.
// A unit with no visits has one row with y = -99. Its emission log probability
// is zero under either state, but its occupancy transition still contributes.
// Predictor inputs are logits: mu for detection, occ for initial occupancy,
// colo for colonization, and ex for extinction. These inputs are not
// probabilities or log-probabilities.
real occupancy_multi_loglik(array[] int y, vector mu, vector occ,
                            vector colo, vector ex,
                            array[] int series_size, array[] int unit_size,
                            array[] int Q) {
  int n = num_elements(y);
  int series_start = 1;
  real lp = 0;
  if (n < 1 || num_elements(mu) != n || num_elements(occ) != n
      || num_elements(colo) != n || num_elements(ex) != n
      || num_elements(series_size) != n || num_elements(unit_size) != n
      || num_elements(Q) != n) {
    reject("Multiseason inputs must have the same nonzero packed-row length.");
  }
  while (series_start <= n) {
    int series_end = series_start + series_size[series_start] - 1;
    int unit_start = series_start;
    real alpha_0 = 0;
    real alpha_1 = 0;
    if (series_start < 1 || series_end < series_start || series_end > n) {
      reject("Invalid series block.");
    }
    while (unit_start <= series_end) {
      int unit_end = unit_start + unit_size[unit_start] - 1;
      real emission = 0;
      if (unit_start < series_start || unit_end < unit_start || unit_end > series_end) {
        reject("Invalid unit block.");
      }
      if (y[unit_start] == -99) {
        if (unit_end != unit_start || Q[unit_start] != 0) {
          reject("An unvisited unit must have one placeholder row and Q = 0.");
        }
      } else {
        emission = bernoulli_logit_lpmf(y[unit_start:unit_end] |
                                       mu[unit_start:unit_end]);
      }
      if (unit_start == series_start) {
        alpha_0 = Q[unit_start] == 1 ? negative_infinity()
                                     : bernoulli_logit_lpmf(0 | occ[unit_start]);
        alpha_1 = bernoulli_logit_lpmf(1 | occ[unit_start]) + emission;
      } else {
        real previous_0 = alpha_0;
        real previous_1 = alpha_1;
        if (Q[unit_start] == 1) {
          alpha_0 = negative_infinity();
        } else {
          alpha_0 = log_sum_exp(previous_0 + bernoulli_logit_lpmf(0 | colo[unit_start]),
                                previous_1 + bernoulli_logit_lpmf(1 | ex[unit_start]));
        }
        alpha_1 = log_sum_exp(previous_0 + bernoulli_logit_lpmf(1 | colo[unit_start]),
                              previous_1 + bernoulli_logit_lpmf(0 | ex[unit_start]))
                  + emission;
      }
      unit_start = unit_end + 1;
    }
    lp += log_sum_exp(alpha_0, alpha_1);
    series_start = series_end + 1;
  }
  return lp;
}

// Predictor inputs are logits: mu for detection, occ for initial occupancy,
// colo for colonization, and ex for extinction. These inputs are not
// probabilities or log-probabilities.
real occupancy_multi_colex_lpmf(array[] int y, vector mu, vector occ,
                                vector colo, vector ex,
                                array[] int series_size, array[] int unit_size,
                                array[] int Q) {
  return occupancy_multi_loglik(y, mu, occ, colo, ex,
                                series_size, unit_size, Q);
}

// Predictor inputs are logits: mu for detection, colo for colonization, and
// ex for extinction. Initial occupancy is derived internally at equilibrium.
// These inputs are not probabilities or log-probabilities.
real occupancy_multi_colex_eq_lpmf(array[] int y, vector mu,
                                   vector colo, vector ex,
                                   array[] int series_size, array[] int unit_size,
                                   array[] int Q) {
  // Equilibrium odds are colonization / extinction; compute on the log scale.
  vector[num_elements(colo)] occ = log_inv_logit(colo) - log_inv_logit(ex);
  return occupancy_multi_loglik(y, mu, occ, colo, ex,
                                series_size, unit_size, Q);
}

// mu, occ, and colo are logits of detection, initial occupancy, and
// colonization probabilities. autologistic is an additive log-odds effect of
// previous occupancy: colo + autologistic is the logit persistence probability.
// No predictor input is a probability or a log-probability.
real occupancy_multi_autologistic_lpmf(array[] int y, vector mu, vector occ,
                                       vector colo, vector autologistic,
                                       array[] int series_size, array[] int unit_size,
                                       array[] int Q) {
  vector[num_elements(colo)] ex = -(colo + autologistic);
  return occupancy_multi_loglik(y, mu, occ, colo, ex,
                                series_size, unit_size, Q);
}

// mu and colo are logits of detection and colonization probabilities.
// autologistic is an additive log-odds effect of previous occupancy:
// colo + autologistic is the logit persistence probability. Initial occupancy
// is derived internally at equilibrium. No predictor input is a probability
// or a log-probability.
real occupancy_multi_autologistic_eq_lpmf(array[] int y, vector mu,
                                          vector colo, vector autologistic,
                                          array[] int series_size, array[] int unit_size,
                                          array[] int Q) {
  vector[num_elements(colo)] ex = -(colo + autologistic);
  vector[num_elements(colo)] occ = log_inv_logit(colo) - log_inv_logit(ex);
  return occupancy_multi_loglik(y, mu, occ, colo, ex,
                                series_size, unit_size, Q);
}
