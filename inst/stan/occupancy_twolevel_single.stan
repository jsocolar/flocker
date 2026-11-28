// Each group is a self-contained one-level occupancy dataset, conditional on
// the group being occupied. Its outer occupancy indicator is marginalized here.
// Q and group_Q repeat the unit and group detection indicators at every row.
// Predictor inputs are logits: mu for detection conditional on unit occupancy,
// occ for unit occupancy conditional on group occupancy, and Omega for group
// occupancy. These inputs are not probabilities or log-probabilities.
real occupancy_twolevel_single_lpmf(array[] int y, vector mu, vector occ,
                                    vector Omega, array[] int group_size,
                                    array[] int unit_size, array[] int Q,
                                    array[] int group_Q) {
  int n = num_elements(y);
  int group_start = 1;
  real lp = 0;
  if (n < 1 || num_elements(mu) != n || num_elements(occ) != n
      || num_elements(Omega) != n || num_elements(group_size) != n
      || num_elements(unit_size) != n || num_elements(Q) != n
      || num_elements(group_Q) != n) {
    reject("Two-level inputs must have the same nonzero packed-row length.");
  }
  while (group_start <= n) {
    int group_end = group_start + group_size[group_start] - 1;
    real lp_present;
    if (group_start < 1 || group_end < group_start || group_end > n) {
      reject("Invalid group block.");
    }
    lp_present = bernoulli_logit_lpmf(1 | Omega[group_start])
      + occupancy_single_lpmf(y[group_start:group_end] |
          mu[group_start:group_end], occ[group_start:group_end],
          unit_size[group_start:group_end],
          Q[group_start:group_end]);
    if (group_Q[group_start] == 1) {
      lp += lp_present;
    } else {
      lp += log_sum_exp(bernoulli_logit_lpmf(0 | Omega[group_start]), lp_present);
    }
    group_start = group_end + 1;
  }
  return lp;
}
