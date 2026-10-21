// Each group is a self-contained one-level occupancy dataset, conditional on
// the group being occupied. Its outer occupancy indicator is marginalized here.
// Q and group_Q repeat the unit and group detection indicators at every row.
real occupancy_twolevel_single_lpmf(array[] int y, vector mu, vector occ,
                                    vector Omega, array[] int n_group,
                                    array[] int group_starts, array[] int n_unit,
                                    array[] int unit_starts, array[] int Q,
                                    array[] int group_Q) {
  int n = num_elements(y);
  int groups;
  real lp = 0;
  if (n < 1 || num_elements(mu) != n || num_elements(occ) != n
      || num_elements(Omega) != n || num_elements(n_group) != n
      || num_elements(group_starts) != n || num_elements(n_unit) != n
      || num_elements(unit_starts) != n || num_elements(Q) != n
      || num_elements(group_Q) != n) {
    reject("Two-level inputs must have the same nonzero packed-row length.");
  }
  groups = n_group[1];
  if (groups < 1 || groups > n || group_starts[1] != 1) {
    reject("Invalid group count or first start.");
  }
  for (g in 1:groups) {
    int group_start = group_starts[g];
    int group_end = g < groups ? group_starts[g + 1] - 1 : n;
    real lp_present;
    if (group_start < 1 || group_end < group_start || group_end > n) {
      reject("Invalid group block.");
    }
    lp_present = bernoulli_logit_lpmf(1 | Omega[group_start])
      + occupancy_single_lpmf(y[group_start:group_end] |
          mu[group_start:group_end], occ[group_start:group_end],
          n_unit[group_start:group_end], unit_starts[group_start:group_end],
          Q[group_start:group_end]);
    if (group_Q[group_start] == 1) {
      lp += lp_present;
    } else {
      lp += log_sum_exp(bernoulli_logit_lpmf(0 | Omega[group_start]), lp_present);
    }
  }
  return lp;
}
