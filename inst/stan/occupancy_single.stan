// One-level, single-season occupancy with event-varying detection.
// Each unit occupies one contiguous block of observed events, in visit order.
// y contains binary outcomes; mu and occ contain event-row linear predictors
// on the logit scale. n_unit[1] is the number of units; unit_start[1:n_unit[1]]
// contains their event-block starts. Both auxiliary arrays have event-row
// length, with unused entries padded with -99. Only occ at block starts is
// read. Units without observed events are not represented in this contract.
// Unit covariates are repeated at their actual values on every event row:
// brms therefore centers predictors and constructs bases at the visit level.
real occupancy_single_lpmf(array[] int y, vector mu, vector occ,
                           array[] int n_unit, array[] int unit_start) {
  int n = num_elements(y);
  int units;
  real lp = 0;

  if (num_elements(mu) != n || num_elements(occ) != n
      || num_elements(n_unit) != n || num_elements(unit_start) != n) {
    reject("occupancy_single: inputs must have the same event-row length");
  }
  if (n < 1) {
    reject("occupancy_single: at least one observed event is required");
  }
  units = n_unit[1];
  if (units < 1 || units > n || unit_start[1] != 1) {
    reject("occupancy_single: invalid unit count or first block start");
  }

  for (u in 1:units) {
    int start = unit_start[u];
    int end = u < units ? unit_start[u + 1] - 1 : n;
    real lp_present;
    if (start < 1 || end < start || end > n) {
      reject("occupancy_single: invalid event block for unit ", u);
    }
    lp_present = bernoulli_logit_lpmf(1 | occ[start])
                 + bernoulli_logit_lpmf(y[start:end] | mu[start:end]);
    if (sum(y[start:end]) > 0) {
      lp += lp_present;
    } else {
      lp += log_sum_exp(bernoulli_logit_lpmf(0 | occ[start]), lp_present);
    }
  }
  return lp;
}
