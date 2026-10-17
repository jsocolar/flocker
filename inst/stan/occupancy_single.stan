// One-level, single-season occupancy with event-varying detection.
// Each unit occupies one contiguous block of observed events, in visit order.
// y contains binary outcomes; mu and occ contain event-row linear predictors
// on the logit scale. n_unit[1] is the number of units; unit_starts[1:n_unit[1]]
// contains their event-block starts. Both auxiliary arrays have event-row
// length, with unused entries padded with -99. Only occ at block starts is
// read. Units without observed events are not represented in this contract.
// Unit covariates are repeated at their actual values on every event row:
// brms therefore centers predictors and constructs bases at the visit level.
real occupancy_single_lpmf(array[] int y, vector mu, vector occ,
                           array[] int n_unit, array[] int unit_starts) {
  int n = num_elements(y);
  int units;
  real lp = 0;

  if (num_elements(mu) != n || num_elements(occ) != n
      || num_elements(n_unit) != n || num_elements(unit_starts) != n) {
    reject("occupancy_single: inputs must have the same event-row length");
  }
  if (n < 1) {
    reject("occupancy_single: at least one observed event is required");
  }
  units = n_unit[1];
  if (units < 1 || units > n || unit_starts[1] != 1) {
    reject("occupancy_single: invalid unit count or first block start");
  }

  for (u in 1:units) {
    int unit_start = unit_starts[u];
    int unit_end = u < units ? unit_starts[u + 1] - 1 : n;
    real lp_present;
    if (unit_start < 1 || unit_end < unit_start || unit_end > n) {
      reject("occupancy_single: invalid event block for unit ", u);
    }
    lp_present = bernoulli_logit_lpmf(1 | occ[unit_start])
                 + bernoulli_logit_lpmf(y[unit_start:unit_end] | mu[unit_start:unit_end]);
    if (sum(y[unit_start:unit_end]) > 0) {
      lp += lp_present;
    } else {
      lp += log_sum_exp(bernoulli_logit_lpmf(0 | occ[unit_start]), lp_present);
    }
  }
  return lp;
}
