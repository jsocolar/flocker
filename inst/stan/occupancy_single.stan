// One-level, single-season occupancy with event-varying detection.
// Each unit occupies one contiguous block of observed events, in visit order.
// y contains binary outcomes; mu and occ contain event-row linear predictors
// on the logit scale. n_event repeats the unit's positive event count on every
// row of its block. Only occ and n_event at the first row of each block are
// read. Units without observed events are not represented in this contract.
// Unit covariates are repeated at their actual values on every event row:
// brms therefore centers predictors and constructs bases at the visit level.
real occupancy_single_lpmf(array[] int y, vector mu, vector occ,
                           array[] int n_event) {
  int n = size(y);
  int start = 1;
  real lp = 0;

  if (num_elements(mu) != n || num_elements(occ) != n || size(n_event) != n) {
    reject("occupancy_single: inputs must have the same event-row length");
  }

  while (start <= n) {
    int count = n_event[start];
    int end;
    real lp_present;
    if (count < 1 || count > n - start + 1) {
      reject("occupancy_single: invalid event count at block start ", start);
    }
    end = start + count - 1;
    lp_present = bernoulli_logit_lpmf(1 | occ[start])
                 + bernoulli_logit_lpmf(y[start:end] | mu[start:end]);
    if (sum(y[start:end]) > 0) {
      lp += lp_present;
    } else {
      lp += log_sum_exp(bernoulli_logit_lpmf(0 | occ[start]), lp_present);
    }
    start = end + 1;
  }
  return lp;
}
