// Each unit starts with its row count in unit_size; other rows contain -99.
// Any contiguous slice of complete units is a valid standalone input.
// Q repeats the unit detection indicator at every visit.
// Predictor inputs: mu is logit detection probability; occ is logit occupancy
// probability. Neither input is a probability or a log-probability.
real occupancy_single_lpmf(array[] int y, vector mu, vector occ,
                           array[] int unit_size,
                           array[] int Q) {
  int n = num_elements(y);
  int unit_start = 1;
  real lp = 0;

  if (num_elements(mu) != n || num_elements(occ) != n
      || num_elements(unit_size) != n
      || num_elements(Q) != n) {
    reject("occupancy_single: inputs must have the same event-row length");
  }
  if (n < 1) {
    reject("occupancy_single: at least one observed event is required");
  }
  while (unit_start <= n) {
    int unit_end = unit_start + unit_size[unit_start] - 1;
    real lp_present;
    if (unit_start < 1 || unit_end < unit_start || unit_end > n) {
      reject("occupancy_single: invalid event block");
    }
    lp_present = bernoulli_logit_lpmf(1 | occ[unit_start])
                 + bernoulli_logit_lpmf(y[unit_start:unit_end] | mu[unit_start:unit_end]);
    if (Q[unit_start] == 1) {
      lp += lp_present;
    } else {
      lp += log_sum_exp(bernoulli_logit_lpmf(0 | occ[unit_start]), lp_present);
    }
    unit_start = unit_end + 1;
  }
  return lp;
}
