// Predictor inputs: mu is logit detection probability; occ is logit occupancy
// probability. Neither input is a probability or a log-probability.
real occupancy_single_C_lpmf(int y, real mu, real occ, int trials) {
  if (y == 0) {
    return log_sum_exp(bernoulli_logit_lpmf(0 | occ),
                       bernoulli_logit_lpmf(1 | occ) +
                         binomial_logit_lpmf(0 | trials, mu));
  } else {
    return bernoulli_logit_lpmf(1 | occ) +
      binomial_logit_lpmf(y | trials, mu);
  }
}
