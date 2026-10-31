// Slice whole series starts, never visits or seasons within a series.
// Predictor inputs mu, occ, colo, and ex are logits, as in the serial family.
real occupancy_multi_colex_partial_sum(array[] int starts_slice,
                                       int start, int end,
                                       array[] int y, vector mu, vector occ,
                                       vector colo, vector ex,
                                       array[] int series_starts, int n_series,
                                       array[] int n_unit, array[] int unit_starts,
                                       array[] int Q) {
  real lp = 0;
  for (i in 1:num_elements(starts_slice)) {
    int s = start + i - 1;
    int series_end = s < n_series ? series_starts[s + 1] - 1 : num_elements(y);
    lp += occupancy_multi_series_loglik(y, mu, occ, colo, ex,
                                        n_unit, unit_starts, Q,
                                        starts_slice[i], series_end);
  }
  return lp;
}

// Predictor inputs are logits. Compilation/runtime threading is enabled by
// brms::threading(force = TRUE); brms must not partition the packed visit rows.
real occupancy_multi_colex_threaded_lpmf(array[] int y, vector mu, vector occ,
                                         vector colo, vector ex,
                                         array[] int n_series, array[] int series_starts,
                                         array[] int n_unit, array[] int unit_starts,
                                         array[] int Q) {
  int n = num_elements(y);
  int series;
  if (n < 1 || num_elements(mu) != n || num_elements(occ) != n
      || num_elements(colo) != n || num_elements(ex) != n
      || num_elements(n_series) != n || num_elements(series_starts) != n
      || num_elements(n_unit) != n || num_elements(unit_starts) != n
      || num_elements(Q) != n) {
    reject("Multiseason inputs must have the same nonzero packed-row length.");
  }
  series = n_series[1];
  if (series < 1 || series > n || series_starts[1] != 1) {
    reject("Invalid multiseason block counts or first starts.");
  }
  // A minimum grainsize of one lets the dynamic scheduler combine series.
  return reduce_sum(occupancy_multi_colex_partial_sum,
                    series_starts[1:series], 1,
                    y, mu, occ, colo, ex, series_starts, series,
                    n_unit, unit_starts, Q);
}
