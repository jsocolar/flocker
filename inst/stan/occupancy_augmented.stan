// Observed and augmented species are the groups in a two-level model.
real occupancy_augmented_lpmf(array[] int y, vector mu, vector occ,
                              vector Omega, array[] int n_group,
                              array[] int group_starts, array[] int n_unit,
                              array[] int unit_starts, array[] int Q,
                              array[] int group_Q) {
  return occupancy_twolevel_single_lpmf(y | mu, occ, Omega, n_group,
                                        group_starts, n_unit, unit_starts, Q, group_Q);
}
