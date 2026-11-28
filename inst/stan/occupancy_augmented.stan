// Observed and augmented species are the groups in a two-level model.
// Predictor inputs are logits: mu for detection conditional on site occupancy,
// occ for site occupancy conditional on species inclusion, and Omega for species
// inclusion. These inputs are not probabilities or log-probabilities.
real occupancy_augmented_lpmf(array[] int y, vector mu, vector occ,
                              vector Omega, array[] int group_size,
                              array[] int unit_size, array[] int Q,
                              array[] int group_Q) {
  return occupancy_twolevel_single_lpmf(y | mu, occ, Omega, group_size,
                                        unit_size, Q, group_Q);
}
