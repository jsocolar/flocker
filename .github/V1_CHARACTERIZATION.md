# 1.1 packing characterization

These tests establish a baseline for the migration described in V2_DESIGN.md;
that document remains authoritative. They do not prescribe the v2 row layout.

- `test-packing-characterization.R` checks alignment of observations and unique
  event/unit/group covariates through nontrivial two-level permutations, including
  a singleton group and singleton visits. It also covers aggregated single-season
  sufficient statistics, singleton event-map dimensions, named augmented species
  with site mappings, and internal versus trailing dynamic units.
- `test-utils.R` already covers deterministic ragged reconstruction across the
  event-varying single-season, two-level, augmented, and dynamic paths.
- `test-exact-likelihoods.R` supplies independent marginal-probability examples
  and brute-force dynamic histories. These R checks alone do not establish that
  generated Stan programs implement the same likelihood.
- `test-singleton-draws.R` and `test-trailing-multiseason.R` cover fitted output
  dimensions, group names, coherent two-level samples, and trailing-unit masks.

Keep the ecological expected values when replacing the reconstruction adapter.
Do not preserve the old variable-width index columns or representative-first
packing merely to satisfy characterization tests.

Intentional differences to account for during migration:

- 1.1 does not generally preserve input dimension names in reconstructed outputs;
  fitted outputs manufacture site/visit/timestep names. The new tests check
  existing group/species names but do not require loss of other names to persist.
- The aggregated `single_C` reconstruction repeats unit values over a rectangular
  visit grid, including unsampled visits. Its test records current behavior, not
  v2 availability semantics.
- Dynamic internal gaps retain unit predictors through dummy event rows; trailing
  units are omitted even with complete covariates. The tests assert reconstructed
  ecological values, not dummy-row values or their preprocessing effects.
- `newdata_checks` is not v2's explicit data role. Prediction-only unknown outcomes
  and trailing forecasts need new v2 tests rather than inherited 1.1 expectations.

No new MCMC fixtures are required by these additions. Existing fixture setup is
still eager; restructuring it into lazy per-family getters is separate work.
