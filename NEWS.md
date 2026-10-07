# flocker 1.1-0

* Added the `twolevel_single` model family for single-season occupancy models
  in which groups of closure units may themselves be structurally absent.
  Users specify the grouping factor with `level2_group`, may supply a separate
  table of group-level predictors with `level2_covs`, and model group-level
  occupancy with `f_meta`.
* Refactored the data-augmented model as a special case of the general
  two-level model, with shared data formatting, likelihood, fitting, and
  post-processing machinery.
* Breaking change for augmented models: every species supplied in the original
  `obs` array must have at least one detection. Never-observed pseudospecies
  must instead be added through `n_aug`.
* Breaking change for augmented models: `get_Z()` now returns a list containing
  `unit` states for species-site combinations and `level2` states for species.
  Generic two-level models use the same return structure. Joint samples from
  the two levels are mutually consistent.
* Two-level outputs from `get_Z()`, `log_lik_flocker()`, and `loo_flocker()` now
  retain user-supplied group names. Augmented models similarly retain supplied
  species names and assign stable names to augmented pseudospecies.
* Improved the numerical stability of two-level log-likelihood calculations by
  evaluating occupancy mixtures in log space.
* Fixed reconstruction of maximally ragged observation and covariate arrays,
  including multiseason data with a final season that is unobserved in every
  series, and made post-processing robust when exactly one posterior draw is
  requested.
* Fitted models and formatted data objects are now tagged with the version of
  `flocker` that created them. Two-level packing metadata is stored separately
  from the underlying `brmsfit` data.

# flocker 1.0-2
* fixed a bug that resulted in incorrect log likelihoods returned by log_lik_flocker for the data augmented model. Model fitting was not affected, and posteriors were correct.
* log-likelihood computations now enabled over new data.
* dramatically improved efficiency in certain post-processing functions for multiseason models
* CI runs properly
* backend hygiene improvements
* updated compatibility with upcoming changes to `loo_compare()` output
  structure in the `loo` package (> 2.9.0), which now returns a data frame
  instead of a matrix and includes additional diagnostic columns.

# flocker 1.0-1
* fixed a bug that was throwing an uninformative error when making predictions in multiseason models where history_condition is TRUE
* dramatic efficiency improvements to `get_positions` as applied to multiseason models.
* groundwork laid for enabling log-likelihood computations over new data.


# flocker 1.0-0

* Initial CRAN submission.
* Fixed bug in mixed predictive checking.
* Fixed bug that prevented post-processing with new data with different 
numbers of visits/seasons than original data.
* Substantially more thorough unit testing.
* Removed testing dependency on `cmdstanr` (not on CRAN)
* Refactored to avoid shipping large objects as package data
