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
