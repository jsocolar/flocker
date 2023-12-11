# flocker 1.0-0

* Initial CRAN submission.
* Fixed bug in mixed predictive checking.
* Fixed bug that prevented post-processing with new data with different 
numbers of visits/seasons than original data.
* Substantially more thorough unit testing.
* Removed testing dependency on `cmdstanr` (not on CRAN)
* Refactored to avoid shipping large objects as package data