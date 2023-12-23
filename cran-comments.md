## R CMD check results

0 errors | 0 warnings | 1 note

* This is a second re-submission of a new release. Compared to the first
re-submission:
* "in R" removed from the package title
* `dontrun` replaced with `donttest` for long-running examples
* removed `dontrun` for fast-running examples
* working directory returned to original state at end of 
`inst/cache_vignettes.R`