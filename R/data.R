#' Example flocker data
#'
#' Example data to be formatted for use in a single-season occupancy model
#'
#' @format ## `example_flocker_data`
#' A list with elements appropriate to pass to make_flocker_data:
#' \describe{
#'   \item{obs}{A detection/nondetection matrix. Rows are units and colums are
#'   visits}
#'   \item{unit_covs}{A data.frame of four unit covs. `uc1` and `uc2` are 
#'   numeric; `grp` and `species` are categorical. Every level of `species` is
#'   present for every value of `uc1` and `uc2`, as would be the case in a 
#'   typical multispecies model}
#'   \item{event_covs}{A named list of two matrices, each representing a numeric
#'   covariate that varies by visit but not by species}
#' }
#' @source data-raw/example_flocker_data.R
"example_flocker_data"


#' Example single-season flocker model
#'
#' A fitted single-season occupancy model from flocker
#'
#' @format ## `example_flocker_model_single`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_single"


#' Example multi-season flocker model (colex_explicit)
#'
#' A fitted multi-season occupancy model from flocker.
#' One species 200 point, six season, four rep, missing seasons, missing reps.
#' colex parameterization with explicit init.
#'
#' @format ## `example_flocker_model_multi`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_multi"


#' Example multi-season flocker model (auto_eq)
#'
#' A fitted multi-season occupancy model from flocker.
#' One species 200 point, six season, four rep, missing seasons, missing reps.
#' autologistic parameterization with equilibrium init.
#'
#' @format ## `example_flocker_model_multi_auto_eq`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_multi_auto_eq"

