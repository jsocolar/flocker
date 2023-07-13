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

#' Example single-season rep-constant flocker model
#'
#' A fitted single-season rep-constant occupancy model from flocker
#'
#' @format ## `example_flocker_model_single_C`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_single_C"


#' Example data-augmented flocker model
#'
#' A fitted data-augmented occupancy model from flocker
#'
#' @format ## `example_flocker_model_aug`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_aug"


#' Example multi-season flocker model (colex explicit)
#'
#' A fitted multi-season occupancy model from flocker.
#' One species 200 point, six season, four rep, missing seasons, missing reps.
#' colex parameterization with explicit init.
#'
#' @format ## `example_flocker_model_multi_colex_ex`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_multi_colex_ex"


#' Example multi-season flocker model (colex equilibirum)
#'
#' A fitted multi-season occupancy model from flocker.
#' One species 200 point, six season, four rep, missing seasons, missing reps.
#' colex parameterization with equilibrium init.
#'
#' @format ## `example_flocker_model_multi_colex_eq`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_multi_colex_eq"


#' Example multi-season flocker model (autologistic explicit)
#'
#' A fitted multi-season occupancy model from flocker.
#' One species 200 point, six season, four rep, missing seasons, missing reps.
#' autologistic parameterization with explicit init.
#'
#' @format ## `example_flocker_model_multi_auto_ex`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_multi_auto_ex"

#' Example multi-season flocker model (autologistic equilibrium)
#'
#' A fitted multi-season occupancy model from flocker.
#' One species 200 point, six season, four rep, missing seasons, missing reps.
#' autologistic parameterization with equilibrium init.
#'
#' @format ## `example_flocker_model_multi_auto_eq`
#' A flocker_fit and brmsfit object
#' @source data-raw/example_flocker_model.R
"example_flocker_model_multi_auto_eq"

