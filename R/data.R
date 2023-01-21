#' Example flocker data
#'
#' Example data
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