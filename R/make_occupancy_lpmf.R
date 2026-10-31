##### single #####

#' Read the fixed-input single-season likelihood source
#' @return Character string of Stan code
#' @noRd
make_occupancy_single_lpmf <- function() {
  path <- system.file("stan", "occupancy_single.stan",
                      package = "flocker", mustWork = TRUE)
  paste(readLines(path), collapse = "\n")
}

##### single_C #####

#' Create Stan code for likelihood function occupancy_single_C_lpmf.
#' The purpose of defining this custom family, rather than using brms's zero-inflated 
#' binomial, is to ensure that the occupancy parameters are interpretable with
#' values of 1 in the marginalized state reflecting occupancy rather than non-
#' occupancy.
#' @return Character string of Stan code corresponding to occupancy_single_C_lpmf
#' @noRd
make_occupancy_single_C_lpmf <- function () {
  path <- system.file("stan", "occupancy_single_C.stan",
                      package = "flocker", mustWork = TRUE)
  paste(readLines(path), collapse = "\n")
}

##### twolevel_single #####

#' Create Stan code for likelihood function occupancy_twolevel_single_lpmf for 
#' a rep-varying two-level single-season model. 
#' @return Character string of Stan code corresponding to 
#'    occupancy_twolevel_single_lpmf
#' @noRd
make_occupancy_twolevel_single_lpmf <- function() {
  path <- system.file("stan", "occupancy_twolevel_single.stan",
                      package = "flocker", mustWork = TRUE)
  paste(make_occupancy_single_lpmf(), paste(readLines(path), collapse = "\n"), sep = "\n")
}

#' Create the augmented wrapper around the two-level likelihood
#' @return Stan functions for augmented occupancy
#' @noRd
make_occupancy_augmented_lpmf <- function() {
  path <- system.file("stan", "occupancy_augmented.stan",
                      package = "flocker", mustWork = TRUE)
  paste(make_occupancy_twolevel_single_lpmf(), paste(readLines(path), collapse = "\n"), sep = "\n")
}

# All four multiseason families share one stable forward likelihood source.
make_occupancy_multi_lpmf <- function(threaded = FALSE) {
  path <- system.file("stan", "occupancy_multi.stan", package = "flocker", mustWork = TRUE)
  code <- paste(readLines(path), collapse = "\n")
  if (threaded) {
    path <- system.file("stan", "occupancy_multi_colex_threaded.stan",
                        package = "flocker", mustWork = TRUE)
    code <- paste(code, paste(readLines(path), collapse = "\n"), sep = "\n")
  }
  code
}
