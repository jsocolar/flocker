#' Define the occupancy model family with visit-varying covariates
#' Primarily for internal use in \code{flock()}.
#' @param max_visit the maximum number of repeat visits to a unit
#' @return a "customfamily" "brmsfamily" object from brms

occupancy_V <- function(max_visit) {
  brms::custom_family(
    "occupancy_V", dpars = c("mu", "occ"),
    links = c("identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_visit, Q, visit_index1...
    vars = c("vint1", "vint2", "vint3", paste0("vint", 3 + (1:max_visit))),
    loop = FALSE)
}


#' Define the occupancy model family with visit-constant covariates
#' Primarily for internal use in \code{flock()}.
#' @return a "customfamily" "brmsfamily" object from brms

occupancy_C <- function() {
  brms::custom_family(
    "occupancy_C", dpars = c("mu", "occ"),
    links = c("identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for trials
    vars = c("vint1[n]"),
    loop = TRUE)
}

