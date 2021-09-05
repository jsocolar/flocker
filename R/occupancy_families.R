#' Define the rep-varying occupancy family
#' Primarily for internal use in \code{flock()}.
#' @param max_rep the maximum number of repeat sampling events at a unit
#' @return a "customfamily" "brmsfamily" object from brms

occupancy_V <- function(max_rep) {
  brms::custom_family(
    "occupancy_V", dpars = c("mu", "occ"),
    links = c("identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_rep, Q, rep_index1...
    vars = c("vint1", "vint2", "vint3", paste0("vint", 3 + (1:max_rep))),
    loop = FALSE)
}


#' Define the rep-constant occupancy family
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

