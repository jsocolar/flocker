#' Define the rep-varying occupancy family
#' @param max_rep the maximum number of repeat sampling events at a unit
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_single <- function(max_rep) {
  brms::custom_family(
    "occupancy_single", dpars = c("mu", "occ"),
    links = c("identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_rep, Q, rep_index1...
    vars = c(paste0("vint", seq(3 + max_rep))),
    loop = FALSE)
}

#' Define the rep-varying occupancy family with threading
#' @param max_rep the maximum number of repeat sampling events at a unit
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_single_threaded <- function (max_rep) {
  brms::custom_family(
    "occupancy_single_threaded", dpars = c("mu", "occ"),
    links = c("identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_rep, Q, rep_index1...
    vars = c(paste0("vint", seq(3 + max_rep))),
    loop = FALSE)
}

#' Define the rep-constant occupancy family
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_single_C <- function() {
  brms::custom_family(
    "occupancy_single_C", dpars = c("mu", "occ"),
    links = c("identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for trials
    vars = c("vint1[n]"),
    loop = TRUE)
}

#' Define the rep-varying augmented occupancy family
#' @param max_rep the maximum number of repeat sampling events at a unit
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_augmented <- function(max_rep) {
  brms::custom_family(
    "occupancy_augmented", dpars = c("mu", "occ", "Omega"),
    links = c("identity", "identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_rep, Q, ... rep_index1...
    vars = c(paste0("vint", seq(6 + max_rep))),
    loop = FALSE)
}

#' Define the colonization-extinction family
#' @param max_year the maximum number of seasons at a colex unit
#' @param max_rep the maximum number of repeat sampling events at a closure unit
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_multi_colex <- function(max_year, max_rep) {
  brms::custom_family(
    "occupancy_multi_colex", dpars = c("mu", "occ", "colo", "ex"),
    links = c("identity", "identity", "identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_rep, Q, rep_index1...
    vars = c(paste0("vint", seq(5 + max_year + max_rep))),
    loop = FALSE)
}

#' Define the colonization-extinction equilibrium family
#' @param max_year the maximum number of seasons at a colex unit
#' @param max_rep the maximum number of repeat sampling events at a closure unit
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_multi_colex_eq <- function(max_year, max_rep) {
  brms::custom_family(
    "occupancy_multi_colex_eq", dpars = c("mu", "colo", "ex"),
    links = c("identity", "identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_rep, Q, rep_index1...
    vars = c(paste0("vint", seq(5 + max_year + max_rep))),
    loop = FALSE)
}

#' Define the autologistic family
#' @param max_year the maximum number of seasons at a colex unit
#' @param max_rep the maximum number of repeat sampling events at a closure unit
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_multi_autologistic <- function(max_year, max_rep) {
  brms::custom_family(
    "occupancy_multi_autologistic", dpars = c("mu", "occ", "colo", "autologistic"),
    links = c("identity", "identity", "identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_rep, Q, rep_index1...
    vars = c(paste0("vint", seq(5 + max_year + max_rep))),
    loop = FALSE)
}

#' Define the autologistic_eq family
#' @param max_year the maximum number of seasons at a colex unit
#' @param max_rep the maximum number of repeat sampling events at a closure unit
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_multi_autologistic_eq <- function(max_year, max_rep) {
  brms::custom_family(
    "occupancy_multi_autologistic_eq", dpars = c("mu", "colo", "autologistic"),
    links = c("identity", "identity", "identity"),
    type = "int", 
    # Integer aterms (vint) for n_unit, n_rep, Q, rep_index1...
    vars = c(paste0("vint", seq(5 + max_year + max_rep))),
    loop = FALSE)
}
