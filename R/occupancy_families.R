#' Define the rep-varying occupancy family
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_single <- function() {
  brms::custom_family(
    "occupancy_single", dpars = c("mu", "occ"),
    links = c("identity", "identity"),
    type = "int", 
    # Unit count, event-block starts, and detection indicators, padded to event rows.
    vars = c("vint1", "vint2", "vint3"),
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

#' Define the two-level single-season occupancy family
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_twolevel_single <- function() {
  brms::custom_family(
    "occupancy_twolevel_single", dpars = c("mu", "occ", "Omega"),
    links = c("identity", "identity", "identity"),
    type = "int", 
    vars = paste0("vint", seq_len(6)),
    loop = FALSE)
}

#' Define the rep-varying augmented occupancy family
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_augmented <- function() {
  brms::custom_family(
    "occupancy_augmented", dpars = c("mu", "occ", "Omega"),
    links = c("identity", "identity", "identity"),
    type = "int", 
    vars = paste0("vint", seq_len(6)),
    loop = FALSE)
}

#' Define the colonization-extinction family
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_multi_colex <- function() {
  brms::custom_family(
    "occupancy_multi_colex", dpars = c("mu", "occ", "colo", "ex"),
    links = c("identity", "identity", "identity", "identity"),
    type = "int", 
    # Series count/start, unit count/start, and unit detection indicators.
    vars = paste0("vint", seq_len(5)),
    loop = FALSE)
}

#' Define the colonization-extinction equilibrium family
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_multi_colex_eq <- function() {
  brms::custom_family(
    "occupancy_multi_colex_eq", dpars = c("mu", "colo", "ex"),
    links = c("identity", "identity", "identity"),
    type = "int", 
    # Series count/start, unit count/start, and unit detection indicators.
    vars = paste0("vint", seq_len(5)),
    loop = FALSE)
}

#' Define the autologistic family
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_multi_autologistic <- function() {
  brms::custom_family(
    "occupancy_multi_autologistic", dpars = c("mu", "occ", "colo", "autologistic"),
    links = c("identity", "identity", "identity", "identity"),
    type = "int", 
    # Series count/start, unit count/start, and unit detection indicators.
    vars = paste0("vint", seq_len(5)),
    loop = FALSE)
}

#' Define the autologistic_eq family
#' @return a "customfamily" "brmsfamily" object from brms
#' @noRd
occupancy_multi_autologistic_eq <- function() {
  brms::custom_family(
    "occupancy_multi_autologistic_eq", dpars = c("mu", "colo", "autologistic"),
    links = c("identity", "identity", "identity"),
    type = "int", 
    # Series count/start, unit count/start, and unit detection indicators.
    vars = paste0("vint", seq_len(5)),
    loop = FALSE)
}
