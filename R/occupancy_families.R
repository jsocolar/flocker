#' Define the occupancy model family with visit-varying covariates
#' Primarily for internal use in \code{flock()}.
#' @param max_visit the maximum number of repeat visits to a site
#' @return a "customfamily" "brmsfamily" object from brms

occupancy_vv <- function(max_visit) {
  brms::custom_family(
    "occupancy_vv", dpars = c("mu", "occ"),
    links = c("logit", "logit"),
    type = "int", 
    # Integer aterms (vint) for nsite, nvisit, Q, visit_index1...
    vars = c("vint1", "vint2", "vint3", paste0("vint", 3 + (1:max_visit))),
    loop = FALSE)
}


#' Define the occupancy model family with visit-constant covariates
#' Primarily for internal use in \code{flock()}.
#' @param max_visit the maximum number of repeat visits to a site
#' @return a "customfamily" "brmsfamily" object from brms

occupancy_vc <- function() {
  brms::custom_family(
    "occupancy_vc", dpars = c("mu", "occ"),
    links = c("logit", "logit"),
    type = "int", 
    # Integer aterms (vint) for trials
    vars = c("vint1[n]"),
    loop = TRUE)
}