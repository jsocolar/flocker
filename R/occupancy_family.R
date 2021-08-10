#' Define the occupancy model family
#' Primarily for internal use in \code{flock()}.
#' @param max_visit the maximum number of repeat visits to a site
#' @return a "customfamily" "brmsfamily" object from brms
#' @export

occupancy_family <- function(max_visit) {
  brms::custom_family(
    "occupancy", dpars = c("mu", "occ"),
    links = c("logit", "logit"),
    type = "int", 
    # Integer aterms (vint) for nsite, nvisit, Q, visit_index1...
    vars = c("vint1", "vint2", "vint3", paste0("vint", 3 + (1:max_visit))),
    loop = FALSE)
}


