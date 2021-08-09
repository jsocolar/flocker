#' Create Stan code for the occupancy_lpmf function
#' Primarily for internal use in \code{flock()}.
#' @param links links for detection and occupancy respectively, passed to
#' \code{brms::custom_family(links)}
#' @return a "customfamily" "brmsfamily" object from brms
#' @export

occupancy_family <- function(links = c("logit", "logit")) {
  brms::custom_family(
    "occupancy", dpars = c("mu", "occ"),
    links = links,
    type = "real", 
    vars = c(".nsite", ".nvisit", ".Q", paste0(".visit_index", 1:.max_visit)),
    loop = FALSE)
}
